'''
FUNCTION ftest.py

Function to perform a formal F-test on two sets of residuals of fits to data
(res1, res2) that have respectively been fitted using parms1, parms2 (where 
these are integers equal to the number of parameters used in each fit).

IMPORTANTLY: 
  This will test whether the second fit (yielding res2) is statistically
  superior to the first fit, but using more parms.
    I.e.:          sum(abs(res2)) < sum(abs(res1)) 
    and            parms2 > parms1
  The degrees of freedom for each fit are therefore:
    1) N - parms1
    2) N - parms2
  where N is the number of data, or length(res*). The residuals are just equal 
  to dobs - dpred, so we square and sum these to get the chi^2 values (Ea)

Z. Eilon

J. Russell :
  This version takes the degrees of freedom as input (instead of model 
  parameters) for the case where v = N - M is not accurate.
  P > 1 : Model 1 actually fits the data better than Model 2 (model 1 smaller chi^2)
  P = 1 : Model 1 fits the data same as model 2 (same chi^2)
  P < 1 : Model 2 fits the data better than model 1 (model 2 smaller chi^2)
  P < 0.05 : Model 2 fits the data better model 1 with 95% confidence

S. Mosher: 
  Translation from MATLAB to Python
'''
# Import modules and functions
import numpy as np
from scipy.stats import f
from funcs import coord_txs, calc

def dof(res1, v1, res2, v2):
  # Calculate chi**2 sums
  Ea_1 = np.sum(res1**2)
  Ea_2 = np.sum(res2**2)
  
  Fobs = (Ea_1/v1)/(Ea_2/v2)
  P = 1 - ( f.cdf(Fobs, v1, v2) - f.cdf(1/Fobs, v1, v2) )

  return P

def test(res, xs, ys, zs, lat0, lon0, vpw):
  # Grab intermediate variables from results
  xM = np.mean(res.xs)
  yM = np.mean(res.ys)
  zM = np.mean(res.zs)
  tat = res.tats
  tatM = np.mean(tat)
  dvpM = np.mean(res.dvps)
  v_effM = np.mean(res.v_effs)
  
  # Set grid size and grids
  ngridpts = 40
  D = max(np.std(res.xs), np.std(res.ys), np.std(res.zs)) * 4
  Dx = D #std(x)*6; #10; # std(x)*2
  Dy = D #std(y)*6; #10; # std(y)*2
  Dz = D #std(z)*6; #10; # std(z)*2
  dx = 2*Dx/ngridpts #0.5; # 2*Dx/ngridpts;
  dy = 2*Dy/ngridpts #0.5; # 2*Dy/ngridpts;
  dz = 2*Dz/ngridpts #0.5; # 2*Dz/ngridpts;
  xg = np.linspace(xM - Dx, xM + Dx, ngridpts)
  yg = np.linspace(yM - Dy, yM + Dy, ngridpts)
  zg = np.linspace(zM - Dz, zM + Dz, ngridpts)
  lat_grid, lon_grid = coord_txs.xy2latlon(xg, yg, lat0, lon0)
  Nx = len(xg)
  Ny = len(yg)
  Nz = len(zg)
  
  # Bootstrap residuals
  twt_pre = calc.twt(xM, yM, zM, xs, ys, zs, vpw, dvpM, tatM)
  resids = res.twts - twt_pre
  
  # Determine the eigenvectors for z, vpw, and tat
  #from IPython.core.debugger import Tracer
  #Tracer()()
  #X = np.concatenate([zg, vpws, tat])
  #X = X.reshape(X.shape[0], 1)
  #V = LA.eigh(np.dot(X, X.T))[1]
  #
  #eigvec1 = V[:,0] # Closest to TAT axis
  #eigvec2 = V[:,1] # Closest to V_w axis
  #eigvec3 = V[:,2] # Closest to z axis
  #eig3_z = eigvec3[0]
  #eig3_vw = eigvec3[1]
  #eig3_tat = eigvec3[2]
  #print(eig3_z, eig3_vw, eig3_tat)
  
  # Grid search
  P = np.ndarray(shape=(Nx, Ny, Nz))
  E_gs = np.ndarray(shape=(Nx, Ny, Nz))
  Xg, Yg, Zg = np.meshgrid(xg, yg, zg)
  LONgrd, LATgrd, Zgrd = np.meshgrid(lon_grid, lat_grid, zg)
  
  # NESTED LOOPS SLOW IN PYTHON
  for i in range(Nx):
    for j in range(Ny):
      for k in range(Nz):
        # Apply scaling to vpw and tat to account for tradeoffs with Z
        #dz = zg[k] - np.mean(res.zs)
        #dvw = (eig3_vw/eig3_z) * dz # pertb to water velocity to account for dz
        #dtat = (eig3_tat/eig3_z) * dz # pertb to tat to account for dz
  
        # Grid search residual
        #twt_pre_gs = \
        #  calc.twt(xg[i], yg[j], zg[k], xs, ys, zs, vpw, dvpM+dvw, tatM + dtat)
        twt_pre_gs = calc.twt(xg[i], yg[j], zg[k], xs, ys, zs, vpw, dvpM, tatM)
        resid_gs = res.twts - twt_pre_gs
        
        # Calculate P statistic
        P[i,j,k] = dof(resid_gs, v_effM, resids, v_effM)  
        E_gs[i,j,k] = np.sqrt(np.mean(resid_gs**2))
  
  Pz_max= np.amax(P) 
  zmax = np.where(P == Pz_max)[0][0]

  Py_max = np.amax(P[:,:,zmax])
  ymax = np.where(P[:,:,zmax] == Py_max)[0][0]
  
  Px_max = np.amax(P[:,ymax,zmax])
  xmax = np.where(P[:,ymax,zmax] == Px_max)[0][0]

  return xg, yg, zg, Xg, Yg, Zg, P, xmax, ymax, zmax