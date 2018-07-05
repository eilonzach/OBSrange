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
  where N is the number of data, or len(res*). The residuals are just equal to
  dobs - dpred, so we square and sum these to get the chi^2 values (Ea)

Z. Eilon

J. Russell:
  This version takes the degrees of freedom as input (instead of model 
  parameters) for the case where v = N - M is not accurate.
  P > 1 : Model 1 actually fits the data better than Model 2 (model 1 smaller
          chi^2)
  P = 1 : Model 1 fits the data same as model 2 (same chi^2)
  P < 1 : Model 2 fits the data better than model 1 (model 2 smaller chi^2)
  P < 0.05 : Model 2 fits the data better model 1 with 95% confidence

S. Mosher: 
  Translation from MATLAB to Python
'''
# Import modules and functions
import numpy as np
import scipy.linalg as LA
from scipy.stats import f
from funcs import coord_txs, calc

def dof(res1, v1, res2, v2):
  # Calculate chi**2 sums.
  Ea_1 = np.sum(res1**2)
  Ea_2 = np.sum(res2**2)
  
  Fobs = (Ea_1/v1)/(Ea_2/v2)
  P = 1 - ( f.cdf(Fobs, v1, v2) - f.cdf(1/Fobs, v1, v2) )
  
  return P

def test(R, coords, lat0, lon0, vpw):
  # Grab coords.
  xs = coords[1][0]
  ys = coords[1][1]
  zs = coords[1][2]

  # Grab intermediate variables from results.
  xM = np.mean(R.xs)
  yM = np.mean(R.ys)
  zM = np.mean(R.zs)
  tatM = np.mean(R.tats)
  dvpM = np.mean(R.dvps)
  v_effM = np.mean(R.v_effs)
  #v_effM = np.mean(R.v_effs[0:-1])
  
  # Set grid size and grids.
  ngridpts = 40
  D = max(np.std(R.xs), np.std(R.ys), np.std(R.zs)) * 4
  Dx = D 
  Dy = D 
  Dz = D 
  dx = 2*Dx/ngridpts 
  dy = 2*Dy/ngridpts 
  dz = 2*Dz/ngridpts
  xg = np.linspace(xM - Dx, xM + Dx, ngridpts)
  yg = np.linspace(yM - Dy, yM + Dy, ngridpts)
  zg = np.linspace(zM - Dz, zM + Dz, ngridpts)
  assert len(xg) == len(yg) == len(zg)

  lat_grid, lon_grid = coord_txs.xy2latlon(xg, yg, lat0, lon0)
  Nx = len(xg)
  Ny = len(yg)
  Nz = len(zg)

  # Bootstrap residuals.
  twt_pre = calc.twt(xM, yM, zM, xs, ys, zs, vpw, dvpM, tatM)
  resid_bs = R.twts - twt_pre
  
  # Determine the eigenvectors for z, vpw, and tat.
  X = np.vstack([R.zs, R.vpws, R.tats]).T
  V = LA.eigh(np.dot(X.T, X))[1]
  
  eigvec1 = V[:,0] # Closest to TAT axis
  eigvec2 = V[:,1] # Closest to V_w axis
  eigvec3 = V[:,2] # Closest to z axis
  eig3_z = eigvec3[0]
  eig3_vw = eigvec3[1]
  eig3_tat = eigvec3[2]
  
  # Grid search
  P = np.ndarray(shape=(Nx, Ny, Nz))
  E_gs = np.ndarray(shape=(Nx, Ny, Nz))
  Xg, Yg, Zg = np.meshgrid(xg, yg, zg)
  LONgrd, LATgrd, Zgrd = np.meshgrid(lon_grid, lat_grid, zg)
  
  for i in range(Nx):
    for j in range(Ny):
      for k in range(Nz):
        # Apply scaling to vpw and tat to account for tradeoffs with z.
        dz = Zg[i,j,k] - zM
        # Perturbation to water velocity to account for dz.
        dvw = (eig3_vw/eig3_z) * dz
        # Perturbation to tat to account for dz. 
        dtat = (eig3_tat/eig3_z) * dz 
        
        # Grid search residual.
        twt_pre_gs = calc.twt(Xg[i,j,k], Yg[i,j,k], Zg[i,j,k], xs, ys, zs,
                              vpw, dvpM + dvw, tatM + dtat)

        resid_gs = R.twts - twt_pre_gs
        
        # Calculate P statistic.
        P[i,j,k] = dof(resid_gs, v_effM, resid_bs, v_effM)  
        E_gs[i,j,k] = np.sqrt( np.sum(resid_gs**2) / len(resid_gs))

  xmax, ymax, zmax = np.where( P == np.amax(P) )

  return xg, yg, zg, Xg, Yg, Zg, P, xmax[0], ymax[0], zmax[0]