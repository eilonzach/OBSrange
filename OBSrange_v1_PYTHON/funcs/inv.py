import numpy as np
import numpy.linalg as LA
from funcs import results, calc, coord_txs

from IPython.core.debugger import Tracer

def invert_till_stable(X, Y, Z, V, T, parameters, m0_strt, resl, xs, ys, x0, y0, z0, lat0, lon0):
  # Grab required parameters
  vpw0     = parameters[0]    
  dvp0     = parameters[1]    
  tat0     = parameters[2]    
  N_bs     = parameters[3]    
  E_thresh = parameters[4]
  twtcorr  = parameters[5]   
  dampx    = parameters[7]   
  dampy    = parameters[8]   
  dampz    = parameters[9]   
  damptat  = parameters[10]
  dampdvp  = parameters[11]
  eps      = parameters[12]   
  M        = parameters[13]
  Nobs     = X.shape[0]

  for i in range(N_bs):

    # Create a dictionary for the output models from each bootstrap iteration. 
    resl.models[str(i)] = {}
    
    # Grab a set of randomly sampled model parameters for current bootstrap it. 
    xbs = X[:,i]
    ybs = Y[:,i]
    zbs = Z[:,i]
    vbs = V[:,:,i]
    twtbs = T[:,i]
    
    # Reset starting model vector. Re-initialize RMS. 
    m0 = np.zeros( len(m0_strt) )
    m0[0] = m0_strt[0]
    m0[1] = m0_strt[1]
    m0[2] = m0_strt[2]
    m0[3] = m0_strt[3]
    m0[4] = m0_strt[4]
    dE = 1000
    
    # Count the number of iterations until dE becomes < E_thresh.
    j = 0 
    
    # Iterate models until RMS stabilizes.
    while dE > E_thresh:
      # Set current model parameters
      x = m0[0]
      y = m0[1]
      z = m0[2]
      tat = m0[3]
      dvp = m0[4]
      
      # Apply correction to two-way travel-times due to ship velocity.
      corrected, corrections, vr = \
                   calc.tt_corr(x, y, z, xbs, ybs, zbs, vbs, vpw0, dvp, twtbs)
      if twtcorr:
        twtbs = corrected
      
      # Build the G matrix.
      G = calc.G(x, y, z, dvp, xbs, ybs, zbs, vpw0, Nobs, M)
      
      # Set up norm damping for each parameter.
      H = np.eye(M) * np.diag([dampx, dampy, dampz, damptat, dampdvp])
      h = np.zeros(M)
      
      # Calculate predicted travel-times for this iteration and residuals.
      twt_pre = calc.twt(x, y, z, xbs, ybs, zbs, vpw0, dvp, tat)
      dtwt = twtbs - twt_pre
      # Calculate RMS error.
      E = np.sqrt( np.sum(dtwt**2) / Nobs)
      
      # Least squares solution. Model update.   
      F = np.concatenate([G, H])
      f = np.concatenate([dtwt, h])
      Finv = LA.solve(np.dot((F.T), F) + (eps * np.eye(M)), F.T)
      m1 = m0 + np.dot(Finv, f)
      
      # Calculate the effective degrees of freedom (for F-test)
      v = len(f) - np.trace(np.dot(F, Finv))
      # Record output of current iteration in m0, E, and dtwt.
      resl.models[str(i)][str(j)] = {'m0': m0, 'E': E, 'dtwt': dtwt}
      
      # Update inversion rms.
      if j > 0:
        dE = resl.models[str(i)][str(j-1)]['E']-resl.models[str(i)][str(j)]['E']
      
      # Update model for next iteration and RMS counter.
      m0 = m1
      j += 1
      
    # Update results object with stabilized results of current bootstrap it.
    resl.update(i,m0,vpw0,E,v,dtwt,twtbs,corrections,vr,x0,y0,z0,xs,ys)
    
    # Convert stabilized result coords of current iteration back to lat lon.
    resl.lats[i], resl.lons[i] = \
                     coord_txs.xy2latlon(resl.xs[i], resl.ys[i], lat0, lon0)