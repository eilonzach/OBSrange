'''
CLASS results.py

A class to initialize various vectors and matrices, as well as an output
dictionary, all for holding results of the bootstrap inversion.

Stephen M. 04/23/18
'''
# Import modules and functions
import numpy as np

class res:
  '''
  Initialization
  '''
  def __init__(self, N=1, M=1):  
    # Vectors
    self.lats    = np.zeros(N) # sensor latitudes
    self.lons    = np.zeros(N) # sensor longitudes
    self.xs      = np.zeros(N) # sensor x positions
    self.ys      = np.zeros(N) # sensor y positions
    self.zs      = np.zeros(N) # sensor z positions
    self.dzs     = np.zeros(N) # changes in depth
    self.tats    = np.zeros(N) # sensor turn-around times
    self.vpws    = np.zeros(N) # sound speeds for water
    self.dvps    = np.zeros(N) # perturbations to water sound speeds
    self.drifts  = np.zeros(N) # sensor drift distance
    self.dxdrfts = np.zeros(N) # changes to drift in x-direction
    self.dydrfts = np.zeros(N) # changes to drift in y-direction
    self.azs     = np.zeros(N) # sensor drift azimuth
    self.E_rms   = np.zeros(N) # RMS error terms
    self.v_effs  = np.zeros(N) # effective degrees of freedom for F-test
    # Matrices
    self.dists   = np.ndarray(shape=(M, N)).T # ship distances from OBS
    self.az_locs = np.ndarray(shape=(M, N)).T # ship azimuthz from OBS
    self.twts    = np.ndarray(shape=(M, N))   # travel-times
    self.corrs   = np.ndarray(shape=(M, N))   # travel-time corrections
    self.dtwts   = np.ndarray(shape=(M, N))   # travel-time residuals
    self.vrs      = np.ndarray(shape=(M, N))  # radial velocities of ship.
    # Output dictionary
    self.models  = {}                         # holds each model update
  '''
  A function to update the results object.
  '''
  def update(self, i, m0, vpw, E, v, dt, tts, corr, vr, x0, y0, z0, xs, ys):
    self.xs[i]      = m0[0]
    self.ys[i]      = m0[1]
    self.zs[i]      = m0[2]
    self.tats[i]    = m0[3]
    self.dvps[i]    = m0[4]
    self.vpws[i]    = vpw + self.dvps[i]
    self.E_rms[i]   = E
    self.v_effs[i]  = v
    self.dtwts[:,i] = dt
    self.twts[:,i]  = tts
    self.corrs[:,i] = corr
    self.vrs[:,i]   = vr
    # Calculate OBS drift distance and azimuth.
    self.dxdrfts[i] = m0[0] - x0
    self.dydrfts[i] = m0[1] - y0
    self.dzs[i] = m0[2] - z0
    self.drifts[i] = np.sqrt(self.dxdrfts[i]**2 + self.dydrfts[i]**2)
    az = -(np.degrees( np.arctan2(self.dydrfts[i], self.dxdrfts[i]) )) + 90
    if az < 0:
      az += 360
    self.azs[i] = az
    # Calculate ship distance and azimuths from OBS.
    dxs = xs - m0[0]
    dys = ys - m0[1]
    self.dists[i,:] = np.sqrt(dxs**2 + dys**2)
    az_locs = - (np.degrees(np.arctan2(dys , dxs))) + 90
    az_locs[az_locs < 0 ] += 360
    self.az_locs[i,:] = az_locs