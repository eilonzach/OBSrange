'''
FUNCTION locate.py

This is the main function of OBSrange.

Uses two-way travel time information and ship coordinates from an OBS survey to
invert for station location on the seafloor (Lat, Lon, Depth), turn-around time
(TAT), static correction to the sound velocity through the water column (dvp),
and the velocity of the ship in the radial direction of the survey circle (vr).

Josh R. & Zach E. & Stephen M. 4/23/18
'''
# Import modules and functions
import numpy as np
import numpy.linalg as LA
from funcs import pings, coord_txs, bootstrap, results, calc, ftest, plots, vel

def instruments(datafile, parameters):
  ################### Independent Parameter Initializations ####################
 
  print('\n Initializing independent parameters ...')
  vpw0     = parameters[0]    
  dvp0     = parameters[1]    
  tat0     = parameters[2]    
  N_bs     = parameters[3]    
  E_thresh = parameters[4]
  twtcorr  = parameters[5] 
  npts     = parameters[6]    
  dampx    = parameters[7]   
  dampy    = parameters[8]   
  dampz    = parameters[9]   
  damptat  = parameters[10]
  dampdvp  = parameters[11]
  eps      = parameters[12]    
  M        = parameters[13]      
  
  ######################### Load and Clean Input Data ##########################
  
  # Load data.
  print('\n Loading file ' + datafile.split('/')[-1] +' ...')
  data = pings.load(datafile)
  
  # Perform quality control on loaded data.
  print('\n Performing quality control ...')
  data, data_bad = pings.qc(data, vpw0, thresh=500)
  N_badpings = len(data_bad['twts'])
  print(' Number of pings removed: ' + str(N_badpings))
  
  ##################### Intermediate Variable Declarations #####################
  
  lon0 = data['lon_drop']   # Drop latitude                  -- decimal degrees
  lat0 = data['lat_drop']   # Drop longitude                 -- decimal degrees
  z0 = data['z_drop']       # Drop depth                     -- meters
  lats = data['lats']       # Latitudes of survey points     -- decimal degress
  lons = data['lons']       # Longitudes of survey points    -- decimal degrees
  ts = data['ts']           # Times survey points recieved   -- seconds
  twts = data['twts']       # Two-way travel-times of points -- msecs
  Nobs = len(data['twts'])  # Number of survey points        -- integer 
  sta = data['sta']         # Station name                   -- string
  
  # Convert coordinates of survey points to x-y and save in separate arrays. 
  xs = np.zeros(len(lons))
  ys = np.zeros(len(lats))
  zs = np.zeros(len(lons))
  xs, ys = coord_txs.latlon2xy(lat0, lon0, lats, lons)
  x0, y0 = coord_txs.latlon2xy(lat0, lon0, lat0, lon0)

  # Package coordinates for passing to other functions.
  drop_coords = [x0, y0, z0, lat0, lon0]
  ship_coords = [xs, ys, zs]
  coords = [drop_coords, ship_coords]
  
  # Calculate velocity vector of ship at each survey point. Apply smoothing.
  vs = vel.vector(xs, ys, zs, ts)
  vs = vel.smooth(vs, npts)

  ############################ Bootstrap Resampling ############################
  
  print('\n Performing bootstrap resampling ...')
  
  # Randomly resample model data. First columns are unpermutted input data.
  X, Y, Z, V, TWT, indxs = bootstrap.sampling(xs, ys, zs, vs, twts, N_bs)

  ################################## Inversion #################################
  
  print('\n Running bootstrap inversion ...')
  
  # Initialize starting model.
  m0_strt = np.array([x0, y0, z0, tat0, dvp0])
  
  # Initialize a results object to hold various results.
  R = results.resl(N_bs, Nobs)

  # Perform bootstrap inversion.
  R = bootstrap.inv(X, Y, Z, V, TWT, R, parameters, m0_strt, coords)

  # Unscramble randomly sampled data for plotting and evaluation.
  R.dtwts = np.mean(bootstrap.unscramble(R.dtwts, indxs), axis=1)
  R.twts = np.mean(bootstrap.unscramble(R.twts, indxs), axis=1)
  R.corrs = np.mean(bootstrap.unscramble(R.corrs, indxs), axis=1)
  R.vrs = np.mean(bootstrap.unscramble(R.vrs, indxs), axis=1)

  ################# F-test for uncertainty using a grid search #################

  print('\n Performing F-test ...')
  xg, yg, zg, Xg, Yg, Zg, P, mx, my, mz= ftest.test(R, coords, lat0, lon0, vpw0)
  
  ################################### Plots ####################################
  
  print('\n Creating Figures ...')
  
  # Histograms of model parameters.
  Nbins = 15
  fig1 = plots.model_histos(R, Nbins)
  
  # Survey map.
  fig2 = plots.survey_map(lat0, lon0, z0, lats, lons, zs, R, data_bad)
  
  # Model misfit histogram.
  fig3 = plots.misfit(R.E_rms, Nbins)
  
  # Model residuals at each site.
  fig4 = plots.residuals(lats, lons, xs, ys, vs, R, Nobs)
  
  # F-test plots
  fig5 = plots.ftest(xg, yg, zg, Xg, Yg, Zg, P, mx, my, mz, R)

  # Resolution and covariance plots
  #fig6 = plots.resolution_covariance()
  
  ######################### Package and Return Results #########################
  
  final_results = {'sta': sta,          # station
                   'xs': R.xs,          # sensor x-coordinate after each itr 
                   'ys': R.ys,          # sensor y-coordinate after each itr 
                   'zs': R.zs,          # sensor z-coordinate after each itr 
                   'lats': R.lats,      # sensor latitude after each itr
                   'lons': R.lons,      # sensor longitude after each itr
                   'tats': R.tats,      # turn-around times after each itr
                   'vpws': R.vpws,      # perturbations to vpw after each itr
                   'x_sta': R.xs,       # sta x-coordinate after each itr
                   'y_sta': R.ys,       # sta y-coordinate after each itr
                   'z_sta': R.zs,       # sta z-coordinate after each itr
                   'dzs': R.dzs,        # change to sensor depth after each itr
                   'E_rms': R.E_rms,    # rms after each itr
                   'twts': R.twts,      # two-way travel times
                   'dtwts': R.dtwts,    # travel-time residuals 
                   'corrs': R.corrs,    # travel-time corrections (if computed)
                   'xdrfts': R.dxdrfts, # drift in x at each survey point
                   'ydrfts': R.dydrfts, # drift in y at each survey point
                   'drifts': R.drifts,  # sensor drift distance
                   'azs': R.azs,        # sensor drift azimuth
                   'vrs': R.vrs,        # radial velocities of ship.
                   'Nbad': N_badpings,  # number of pings removed
                   'svy_lats': lats,    # lats of survey points
                   'svy_lons': lons,    # lons of survey points
                   'svy_xs': xs,        # x-coordinates of survey points
                   'svy_ys': ys,        # y-coordinates of survey points
                   'svy_zs': zs         # z-coordinates of survey points
                   }

  figs = [fig1, fig2, fig3, fig4, fig5]

  return final_results, figs

  #################################### FIN #####################################