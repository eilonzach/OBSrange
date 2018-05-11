'''
SCRIPT locate_OBS_bootsrap.py

Uses two-way travel time information and ship coordinates from an OBS survey to
invert for station location on the seafloor (Lat, Lon, Depth), turn-around time
(TAT), static correction to the sound velocity through the water column (dvp),
and the velocity of the ship in the radial direction of the survey circle (vr0).

Josh R. & Zach E. & Stephen M. 4/23/18
'''
# Import modules and functions
import numpy as np
import numpy.linalg as LA
from funcs import pings, coord_txs, bootstrap, results, calc, ftest, plots, vel

from IPython.core.debugger import Tracer
#################### Independent Parameter Initializations #####################

print('\n Initializing independent parameters ...')
vpw0 = 1500.0     # Velocity of sound in water (m/s)
dvp0 = 0          # Inital perturbation to vpw (m/s)
tat0 = 0.013      # Sensor turn-around time (s)
N_bs = 2000       # Number of bootstrap iterations
E_thresh = 0.001  # RMS reduction threshold for inversion
twtcorr = False   # Option to apply a travel-time correction for ship velocity
npts = 5          # Npts in moving avg smoothing of ship vel. (1 = no smoothing)
dampx = 0         # Norm damping for each model parameter             
dampy = 0         #             .
dampz = 0         #             .    
damptat = 2e-1    #             .     
dampdvp = 5e-8    #             .      
eps = 1e-10       # Global norm damping for stabilization
M = 5             # Number of model parameters

########################### Load and Clean Input Data ##########################

# Specify location of survey data file and where to save results.
datafile = '../OrcaAcousticSurvey/EE04.txt'

# Load data.
print('\n Loading file ' + datafile.split('/')[-1] +' ...')
data = pings.load(datafile)

# Perform quality control on loaded data.
print('\n Performing quality control ...')
data, data_bad = pings.qc(data, vpw0, thresh=500)
N_badpings = len(data_bad['twts'])
print(' Number of pings removed: ' + str(N_badpings))

###################### Intermediate Variable Declarations ######################

lon0 = data['lon_drop']   # Drop latitude                   -- decimal degrees
lat0 = data['lat_drop']   # Drop longitude                  -- decimal degrees
z0 = data['z_drop']       # Drop depth                      -- meters
lats = data['lats']       # Latitudes of survey points      -- decimal degress
lons = data['lons']       # Longitudes of survey points     -- decimal degrees
ts = data['ts']           # Times survey points recieved    -- seconds
twts = data['twts']       # Two-way travel-times of points  -- msecs
Nobs = len(data['twts'])  # Number of survey points         -- integer 

# Convert coordinates of survey points to x-y and save in separate arrays. 
xs = np.zeros(len(lons))
ys = np.zeros(len(lats))
zs = np.zeros(len(lons))
xs, ys = coord_txs.latlon2xy(lat0, lon0, lats, lons)
x0, y0 = coord_txs.latlon2xy(lat0, lon0, lat0, lon0)

# Calculate velocity vector of ship at each survey point. Apply smoothing.
vs = vel.vector(xs, ys, zs, ts)
vs = vel.smooth(vs, npts)

################### Bootstrap Resampling and Initializations ###################

print('\n Preparing for bootsrap inversion ...')

# Initialize starting model.
m0_strt = np.array([x0, y0, z0, tat0, dvp0])

# Initialize a results object to hold various inversion results.
res = results.res(N_bs, Nobs)

# Randomly resample model data.
X, Y, Z, V, TWT, indxs = bootstrap.sampling(xs, ys, zs, vs, twts, N_bs)

# For now add a row of unpermutted data to each bootstrap matrix.
X = np.insert(X, 0, xs, axis=1)
Y = np.insert(Y, 0, ys, axis=1)
Z = np.insert(Z, 0, zs, axis=1)
V = np.insert(V, 0, vs, axis=2)
TWT = np.insert(TWT, 0, twts, axis=1)

############################## Inversion Procedure #############################

print('\n Running bootstrap inversion ...')

# Loop over a random sample of model parameters.
for i in range(N_bs):
  # Create a dictionary for the output models from each bootstrap iteration. 
  res.models[str(i)] = {}
  
  # Grab a set of randomly sampled model parameters for current bootstrap iter. 
  xbs = X[:,i]
  ybs = Y[:,i]
  zbs = Z[:,i]
  vbs = V[:,:,i]
  twtbs = TWT[:,i]

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
    res.models[str(i)][str(j)] = {'m0': m0, 'E': E, 'dtwt': dtwt}
    
    # Update inversion rms.
    if j > 0:
      dE = res.models[str(i)][str(j-1)]['E'] - res.models[str(i)][str(j)]['E']
    
    # Update model for next iteration and RMS counter.
    m0 = m1
    j += 1
  
  # Update results object with stabilized results of current bootstrap iteration
  res.update(i, m0, vpw0, E, v, dtwt, twtbs, corrections, vr, x0, y0, z0, xs, ys)
  
  # Convert stabilized result coords of current iteration back to lat lon.
  res.lats[i], res.lons[i] = \
                           coord_txs.xy2latlon(res.xs[i], res.ys[i], lat0, lon0)

# Unscramble randomly sampled data for plotting and evaluation.
res.dtwts = np.mean(bootstrap.unscramble(res.dtwts, indxs), axis=1)
res.twts = np.mean(bootstrap.unscramble(res.twts, indxs), axis=1)
res.corrs = np.mean(bootstrap.unscramble(res.corrs, indxs), axis=1)
res.vrs = np.mean(bootstrap.unscramble(res.vrs, indxs), axis=1)

################## F-test for uncertainty using a grid search ##################

print('\n Performing F-test ...')
xg, yg, zg, Xg, Yg, Zg, P, xmax, ymax, zmax = \
                                 ftest.test(res, xs, ys, zs, lat0, lon0, vpw0)

#################################### Plots #####################################

print('\n Displaying Figures ...')

# Histograms of model parameters.
Nbins = 15
plots.model_histos(res.lats,res.lons,res.zs,res.tats, res.vpws, res.drifts, Nbins)

# Survey map.
plots.survey_map(lat0,lon0,z0,lats,lons,zs, res.lats,res.lons,res.zs, data_bad)

# Model misfit histogram.
plots.misfit(res.E_rms * 1000, Nbins)

# Model residuals at each site.
plots.residuals(lats, lons, xs, ys, vs, res.corrs, res.dtwts, Nobs)

# F-test plots
plots.ftest(xg, yg, zg, Xg, Yg, Zg, P, xmax, ymax, zmax, res.xs, res.ys, res.zs)

##################################### FIN ######################################