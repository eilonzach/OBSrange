'''
Description
'''
# Import modules and functions
import os
import time
import pickle
import numpy as np
from scipy.stats import hmean
from scipy.interpolate import interp1d
from funcs import get, shootrays

def makegrid(lat, lon, z_sta, stn, t, stn_ssp_dir, ssp_dir):
  from IPython.core.debugger import Tracer
  
  # Take absolute drop depth estimate.
  z_sta = abs(z_sta)

  # This infers whether sation depth is in [m] or [km], should be [km]
  if z_sta > 6:
    z_sta = z_sta/1000

  # Name current station ssp is SUPPOSED to have.
  ssp_fname = stn_ssp_dir + 'SSP_' + stn + '.txt'

  # Check if ssp for current station provided by user. Otherwise create.
  if not os.path.exists(ssp_fname):

    # Compute month of the deployment from t.
    month = time.gmtime(t).tm_mon 

    # Get approptiate ssp/z profile from Levitus database using stn metadata.
    ssp, z = get.lev_based_ssp(lat, lon, month + 1, ssp_dir, ssp_fname)
  
  else:
    ssp, z = [], []
    f = open(ssp_fname, 'r')
    lines = f.readlines()
    for i, line in enumerate(lines):
      if i == 0:
        continue
      z.append(line.split()[0])
      ssp.append(line.split()[1])
    f.close()

  ssp = np.array([float(val) for val in ssp])
  z = np.array([float(val) for val in z])
  v_profile = np.array([z, ssp])/1000

  # Ray-shooting
  ps = np.arange(0.01, 1.501, 0.01)
  zs = z_sta + np.arange(-200, 220, 20)/1000
  zmax = max(zs)

  Np = len(ps)
  Nz = len(zs)

  Dr_ray = np.zeros(shape=(Np, Nz))
  Dx = np.zeros(shape=(Np, Nz))
  t_ray = np.zeros(shape=(Np, Nz))
  Dr_lin = np.zeros(shape=(Np, Nz))
  t_lin = np.zeros(shape=(Np, Nz))
  Dr_ray[:] = np.nan
  Dx[:] = np.nan
  t_ray[:] = np.nan
  Dr_lin[:] = np.nan
  t_lin[:] = np.nan

  for i, ray in enumerate(ps):
    try:
      rx, rz, Dr, rt, rv = shootrays.shootrays(ray, v_profile, zmax)
      # Get intersection of this ray at each depth
      for j, depth in enumerate(zs):
        f_Dr = interp1d(rz, Dr)
        Dr_ray[i,j] = f_Dr(depth)
        f_Dx = interp1d(rz, rx)
        Dx[i,j] = f_Dx(depth)
        f_t = interp1d(rz, rt)
        t_ray[i,j] = f_t(depth)
        Dr_lin[i,j] = np.sqrt(Dx[i,j]**2 + depth**2)
        t_lin[i,j] = Dr_lin[i,j] / hmean(rv[rz <= depth]) 
    except Exception:
      continue

  # difference in ray length( metres)
  dDr_m = 1000 * (Dr_ray - Dr_lin)
  
  # difference in ray travel time (ms) (TWO-way)
  dt_ray_ms = 2 * 1000 * (t_ray - t_lin) 
  # to be clear: this number is POSITIVE if ray is slower than straight line
  # approx. We should therefore SUBTRACT this from the observed travel times
  # to correct them from bent rays to straight lines. Having turned data from
  # rays to lines the code, which assumes lines, can then invert for position
  
  # Add vertical ray
  Dx = np.vstack([np.zeros(Nz), Dx])
  dDr_m = np.vstack([np.zeros(Nz), dDr_m])
  dt_ray_ms = np.vstack([np.zeros(Nz), dt_ray_ms])

  # Interpolate travel-time error onto regular mesh
  minimum = np.min(np.hstack([4, np.nanmax(Dx, axis=0)]))
  Dx_grid = np.arange(0, minimum + 0.005, 0.005) * 1000
  dT_grid = np.zeros(shape=(len(Dx_grid), Nz))
  dT_grid[:] = np.nan

  for i, depth in enumerate(zs):
    indx = ~np.isnan(dt_ray_ms[:,i])
    f_dT = interp1d(Dx[indx == True, i]*1000, dt_ray_ms[indx == True, i])
    dT_grid[:,i] = f_dT(Dx_grid)

  Tracer()()
  # Put results in a dictionary and write to disk
  dT_ray_v_line_grid = {'Dx_grid_m': Dx_grid,
                        'dz_grid_km': zs,
                        'dT_grid_ms': dT_grid}

  correction_file = stn_ssp_dir + 'cor_rayline_' + stn +'.pkl'
  f = open(correction_file, 'wb')
  pickle.dump(dT_ray_v_line_grid, f)
  f.close()

  return dT_ray_v_line_grid