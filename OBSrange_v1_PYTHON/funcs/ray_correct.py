'''
Description
'''
# Import modules and functions
import os
import time
import numpy as np
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
  z = np.array([float(val) for val in z])/1000
  v_profile = np.array([z, ssp])

  Tracer()() 

  # Ray-shooting
  ps = np.arange(0.01, 1.501, 0.01)
  zs = z_sta + np.arange(-200, 220, 20)/1000
  zmax = max(zs)

  Np = len(ps)
  Nz = len(zs)

  Dr_ray = np.zeros(shape=(Np, Nz))
  Dx = np.zeros(shape=(Np, Nz))
  t_ray = np.zeros(shape=(Np, Nz))
  t_lin = np.zeros(shape=(Np, Nz))

  for i, ray in enumerate(ps):
    try:
      shootrays.shootrays(ray, v_profile, zmax)
    except Exception as e:
      Tracer()()
      