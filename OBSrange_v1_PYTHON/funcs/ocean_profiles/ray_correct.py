'''
Description
'''
# Import modules and functions
import os
import time
from funcs.ocean_profiles import get

def makegrid(lat, lon, z, stn, t, ssp_dir):
  from IPython.core.debugger import Tracer

  # Take absolute drop depth estimate.
  z = abs(z)

  # This infers whether z is in [m] or [km], should be [km]
  if z > 6:
    z = z/1000

  # Check if sound-speed profile directory exists. Make if necessary.
  if not os.path.exists(ssp_dir):
    os.mkdir(ssp_dir)

  # Declare name file is SUPPOSED to have.
  ssp_fle = ssp_dir + 'SSP_' + stn + '.txt'

  # If the current ssp file doesn't exist, then make one.
  if not os.path.exists(ssp_fle):

    # Compute month of the deployment from t.
    month = time.gmtime(t).tm_mon 

    # Write ssp file.
    Tracer()()
    get.lev(lat, lon, month + 1, ssp_fle, ssp_dir + 'ssp_09')