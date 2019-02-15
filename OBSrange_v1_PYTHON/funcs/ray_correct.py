'''
Description
'''
# Import modules and functions
import os
import time
from funcs import get

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

  # Ray-shooting
  Tracer()()