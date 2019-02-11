'''
FUNCTION SET get.py

A set of routines to extract SSPs, temperature, salinity, or buoyancy profiles 
along a given point or points along a path (lat, lon). Uses lev_SSPs.py, a
function to load the variables and interpolate onto the desired location.

Gets the World Ocean Atlas files from directory ''
Modify this program to set the directory where files lev_ann.mat and  
lev_latlonZ.mat are located.

These files are preprocessed to save disk space, and have all shortened profiles
of the original data filled in by nearest neighbor data. All profiles extend to
5500 m depth, and there is world wide coverage (yes, including continents)

  See:  http://staff.washington.edu/dushaw/WOA/
  Data originally from:  http://www.nodc.noaa.gov/OC5/indprod.html

Data at points where the original WOA has data should be unchanged from the 
WOA. The annual mean World Ocean Atlas is used.
  
  lat, lon are a single point
 
  typeSSP = 0 indicates annual Levitus,
  typeSSP = 1:12 indicates Levitus for month 1:12
  typeSSP = 13:16 indicates Levitus for Winter, Spring, Summer, or Fall
  typeVAR = 1 indicates sound speed
            2 indicates temperature
            3 indicates salinity
            4 indicates buoyancy

UNITS:  sound speed in [m/s], depth in positive [meters].  

Returns P, a matrix of profiles: 33X(No. of points) and z, the standard 33
World Ocean Atlas depths.

Interpolation from World Ocean Atlas grid points to desired points by cubic
spline interpolation horizontally.

Writes out file "ofile"; if sound speed is called, this can be used in the RAY
or EIGENRAY programs.

The format of file templev.ssp is two columns: 
first,            "-1  range(m)', 
followed by   "depth(n) soundspeed (n)" for n=1:33
and repeated for points lat(m), lon(m)

UNITS: sound speed in [m/s], depth in positive meters.

Modified from original for simpler application (only sound speed) by 

Zach Eilon and Stephen Mosher 01/25/19
'''
# Import modules and functions
from IPython.core.debugger import Tracer


# Interpolate SSPs derived from Levitus database to get ssps at (lats, lons)
def lev_SSPs(type_SSP, lat, lon, dbase_dir):
  Tracer()()

  type_str = ['ann', 'jan', 'feb', 'mar', 'apr', 'may', 'jun', \
              'jul', 'aug', 'sep', 'oct', 'nov', 'dec']

  lev_fname = dbase_dir + '/lev_' + type_str[type_SSP]

  
# Main function.
def lev(lat, lon, type_SSP=0, ofile=False, database_dir='./ssp_09'):

  # Longitude convention conversion.
  if lon < 0:
    lon = lon + 360

  # Set standard data points and month for ssps.
  std_dpts = [0, 10, 20, 30, 50, 75, 100, 125, 200, 250, 300, 400, 500, 600, 
              700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000,
              2500, 3000, 3500, 4000, 4500, 5000, 5500]
  
  type_str = ['Annual', 'January', 'February', 'March', 'April', 'May', 'June',
             'July', 'August', 'September', 'October', 'November', 'December']

  type_str = type_str[type_SSP]

  Tracer()()
  # Obtains ssps from levitus database
  pre_d_SSPs = lev_SSPs(type_SSP, lon, lat, database_dir)