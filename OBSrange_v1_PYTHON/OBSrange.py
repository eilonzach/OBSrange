'''
SCRIPT OBSrange.py

This script sets parameters and loops through several survey files to locate
instruments on the seafloor via a bootstrap inversion procedure. The main 
function that gets called is locate.instruments(). The guts of the code are 
contained in that function and several functions within. The results, written in
an output directory, consist of .pkl and .txt files of the inversion results, as
well as several figures. 

Stephen M. & Zach E. & Josh R. 4/23/18
'''
# Import modules and functions
import os
import pickle
from funcs import fetch, locate, output

############################# Inversion Parameters #############################

vpw = 1500.0      # Assumed velocity of sound in water (m/s)
dvp = 0           # Assumed perturbation to vpw (m/s)
tat = 0.013       # Assumed sensor turn-around time (s)
N_bs = 500        # Number of bootstrap iterations
E_thresh = 1e-5   # RMS reduction threshold for inversion
twtcorr = False   # Option to apply a travel-time correction for ship velocity
npts = 1          # Npts in moving avg smoothing of ship vel. (1 = no smoothing)
dampx = 0         # Norm damping for each model parameter             
dampy = 0         #             .
dampz = 0         #             .    
damptat = 2e-1    #             .     
dampdvp = 5e-8    #             .      
eps = 1e-10       # Global norm damping for stabilization
M = 5             # Number of model parameters being solved for by inversion
QC = True         # Option to perform quality control on survey data points
res_thresh = 500  # Threshold (in ms) beyond which survey data points are tossed
bounds = [0.005, 
          0.025]  # Acceptable bounds for turn-around time solutions.

parameters = [vpw, dvp, tat, N_bs, E_thresh, twtcorr, npts, dampx, dampy, \
              dampz, damptat, dampdvp, eps, M, QC, res_thresh, bounds]

###################### Run Inversion for Each Survey File ######################

# Specify location of survey files.
survey_fles = fetch.data_paths('./survey_files/', matchkey='*.txt')

# Specify head output directory. Make output directories if necessary.
output_dir = './output/'
out_pkls = output_dir+'data_pkls/'
out_plts = output_dir+'plots/'
out_txts = output_dir+'data_txts/'
if not os.path.exists(output_dir):
  os.mkdir(output_dir)
  os.mkdir(out_pkls)
  os.mkdir(out_plts)
  os.mkdir(out_txts)

# Perform locations for each survey site then save output.
for survey_fle in survey_fles:
  results, figs = locate.instruments(survey_fle, parameters)
  output.out(results, figs, out_pkls, out_plts, out_txts)
  
##################################### FIN ######################################