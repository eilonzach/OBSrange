'''
FUNCTION SET bootstrap.py

The first function gets indices for balanced resampling of data. It takes in a
numpy array (data) as well as the number of bootstrap iterations (N). The output
is an array of the same set of indices randomly permuted N times.

Second function re-indexes supplied data arrays according to previously permuted
indices. The overall reult of the first two functions is balanced random 
resampling, where each data point occurs N times in the bootstraping scheme.

The final function inverts the results of the previous two functions (i.e. it 
re-orders an array by its proper index order).

Josh R. & Zach E. & Stephen M. 4/23/18
'''
#Import modules and functions
import numpy as np

def get_bs_indxs(data, N):
  # Create an array of indices.
  Ndata = len(data)
  idxs = np.arange(0, Ndata)

  # Create N copies of indices.
  idxs = np.tile(idxs, N)

  # Randomly permute indices and place in a 2D array.
  rand_idxs = np.random.permutation(idxs)
  rand_idxs = np.reshape(rand_idxs, (Ndata, N))
  
  return rand_idxs

def sampling(x, y, z, v, t, N):
  # Get randomly sampled indices.
  idxs = get_bs_indxs(x, N)

  # Index parameter data by bootstrap indices then return.
  xmatbs = x[idxs]
  ymatbs = y[idxs]
  zmatbs = z[idxs]
  vmatbs = v[idxs]
  tmatbs = t[idxs]

  # Small adjustment to axes of vmatbs, given vmatbs has a different shape.
  vmatbs = np.swapaxes(vmatbs, 1, -1,)
 
  return xmatbs, ymatbs, zmatbs, vmatbs, tmatbs, idxs

def unscramble(data, idxs):
  # Initialze a new array for the unscrambled data.
  unscrambled_data = np.ndarray(shape=data.shape)
  Ninds = idxs.shape[0]
  
  # Re-index by proper index order.
  for index in range(Ninds):
    unscrambled_data[index] = data[index == idxs]
  
  return unscrambled_data