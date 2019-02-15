'''
FUNCTION shootrays.py 

1D ray tracing function all units in km or km/s

Z. Eilon 01/2019 (modified from code by Brandon Schmandt)
Translated to Python and further modified by S. Mosher 01/2019
'''
# Import modules and functions
import numpy as np

from IPython.core.debugger import Tracer


# Given two arrays, this function finds the closest value in the smaller array
# to each value in the larger array.
def find_ilay(r, R):
  
  N = len(r)
  ilay = np.zeros(N, dtype=int)
  
  for i, depth in enumerate(r):
    idx = np.abs(R - depth).argmin()
    ilay[i] = idx

  return ilay

# Linear interpolation function.
def linterp(X, Y, XI):
  
  YI = np.zeros(len(XI))

  ilay = find_ilay(XI, X)
  
  # Linearly interpolate.
  YI = (XI - X[ilay]) * (Y[ilay + 1] - Y[ilay]) / (X[ilay + 1] - X[ilay]) + Y[ilay]

  # Sort out coincident elements.
  olap = np.intersect1d(X, XI)
  for i, val in enumerate(olap):
    YI[XI == olap[i]] = np.mean(Y[X == olap[i]])

  return YI



def shootrays(p, v_profile, zmax, dr=0.001, vdz=0.001):

  # Anonymous functions
  eta = lambda u, p: np.sqrt(u**2 + p**2)
  gradv_dist = lambda b, u1, u2, p: (eta(u1,p)/(b*u1*p)) - (eta(u2,p)/(b*u2*p))
  constv_dist = lambda u, dz, p: (p * dz) / eta(u, p)

  zz1 = np.arange(0, zmax, vdz)
  zz2 = np.arange(vdz, zmax+vdz, vdz)
  

  v1  = linterp(v_profile[0], v_profile[1], zz1)
  v2  = linterp(v_profile[0], v_profile[1], zz2)
  
  Tracer()()

  # assume flat Earth (no radius terms)
  u1 = 1 / v1
  u2 = 1 / v2 
  
  dv = v2 - v1
  
  dz = zz2 - zz1
  b = dv / dz
  const_indx = (b == 0)
  
  X = np.zeros(len(v1))
  X[const_indx == True] = constv_dist( u1[const_indx == True], dz[const_indx == True], p )
  X[const_indx != True] = gradv_dist(b[const_indx != True], u1[const_indx != True],u2[const_indx != True],p)
  X = np.hstack([0, X])
  X[np.imag(X) > 0] = np.nan
  X = X[X != np.isnan(X)]
  Xd = np.cumsum(X)
  rayz = np.hstack([zz1[1], zz2]) 
  #rayz = rayz(1:length(X))
  #dz = dz(1:length(X)-1)
  
  # calc Dr
  Dr1 = np.sqrt(X**2 + dz**2)
  assert(np.imag(Dr1).any() == False)
  Dr1 = np.hstack([0, Dr1])
  Dr2 = np.cumsum(Dr1)
  #Dr = np.unique( np.hstack( [np.linspace(0, max(Dr2) + dr, dr), max(Dr2)])).T
  #rx = interp1(Dr2,Xd,Dr)
  #rz = interp1(Dr2,rayz,Dr)
  #
  #% calculate travel time (JBR)
  #dDr = Dr(2:length(Dr))-Dr(1:length(Dr)-1)
  #dDr = [0; dDr]
  #v = [v1(1; v2)]
  #rv = interp1(Dr2,v,Dr)
  #dt = dDr./rv
  #rt = cumsum(dt)

  #return 