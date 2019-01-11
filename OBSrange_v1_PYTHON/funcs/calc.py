'''
FUNCTION SET calc.py

Functions for several intermediate calculations.

Josh R. & Zach E. & Stephen M. 4/23/18
'''
# Import modules and functions
import numpy as np

'''
A function to calculate two-way travel-times between a series of points
(xs, ys, zs) and a single reference point (x0, y0, z0) given a static sound 
speed for water (vpw), a sound speed perturbation (dvp), and a constant turn-
around-time for the sensor.
'''
def twt(x0, y0, z0, xs, ys, zs, vpw, dvp, tat):
  twt = 2 * np.sqrt((xs-x0)**2 + (ys-y0)**2 + (zs-z0)**2) / (vpw+dvp) + tat
  return twt

'''
A function to correct measured travel-times for ship's radial velocity.
'''
def tt_corr(x0, y0, z0, xs, ys, zs, vs, vp, dvp, tts):
  # Calculate ship's position unit vector.
  r = np.array([xs - x0, ys - y0, zs - z0])
  r_hat = r / np.sqrt( np.sum(r**2, axis=0) )

  # Correct travel-times for ship's radial velocity.
  vr = np.sum(vs.T * r_hat, axis=0)
  dr = vr * tts
  corrections = dr / (vp + dvp)
  corrected = tts + corrections 
  
  # Return. (+/-) if logging ship location at (receive/transmit) time.
  return corrected, corrections, vr

'''
A function to build the G matrix for the inversion.
'''
def G(x, y, z, dvp, x_ship, y_ship, z_ship, vpw, Nobs, M):
  # Initialize G matrix.
  G = np.ndarray(shape=(Nobs, M))
  
  # Anonymous function for calculating distances.
  D = lambda x0, y0, z0, x, y, z : np.sqrt((x0-x)**2 + (y0-y)**2 + (z0-z)**2)

  # First column: dti/dx
  G[:,0] = -(x_ship-x) * 2 / (vpw+dvp) / D(x_ship, y_ship, z_ship, x, y, z)
  # Second column: dti/dy
  G[:,1] = -(y_ship-y) * 2 / (vpw+dvp) / D(x_ship, y_ship, z_ship, x, y, z)
  # Third column: dti/dz
  G[:,2] = -(z_ship-z) * 2 / (vpw+dvp) / D(x_ship, y_ship, z_ship, x, y, z)
  # Fifth column: dti/ddvp
  G[:,3] = -2 * D(x_ship, y_ship, z_ship, x, y, z) / (vpw+dvp)**2

  return G