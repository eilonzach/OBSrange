'''
FUNCTION SET plots.py

A series of functions for plotting various results of the bootstrap inversion.
Not meant to be robust (i.e. reusable in general), just to keep main code
readable.

Stephen M. & Zach E. & Josh R. 4/23/18
'''
# Import modules and functions
import numpy as np
import matplotlib.pyplot as plt

from IPython.core.debugger import Tracer

# Histograms of model parameters
def model_histos(p1, p2, p3, p4, p5, p6, bins):
  # Set up pyplot fig and axes objects.
  fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(9,6))
  ax1 = axes[0,0]
  ax2 = axes[0,1]
  ax3 = axes[0,2]
  ax4 = axes[1,0]
  ax5 = axes[1,1]
  ax6 = axes[1,2]
  c = 'orangered'
  
  # Plot histograms of each parameter
  ax1.hist(p1, bins, edgecolor='k', lw=1.0)
  ax1.axvline(np.mean(p1) - np.std(p1), 0, 1, color=c, ls='--', lw=2.5)
  ax1.axvline(np.mean(p1) + np.std(p1), 0, 1, color=c, ls='--', lw=2.5)
  ax1.axvline(np.mean(p1), 0, 1, color=c, lw=2.5)
  
  ax1.set_ylim(0, 550)
  ax1.set_title('Latitude ($^\circ$)')
  ax1.locator_params(nbins=4)
 
  ax2.hist(p2, bins, edgecolor='k', lw=1.0)
  ax2.axvline(np.mean(p2) - np.std(p2), 0, 1, color=c, ls='--', lw=2.5)
  ax2.axvline(np.mean(p2) + np.std(p2), 0, 1, color=c, ls='--', lw=2.5)
  ax2.axvline(np.mean(p2), 0, 1, color=c, lw=2.5)
  
  ax2.set_ylim(0, 550)
  ax2.set_title('Longitude ($^\circ$)')
  ax2.set_yticks([])
  ax2.locator_params(nbins=4)
  
  ax3.hist(p3, bins, edgecolor='k', lw=1.0)
  ax3.axvline(np.mean(p3) - np.std(p3), 0, 1, color=c, ls='--', lw=2.5)
  ax3.axvline(np.mean(p3) + np.std(p3), 0, 1, color=c, ls='--', lw=2.5)
  ax3.axvline(np.mean(p3), 0, 1, color=c, lw=2.5)
  
  ax3.set_ylim(0, 550)
  ax3.set_title('Depth (m)')
  ax3.set_yticks([])
  
  p4 = p4*1000
  ax4.hist(p4, bins, edgecolor='k', lw=1.0)
  ax4.axvline(np.mean(p4) - np.std(p4), 0, 1, color=c, ls='--', lw=2.5)
  ax4.axvline(np.mean(p4) + np.std(p4), 0, 1, color=c, ls='--', lw=2.5)
  ax4.axvline(np.mean(p4), 0, 1, color=c, lw=2.5)
  
  ax4.set_ylim(0, 550)
  ax4.set_title('TAT (ms)')
  
  ax5.hist(p5, bins, edgecolor='k', lw=1.0)
  ax5.axvline(np.mean(p5) - np.std(p5), 0, 1, color=c, ls='--', lw=2.5)
  ax5.axvline(np.mean(p5) + np.std(p5), 0, 1, color=c, ls='--', lw=2.5)
  ax5.axvline(np.mean(p5), 0, 1, color=c, lw=2.5)
  
  ax5.set_ylim(0, 550)
  ax5.set_title('Water Velocity (m/s)')
  ax5.set_yticks([])

  ax6.hist(p6, bins, edgecolor='k', lw=1.0)
  ax6.axvline(np.mean(p6) - np.std(p6), 0, 1, color=c, ls='--', lw=2.5)
  ax6.axvline(np.mean(p6) + np.std(p6), 0, 1, color=c, ls='--', lw=2.5)
  ax6.axvline(np.mean(p6), 0, 1, color=c, lw=2.5)
  
  ax6.set_ylim(0, 550)
  ax6.set_title('Drift (m)')
  ax6.set_yticks([])
  
  # Attempt to fix ticks for ax1, ax2.
  plt.setp(ax1.get_xticklabels(), rotation=30, horizontalalignment='right')
  plt.setp(ax2.get_xticklabels(), rotation=30, horizontalalignment='right')

  # Display fig.
  plt.tight_layout()
  plt.show()

def survey_map(lat0, lon0, z0, lats1, lons1, zs1, lats2, lons2, zs2, bad):
  # Set up pyplot fig and axes objects.
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8,4))
  ax1 = axes[0]
  ax2 = axes[1]
  
  # Plot survey in map view.
  ax1.scatter(lons1, lats1, c='g', label='Good pings')
  ax1.scatter(bad['lons'], bad['lats'], c='r', label='Bad pings')
  ax1.scatter(lon0, lat0, c='b', label='Drop point')
  ax1.scatter(np.mean(lons2), np.mean(lats2), c='k', label='Final location')
  ax1.legend(ncol=2, shadow=True, fancybox=True, bbox_to_anchor=(1.1, 1.1))
  
  ax1.set_xlabel('Longitude ($^\circ$)')
  ax1.set_ylabel('Latitude ($^\circ$)')
  
  # Plot survey by depth.
  ax2.scatter(lats1, zs1, c='g')
  ax2.scatter(bad['lats'], np.zeros(len(bad['lats'])), c='r')
  ax2.scatter(lat0, z0, c='b')
  ax2.scatter(np.mean(lats2), np.mean(zs2), c='k')
  
  ax2.set_xlabel('Latitude ($^\circ$)')
  ax2.set_ylabel('Depth (m)')
  ax2.yaxis.tick_right()
  ax2.yaxis.set_label_position('right')

  # Display fig.
  #plt.figlegend(ncol=4, shadow=True, fancybox=True, bbox_to_anchor=(0.7, 1.2))
  plt.tight_layout()
  plt.show()

# Plot model misfit histogram
def misfit(data, bins):
  # Set up pyplot fig and axes objects.
  fig, ax = plt.subplots(nrows=1, ncols=1)
  c = 'orangered'
  
  # Plot histogram
  ax.hist(data, bins, edgecolor='k', lw=1.0)
  ax.axvline(np.mean(data) - np.std(data), 0, 1, color=c, ls='--', lw=2.5)
  ax.axvline(np.mean(data) + np.std(data), 0, 1, color=c, ls='--', lw=2.5)
  ax.axvline(np.mean(data), 0, 1, color=c, lw=2.5)

  ax.set_title('Misfit')
  ax.set_xlabel('RMS (ms)')
  
  plt.show()

# Plot various data/residuals at each site
def residuals(lats, lons, xs, ys, vs, corrs, tts, N):
  # Set up pyplot fig and axes objects.
  fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8,6))
  ax1 = axes[0,0]
  ax2 = axes[0,1]
  ax3 = axes[1,0]
  ax4 = axes[1,1]

  # Ship velocities at each site.
  s1 = ax1.scatter(lons, lats, c=np.sum(vs, axis=1))
  ax1.plot(np.mean(lons), np.mean(lats))
  ax1.set_title('Velocity (m/s)')
  ax1.set_xlabel('Longitude ($^\circ$)')
  ax1.set_ylabel('Latitude ($^\circ$)')
  plt.colorbar(s1, ax=ax1)
  
  # Travel-time corrections at each site.
  s2 = ax2.scatter(lons, lats, c=corrs*1000)
  ax2.plot(np.mean(lons), np.mean(lats))
  ax2.set_title('Travel-time corrections (ms)')
  ax2.set_xlabel('Longitude ($^\circ$)')
  ax2.set_ylabel('Latitude ($^\circ$)')
  plt.colorbar(s2, ax=ax2)
  
  # Travel-time residuals at each site.
  s3 = ax3.scatter(xs, ys, c=tts*1000)
  ax3.plot(np.mean(xs), np.mean(ys))
  ax3.set_xlim(-4000, 4000)
  ax3.set_ylim(-4000, 4000)
  ax3.set_title('Travel-time residuals (ms)')
  ax3.set_xlabel('X (m)')
  ax3.set_ylabel('Y (m)')
  plt.colorbar(s3, ax=ax3)

  # Travel-time residuals by observation number and velocity.
  s4 = ax4.scatter(range(N), tts*1000, c=np.sum(vs,  axis=1))
  ax4.set_title('Observation Residuals')
  ax4.set_xlabel('Observation #')
  ax4.set_ylabel('Travel-time residuals (ms)')
  plt.colorbar(s4, ax=ax4, label='Velocity (m/s)')
  
  # Display.
  plt.tight_layout()
  plt.show()

# Error plots for F-test.
def ftest(xg, yg, zg, Xg, Yg, Zg, P, xmax, ymax, zmax, xs, ys, zs):
  # Set up pyplot fig and axes objects.
  fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10,3))
  ax1 = axes[0]
  ax2 = axes[1]
  ax3 = axes[2]

  ang = np.arange(0, 2 * np.pi, 0.1)

  # X-Y plane
  c1 = ax1.contourf(xg, yg, P[:,:,zmax])
  ax1.contour(xg, yg, P[:,:,zmax])
  #ax1.plot(Xg[xmax, ymax, zmax], Yg[xmax, ymax, zmax])
  #ax1.plot(np.mean(xs), np.mean(ys), c='k')
  ax1.plot(np.mean(xs)+2*np.std(xs)*np.cos(ang), np.mean(ys)+2*np.std(ys)*np.sin(ang), c='r', ls='--')
  ax1.set_xlabel('X (m)')
  ax1.set_ylabel('Y (m)')
  ax1.set_title('X-Y')
  plt.colorbar(c1, ax=ax1)
  
  #plot(Xgrd(Ix_max,Iy_max,Iz_max),Ygrd(Ix_max,Iy_max,Iz_max),'ok','markerfacecolor',[1 0 0],'markersize',20,'linewidth',1)
  #plot(mean(x_sta),mean(y_sta),'ok','markerfacecolor',[0 0 0],'markersize',15,'linewidth',1)
  #plot(mean(x_sta)+2*std(x_sta)*cos(ang),mean(y_sta)+2*std(y_sta)*sin(ang),'--r','linewidth',1);

  # X-Z plane
  c2 = ax2.contourf(xg, zg, P[:, ymax, :])
  ax2.contour(xg, zg, P[:, ymax, :])
  #ax2.plot(Xg[xmax, ymax, zmax], Zg[xmax, ymax, zmax])
  #ax2.plot(np.mean(xs), np.mean(zs))
  ax2.plot(np.mean(xs)+2*np.std(xs)*np.cos(ang), np.mean(zs)+2*np.std(zs)*np.sin(ang), c='r', ls='--')
  ax2.set_xlabel('X (m)')
  ax2.set_ylabel('Z (m)')
  ax2.set_title('X-Z')
  plt.colorbar(c2, ax=ax2)
  
  # Y-Z plane
  c3 = ax3.contourf(yg, zg, P[xmax, :, :])
  ax3.contour(yg, zg, P[xmax, :, :])
  #ax3.plot(Yg[xmax, ymax, xmax], Zg[xmax, ymax, xmax])
  #ax3.plot(np.mean(ys), np.mean(zs))
  ax3.plot(np.mean(ys)+2*np.std(ys)*np.cos(ang), np.mean(zs)+2*np.std(zs)*np.sin(ang), c='r', ls='--')
  ax3.set_xlabel('Y (m)')
  ax3.set_ylabel('Z (m)')
  ax3.set_title('Y-Z')
  plt.colorbar(c3, ax=ax3)
  
  # Display
  plt.tight_layout()
  plt.show()