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

# Histograms of model parameters
def model_histos(res, bins):
  # Set up pyplot fig and axes objects.
  fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(9,6))
  ax1 = axes[0,0]
  ax2 = axes[0,1]
  ax3 = axes[0,2]
  ax4 = axes[1,0]
  ax5 = axes[1,1]
  ax6 = axes[1,2]
  c = 'orangered'

  # Parameters to plot.
  p1 = res.lats
  p2 = res.lons
  p3 = res.zs
  p4 = res.tats
  p5 = res.vpws
  p6 = res.drifts

  # Plot histograms of each parameter
  ax1.hist(p1, bins, edgecolor='k', lw=1.0)
  ax1.axvline(np.mean(p1) - np.std(p1), 0, 1, color=c, ls='--', lw=2.5)
  ax1.axvline(np.mean(p1) + np.std(p1), 0, 1, color=c, ls='--', lw=2.5)
  ax1.axvline(np.mean(p1), 0, 1, color=c, lw=2.5)

  # Fix tick labels for lat.
  labels = ax1.get_xticks()
  labels = ['{: .5f}'.format(label) for label in labels]
  ax1.set_xticklabels(labels)
  ax1.locator_params(axis='x', nbins=3)

  ax1.set_ylim(0, 550)
  ax1.set_title('Latitude ($^\circ$)')
  ax1.locator_params(axis='x', nbins=3)
 
  ax2.hist(p2, bins, edgecolor='k', lw=1.0)
  ax2.axvline(np.mean(p2) - np.std(p2), 0, 1, color=c, ls='--', lw=2.5)
  ax2.axvline(np.mean(p2) + np.std(p2), 0, 1, color=c, ls='--', lw=2.5)
  ax2.axvline(np.mean(p2), 0, 1, color=c, lw=2.5)
  
  # Fix tick labels for lon.
  labels = ax2.get_xticks()
  labels = ['{: .5f}'.format(label) for label in labels]
  ax2.set_xticklabels(labels)

  ax2.set_ylim(0, 550)
  ax2.set_title('Longitude ($^\circ$)')
  ax2.set_yticks([])
  ax2.locator_params(axis='x', nbins=2)
  
  ax3.hist(p3, bins, edgecolor='k', lw=1.0)
  ax3.axvline(np.mean(p3) - np.std(p3), 0, 1, color=c, ls='--', lw=2.5)
  ax3.axvline(np.mean(p3) + np.std(p3), 0, 1, color=c, ls='--', lw=2.5)
  ax3.axvline(np.mean(p3), 0, 1, color=c, lw=2.5)
  
  ax3.set_ylim(0, 550)
  ax3.set_title('Depth (m)')
  ax3.set_yticks([])
  ax3.locator_params(axis='x', nbins=4)
  
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
  ax3.locator_params(axis='x', nbins=4)

  ax6.hist(p6, bins, edgecolor='k', lw=1.0)
  ax6.axvline(np.mean(p6) - np.std(p6), 0, 1, color=c, ls='--', lw=2.5)
  ax6.axvline(np.mean(p6) + np.std(p6), 0, 1, color=c, ls='--', lw=2.5)
  ax6.axvline(np.mean(p6), 0, 1, color=c, lw=2.5)
  
  ax6.set_ylim(0, 550)
  ax6.set_title('Drift (m)')
  ax6.set_yticks([])

  # Display fig.
  plt.tight_layout()
  #plt.show()

  # Return fig.
  return fig

def survey_map(lat0, lon0, z0, lats1, lons1, zs1, res, bad):
  # Set up pyplot fig and axes objects.
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8,4))
  ax1 = axes[0]
  ax2 = axes[1]

  # Parameters from results.
  lats2 = res.lats
  lons2 = res.lons
  zs2 = res.zs
  
  # Plot survey in map view.
  ax1.scatter(lons1, 
              lats1,
              s=60,
              marker='o', 
              c='g', 
              edgecolors='k',
              lw=1.0,
              label='Good pings')
  
  ax1.scatter(bad['lons'],
              bad['lats'],
              s=60,
              marker='o',
              c='r',
              edgecolors='k',
              lw=1.0,
              label='Bad pings')
  
  ax1.scatter(lon0,
              lat0,
              s=200,
              marker='v',
              c='cyan',
              edgecolors='k',
              lw=1.0,
              label='Drop point')
  
  ax1.scatter(np.mean(lons2),
              np.mean(lats2),
              s=250,
              marker='*',
              c='yellow',
              edgecolors='k',
              lw=1.0,
              label='Final location')
  
  ax1.legend(ncol=2, shadow=True, fancybox=True, bbox_to_anchor=(1.1, 1.1))

  # Fix tick labels for lon.
  xlabels = ax1.get_xticks()
  xlabels = ['{: .3f}'.format(label) for label in xlabels]
  ax1.set_xticklabels(xlabels)
  ax1.locator_params(axis='x', nbins=3)
  
  ax1.set_xlabel('Longitude ($^\circ$)')
  ax1.set_ylabel('Latitude ($^\circ$)')
  
  # Plot survey by depth. Start with ray-paths.
  for i in range(len(lats1)):
    ax2.plot([lats1[i], np.mean(lats2)*np.ones(len(lats1))[i]],
             [zs1[i], np.mean(zs2)*np.ones(len(lats1))[i]],
             c='gray',
             lw=0.2)

  ax2.scatter(lats1,
              zs1,
              s=60,
              marker='o',
              c='g',
              edgecolors='k',
              lw=1.0)

  ax2.scatter(bad['lats'],
              np.zeros(len(bad['lats'])),
              s=60,
              marker='o',
              c='r',
              edgecolors='k',
              lw=1.0)

  ax2.scatter(lat0,
              z0,
              s=200,
              marker='v',
              c='cyan',
              edgecolors='k',
              lw=1.0)

  ax2.scatter(np.mean(lats2),
              np.mean(zs2),
              s=250,
              marker='*',
              c='yellow',
              edgecolors='k',
              lw=1.0)

  # Fix tick labels for lat.
  xlabels = ax2.get_xticks()
  xlabels = ['{: .3f}'.format(label) for label in xlabels]
  ax2.set_xticklabels(xlabels)
  
  ax2.set_xlabel('Latitude ($^\circ$)')
  ax2.set_ylabel('Depth (m)')
  ax2.yaxis.tick_right()
  ax2.yaxis.set_label_position('right')
  ax2.set_xlim(min(lats1), max(lats1))

  # Display fig.
  plt.tight_layout()
  #plt.show()

  # Return fig.
  return fig

# Plot model misfit histogram
def misfit(data, bins):
  # Set up pyplot fig and axes objects.
  fig, ax = plt.subplots(nrows=1, ncols=1)
  c = 'orangered'

  # Get in units of ms
  data = data * 1000.
  
  # Plot histogram
  ax.hist(data, bins, edgecolor='k', lw=1.0)
  ax.axvline(np.mean(data) - np.std(data), 0, 1, color=c, ls='--', lw=2.5)
  ax.axvline(np.mean(data) + np.std(data), 0, 1, color=c, ls='--', lw=2.5)
  ax.axvline(np.mean(data), 0, 1, color=c, lw=2.5)

  ax.set_title('Misfit')
  ax.set_xlabel('RMS (ms)')
  
  #plt.show()

  # Return fig.
  return fig

# Plot various data/residuals at each site
def residuals(lats, lons, xs, ys, vs, res, N):
  # Set up pyplot fig and axes objects.
  fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8,6))
  ax1 = axes[0,0]
  ax2 = axes[0,1]
  ax3 = axes[1,0]
  ax4 = axes[1,1]

  corrs = res.corrs*1000
  tts = res.dtwts*1000
  azs = np.mean(res.az_locs, axis=0)

  # Ship velocities at each site.
  s1 = ax1.scatter(lons, lats, s=60, c=np.sum(vs, axis=1), edgecolors='k')
  ax1.plot(np.mean(lons), np.mean(lats))
  ax1.set_title('Ship velocity (m/s)')
  ax1.set_xlabel('Longitude ($^\circ$)')
  ax1.set_ylabel('Latitude ($^\circ$)')
  plt.colorbar(s1, ax=ax1)

  # Fix tick labels for lon.
  labels = ax1.get_xticks()
  labels = ['{: .2f}'.format(label) for label in labels]
  ax1.set_xticklabels(labels)
  ax1.locator_params(axis='x', nbins=3)
  
  # Travel-time corrections at each site.
  s2 = ax2.scatter(lons, lats, s=60, c=corrs, edgecolors='k')
  ax2.plot(np.mean(lons), np.mean(lats))
  ax2.set_title('Travel-time corrections (ms)')
  ax2.set_xlabel('Longitude ($^\circ$)')
  ax2.set_ylabel('Latitude ($^\circ$)')
  plt.colorbar(s2, ax=ax2)

  # Fix tick labels for lon.
  labels = ax2.get_xticks()
  labels = ['{: .2f}'.format(label) for label in labels]
  ax2.set_xticklabels(labels)
  ax2.locator_params(axis='x', nbins=3)
  
  # Travel-time residuals at each site.
  s3 = ax3.scatter(xs, ys, s=60, c=tts, edgecolors='k')
  ax3.plot(np.mean(xs), np.mean(ys))
  ax3.set_xlim(-4000, 4000)
  ax3.set_ylim(-4000, 4000)
  ax3.set_title('Travel-time residuals (ms)')
  ax3.set_xlabel('X (m)')
  ax3.set_ylabel('Y (m)')
  plt.colorbar(s3, ax=ax3)

  # Travel-time corrections by azimuth.
  s4 = ax4.scatter(azs, tts, s=60, c=corrs, edgecolors='k')
  ax4.set_title('Travel-time corrections (ms)')
  ax4.set_xlabel('Ship Azimuth ($^\circ$)')
  ax4.set_ylabel('Travel-time residuals (ms)')
  plt.colorbar(s4, ax=ax4)
  
  ax4.set_xlim(-5, 365)

  # Display.
  plt.tight_layout()
  #plt.show()

  # Return fig.
  return fig

# Error plots for F-test.
def ftest(xg, yg, zg, Xg, Yg, Zg, P, xmax, ymax, zmax, res):
  # Set up pyplot fig and axes objects.
  fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10,3))
  ax1 = axes[0]
  ax2 = axes[1]
  ax3 = axes[2]

  # Some parameters.
  ang = np.arange(0, 2*np.pi, 0.1)
  xs = res.xs
  ys = res.ys
  zs = res.zs

  # X-Y plane.
  c1 = ax1.contourf(xg, yg, P[:,:,zmax])
  ax1.contour(xg, yg, P[:,:,zmax])
  ax1.scatter(Xg[xmax,ymax,zmax],Yg[xmax,ymax,zmax],s=150,c='r',edgecolors='k')
  ax1.set_xlabel('X (m)')
  ax1.set_ylabel('Y (m)')
  ax1.set_title('X-Y')
  ax1.set_xlim(min(xg)-10,max(xg)+10)
  ax1.set_ylim(min(yg)-10,max(yg)+10)
  plt.colorbar(c1, ax=ax1)
  
  # X-Z plane.
  c2 = ax2.contourf(xg, zg, P[:,ymax,:])
  ax2.contour(xg, zg, P[:,ymax,:])
  ax2.scatter(Xg[xmax,ymax,zmax],Zg[xmax,ymax,zmax],s=150,c='r',edgecolors='k')
  ax2.set_xlabel('X (m)')
  ax2.set_ylabel('Z (m)')
  ax2.set_title('X-Z')
  ax2.set_xlim(min(xg)-10,max(xg)+10)
  ax2.set_ylim(min(zg)-10,max(zg)+10)
  plt.colorbar(c2, ax=ax2)
  
  # Y-Z plane.
  c3 = ax3.contourf(yg, zg, P[xmax,:,:])
  ax3.contour(yg, zg, P[xmax,:,:])
  ax3.scatter(Yg[xmax,ymax,zmax],Zg[xmax,ymax,zmax],s=150,c='r',edgecolors='k')
  ax3.set_xlabel('Y (m)')
  ax3.set_ylabel('Z (m)')
  ax3.set_title('Y-Z')
  ax3.set_xlim(min(yg)-10,max(yg)+10)
  ax3.set_ylim(min(zg)-10,max(zg)+10)
  plt.colorbar(c3, ax=ax3)
  
  # Display.
  plt.tight_layout()
  #plt.show()

  # Return fig.
  return fig