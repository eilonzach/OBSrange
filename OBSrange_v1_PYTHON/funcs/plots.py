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

  # Final lat, lon.
  fin_lat = str( round(np.mean(res.lats), 5) )
  fin_lon = str( round(np.mean(res.lons), 5) )

  # Parameters to plot.
  p1 = np.mean(res.ys) - res.ys
  p2 = np.mean(res.xs) - res.xs
  p3 = res.zs
  p4 = res.tats * 1000
  p5 = res.vpws
  p6 = res.drifts

  # Plot histograms of each parameter
  ax1.hist(p1, bins, edgecolor='k', lw=1.0)
  ax1.axvline(np.mean(p1) - np.std(p1), 0, 1, color=c, ls='--', lw=2.5)
  ax1.axvline(np.mean(p1) + np.std(p1), 0, 1, color=c, ls='--', lw=2.5)
  ax1.axvline(np.mean(p1), 0, 1, color=c, lw=2.5)
  ax1.set_title('(m) from ' + fin_lat + '$^\circ$N')
 
  ax2.hist(p2, bins, edgecolor='k', lw=1.0)
  ax2.axvline(np.mean(p2) - np.std(p2), 0, 1, color=c, ls='--', lw=2.5)
  ax2.axvline(np.mean(p2) + np.std(p2), 0, 1, color=c, ls='--', lw=2.5)
  ax2.axvline(np.mean(p2), 0, 1, color=c, lw=2.5)
  ax2.set_title('(m) from ' + fin_lon + '$^\circ$E')
  
  ax3.hist(p3, bins, edgecolor='k', lw=1.0)
  ax3.axvline(np.mean(p3) - np.std(p3), 0, 1, color=c, ls='--', lw=2.5)
  ax3.axvline(np.mean(p3) + np.std(p3), 0, 1, color=c, ls='--', lw=2.5)
  ax3.axvline(np.mean(p3), 0, 1, color=c, lw=2.5)
  ax3.set_title('Depth (m)')
  ax3.locator_params(axis='x', nbins=4)
  
  ax4.hist(p4, bins, edgecolor='k', lw=1.0)
  ax4.axvline(np.mean(p4) - np.std(p4), 0, 1, color=c, ls='--', lw=2.5)
  ax4.axvline(np.mean(p4) + np.std(p4), 0, 1, color=c, ls='--', lw=2.5)
  ax4.axvline(np.mean(p4), 0, 1, color=c, lw=2.5)
  ax4.set_title('TAT (ms)')
  
  ax5.hist(p5, bins, edgecolor='k', lw=1.0)
  ax5.axvline(np.mean(p5) - np.std(p5), 0, 1, color=c, ls='--', lw=2.5)
  ax5.axvline(np.mean(p5) + np.std(p5), 0, 1, color=c, ls='--', lw=2.5)
  ax5.axvline(np.mean(p5), 0, 1, color=c, lw=2.5)
  ax5.set_title('Water Velocity (m/s)')
  ax5.locator_params(axis='x', nbins=4)

  ax6.hist(p6, bins, edgecolor='k', lw=1.0)
  ax6.axvline(np.mean(p6) - np.std(p6), 0, 1, color=c, ls='--', lw=2.5)
  ax6.axvline(np.mean(p6) + np.std(p6), 0, 1, color=c, ls='--', lw=2.5)
  ax6.axvline(np.mean(p6), 0, 1, color=c, lw=2.5)
  ax6.set_title('Drift (m)')

  # Display fig.
  plt.tight_layout()
  #plt.show()

  # Return fig.
  return fig

def survey_map(lat0, lon0, z0, lats1, lons1, zs1, res, bad):
  # Set up pyplot fig and axes objects.
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9.6572, 6.75))
  ax1 = axes[0]
  ax2 = axes[1]

  # Parameters from results.
  lats2 = res.lats
  lons2 = res.lons
  zs2 = res.zs
  drop_diff = str( round(np.mean(zs2) - z0, 0) )
  tot_drift = str( round(np.mean(res.drifts), 1) )
  fin_azi = str(round(np.mean(res.azs), 1))
  
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
  
  xlabels = ax1.get_xticks()
  xlabels = ['{: .3f}'.format(label) for label in xlabels]
  ax1.set_xticklabels(xlabels)
  ax1.locator_params(axis='x', nbins=3)
  ax1.set_xlabel('Longitude ($^\circ$)')
  ax1.set_ylabel('Latitude ($^\circ$)')
  ax1.set_title('Drift: ' + tot_drift +' m     Azi: ' + fin_azi +'$^\circ$')
  ax1.legend(ncol=2, shadow=True, fancybox=True, loc=0)
  ax1.axis('square')
  
  # Plot survey by depth. Start with ray-paths.
  for i in range(len(lons1)):
    ax2.plot([lons1[i], np.mean(lons2)*np.ones(len(lons1))[i]],
             [zs1[i], np.mean(zs2)*np.ones(len(lons1))[i]],
             c='gray',
             lw=0.2)

  ax2.scatter(lons1,
              zs1,
              s=60,
              marker='o',
              c='g',
              edgecolors='k',
              lw=1.0)

  ax2.scatter(bad['lons'],
              np.zeros(len(bad['lons'])),
              s=60,
              marker='o',
              c='r',
              edgecolors='k',
              lw=1.0)

  ax2.scatter(lon0,
              z0,
              s=200,
              marker='v',
              c='cyan',
              edgecolors='k',
              lw=1.0)

  ax2.scatter(np.mean(lons2),
              np.mean(zs2),
              s=250,
              marker='*',
              c='yellow',
              edgecolors='k',
              lw=1.0)

  xlabels = ax2.get_xticks()
  xlabels = ['{: .3f}'.format(label) for label in xlabels]
  ax2.set_xticklabels(xlabels)
  ax2.yaxis.tick_right()
  ax2.yaxis.set_label_position('right')
  ax2.set_xlabel('Longitude ($^\circ$)')
  ax2.set_ylabel('Depth (m)')
  ax2.set_title('$Z_{final}$ - $Z_{drop}$ = ' + drop_diff + ' m')
  
  asp = np.diff(ax2.get_xlim())[0] / np.diff(ax2.get_ylim())[0]
  asp /= np.abs(np.diff(ax1.get_xlim())[0] / np.diff(ax1.get_ylim())[0])
  ax2.set_aspect(asp)

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

  corrs = res.corrs * 1000
  resids = res.dtwts * 1000
  azs = np.mean(res.az_locs, axis=1)
  sq_sum_vs = np.sqrt(np.sum(vs**2, axis=1))
  fin_lat = np.mean(res.lats)
  fin_lon = np.mean(res.lons)

  # Ship velocities at each site.
  s1 = ax1.scatter(lons, lats, s=60, c=sq_sum_vs, edgecolors='k')
  ax1.scatter(fin_lon, fin_lat, s=250, marker='*', c='yellow', edgecolors='k')
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
  s2 = ax2.scatter(lons, lats, s=60, c=corrs, edgecolors='k', cmap='RdBu_r')
  ax2.scatter(fin_lon, fin_lat, s=250, marker='*', c='yellow', edgecolors='k')
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
  s3 = ax3.scatter(lons, lats, s=60, c=resids, edgecolors='k', cmap='RdBu_r')
  ax3.scatter(fin_lon, fin_lat, s=250, marker='*', c='yellow', edgecolors='k')
  ax3.set_title('Travel-time residuals (ms)')
  ax3.set_xlabel('Longitude ($^\circ$)')
  ax3.set_ylabel('Latidute ($^\circ$)')
  plt.colorbar(s3, ax=ax3)

  # Fix tick labels for lon.
  labels = ax3.get_xticks()
  labels = ['{: .2f}'.format(label) for label in labels]
  ax3.set_xticklabels(labels)
  ax3.locator_params(axis='x', nbins=3)

  # Travel-time corrections by azimuth.
  s4 = ax4.scatter(azs, resids, s=60, c=corrs, edgecolors='k', cmap='RdBu_r')
  ax4.set_xlabel('Ship Azimuth ($^\circ$)')
  ax4.set_ylabel('Travel-time residuals (ms)')
  plt.colorbar(s4, ax=ax4, label='Travel-time corrections (ms)')
  
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
  x = np.mean(res.xs)
  y = np.mean(res.ys)
  z = np.mean(res.zs)
  xmed = str(round(np.median(xg), 1))
  ymed = str(round(np.median(yg), 1))
  zmed = str(round(np.median(zg), 1))

  # X-Y plane.
  c1 = ax1.contourf(Xg[:,:,zmax] - x,
                    Yg[:,:,zmax] - y,
                    P[:,:,zmax],
                    25)

  ax1.contour(Xg[:,:,zmax] - x,
              Yg[:,:,zmax] - y,
              P[:,:,zmax],
              levels=[0.05, 0.32],
              linewidths=2,
              colors='white')

  ax1.scatter(Xg[xmax,ymax,zmax] - x,
              Yg[xmax,ymax,zmax] - y,
              s=100,
              c='r',
              edgecolors='k',
              lw=1)
  
  ax1.set_xlabel('$\delta$X (m)')
  ax1.set_ylabel('$\delta$Y (m)')
  ax1.set_title('X-Y')
  plt.colorbar(c1, ax=ax1)

  xpos = ax1.get_xlim()[0] + 10
  ypos1 = ax1.get_ylim()[0] + 20
  ypos2 = ax1.get_ylim()[0] + 10
  ax1.text(x=xpos, y=ypos1, s='$\overline{x}$ = ' + xmed + ' m', color='white')
  ax1.text(x=xpos, y=ypos2, s='$\overline{y}$ = ' + ymed + ' m', color='white')

  # X-Z plane.
  c2 = ax2.contourf(Xg[ymax,:,:] - x,
                    Zg[ymax,:,:] - z,
                    P[ymax,:,:],
                    25)
  
  ax2.contour(Xg[ymax,:,:] - x,
              Zg[ymax,:,:] - z,
              P[ymax,:,:],
              levels=[0.05, 0.32],
              linewidths=2,
              colors='white')
  
  ax2.scatter(Xg[xmax,ymax,zmax] - x,
              Zg[xmax,ymax,zmax] - z,
              s=100,
              c='r',
              edgecolors='k',
              lw=1)

  ax2.set_xlabel('$\delta$X (m)')
  ax2.set_ylabel('$\delta$Z (m)')
  ax2.set_title('X-Z')
  plt.colorbar(c2, ax=ax2)
  
  xpos = ax2.get_xlim()[0] + 10
  ypos1 = ax2.get_ylim()[0] + 20
  ypos2 = ax2.get_ylim()[0] + 10
  ax2.text(x=xpos, y=ypos1, s='$\overline{x}$ = ' + xmed + ' m', color='white')
  ax2.text(x=xpos, y=ypos2, s='$\overline{z}$ = ' + zmed + ' m', color='white')
  
  # Y-Z plane.
  c3 = ax3.contourf(Yg[:,xmax,:] - y,
                    Zg[:,xmax,:] - z,
                    P[:,xmax,:],
                    25)

  ax3.contour(Yg[:,xmax,:] - y,
              Zg[:,xmax,:] - z,
              P[:,xmax,:],
              levels=[0.05, 0.32],
              linewidths=2,
              colors='white')
  
  ax3.scatter(Yg[xmax,ymax,zmax] - y,
              Zg[xmax,ymax,zmax] - z,
              s=100,
              c='r',
              edgecolors='k',
              lw=1)
  
  ax3.set_xlabel('$\delta$Y (m)')
  ax3.set_ylabel('$\delta$Z (m)')
  ax3.set_title('Y-Z')
  plt.colorbar(c3, ax=ax3)
  
  xpos = ax3.get_xlim()[0] + 10
  ypos1 = ax3.get_ylim()[0] + 20
  ypos2 = ax3.get_ylim()[0] + 10
  ax3.text(x=xpos, y=ypos1, s='$\overline{y}$ = ' + ymed + ' m', color='white')
  ax3.text(x=xpos, y=ypos2, s='$\overline{z}$ = ' + zmed + ' m', color='white')

  # Display.
  plt.tight_layout()
  #plt.show()

  # Return fig.
  return fig

# Plot model resolution and covariance
def resolution_covariance(R, M):
  # Set up pyplot fig and axes objects.
  fig, axes = plt.subplots(nrows=1, ncols=2)
  ax1 = axes[0]
  ax2 = axes[1]

  # Parameters
  resol = R.resol[:,:,0]
  cov = R.cov[:,:,0]
  model_cov = np.multiply( np.sign(cov), np.log10(np.abs(cov)))
  a = -np.max(np.max(np.log10(np.abs(cov))))
  b =  np.max(np.max(np.log10(np.abs(cov))))

  # Model Resolution
  m1 = ax1.imshow(resol, cmap='RdBu_r', vmin=-1, vmax=1, aspect='equal')
  
  ax1.set_title('Model Resolution')
  xlabels = ax1.get_xticks()
  xlabels = [str(int(label) + 1) for label in xlabels]
  ax1.set_xticklabels(xlabels)
  ylabels = ax1.get_yticks()
  ylabels = [str(int(label) + 1) for label in ylabels]
  ax1.set_yticklabels(ylabels)
  plt.colorbar(m1, fraction=0.046, pad=0.04, ax=ax1)

  # Add grid lines
  for x,y in zip(range(M), range(M)):
    ax1.axvline(x=x+0.5, ymin=0, ymax=1, color='k', lw=2.0)
    ax1.axhline(y=y+0.5, xmin=0, xmax=1, color='k', lw=2.0)

  # Model Covariance
  m2 = ax2.imshow(model_cov, cmap='RdBu_r', vmin=a, vmax=b, aspect='equal')
  
  ax2.set_title('Model Covariance')
  xlabels = ax2.get_xticks()
  xlabels = [str(int(label) + 1) for label in xlabels]
  ax2.set_xticklabels(xlabels)
  ylabels = ax2.get_yticks()
  ylabels = [str(int(label) + 1) for label in ylabels]
  ax2.set_yticklabels(ylabels)
  plt.colorbar(m2, fraction=0.046, pad=0.04, ax=ax2)

  # Add grid lines
  for x,y in zip(range(M), range(M)):
    ax2.axvline(x=x+0.5, ymin=0, ymax=1, color='k', lw=2.0)
    ax2.axhline(y=y+0.5, xmin=0, xmax=1, color='k', lw=2.0)

  # Display
  plt.tight_layout()
  #plt.show()

  # Return fig.
  return fig