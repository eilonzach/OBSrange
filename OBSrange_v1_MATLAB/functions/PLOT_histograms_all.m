%% plot the histograms of model parameters
% histogram bins
Nbins = 15;

% set up figure
f100 = figure(100); clf;
set(gcf,'position',[3 1 1208 697],'color','white');

% set up axes
ax1 = axes('pos',[0.08 0.58 0.25 0.35]);
ax2 = axes('pos',[0.38 0.58 0.25 0.35]);
ax3 = axes('pos',[0.68 0.58 0.25 0.35]);
ax4 = axes('pos',[0.08 0.12 0.25 0.35]);
ax5 = axes('pos',[0.38 0.12 0.25 0.35]);
ax6 = axes('pos',[0.68 0.12 0.25 0.35]);

% make the histograms
[ Ncount_lat, cent_lat ] = plot_hist(ax1,lat_sta,Nbins);
title(ax1,'\textbf{Latitude}','fontsize',18,'interpreter','latex');

[ Ncount_lon, cent_lon ] = plot_hist(ax2,lon_sta,Nbins);
title(ax2,'\textbf{Longitude}','fontsize',18,'interpreter','latex');

[ Ncount_z, cent_z ] = plot_hist(ax3,z_sta,Nbins);
title(ax3,'\textbf{Depth (m)}','fontsize',18,'interpreter','latex');

[ Ncount_TAT, cent_TAT ] = plot_hist(ax4,TAT*1000,Nbins);
title(ax4,'\textbf{TAT (ms)}','fontsize',18,'interpreter','latex');

[ Ncount_vw, cent_vw ] = plot_hist(ax5,V_w,Nbins);
title(ax5,'\textbf{Water Velocity (m/s)}','fontsize',18,'interpreter','latex');

[ Ncount_drift, cent_drift ] = plot_hist(ax6,drift,Nbins);
title(ax6,'\textbf{Drift (m)}','fontsize',18,'interpreter','latex');

%% ticks for the histograms - in m about central value
[cent_x,cent_y] = lonlat2xy(mean(cent_lon),mean(cent_lat),cent_lon,cent_lat);
% ticks every two m
tick_x = [2*floor(min(cent_x)/2):2:2*ceil(max(cent_x)/2)];
tick_y = [2*floor(min(cent_y)/2):2:2*ceil(max(cent_y)/2)];
[~,tick_lat] = xy2lonlat(mean(cent_lon),mean(cent_lat),mean(cent_lon),tick_x);
[tick_lon,~] = xy2lonlat(mean(cent_lon),mean(cent_lat),tick_y,mean(cent_lat));
set(ax1,'xtick',tick_lat,'xticklabel',tick_x,'xlim',[min(tick_lat),max(tick_lat)])
set(ax2,'xtick',tick_lon,'xticklabel',tick_y,'xlim',[min(tick_lon),max(tick_lon)])
xlabel(ax1,sprintf('m from %.5f$^{\\circ}$',mean(cent_lat)),'interpreter','latex');
xlabel(ax2,sprintf('m from %.5f$^{\\circ}$',mean(cent_lon)),'interpreter','latex');

%% drift directions
ax7 = polaraxes('pos',[0.621 0.35 0.16 0.16]); hold on
driftaz = atan2d(x_sta,y_sta);
plot_drift_polar( ax7,driftaz );
