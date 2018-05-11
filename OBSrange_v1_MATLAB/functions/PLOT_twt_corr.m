%% map the ship velocities and the inversion misfits (and twt corrections) 
f3 = figure(3); clf;
set(gcf,'pos',[60 105 1100 800]);

% axes
ax1 = axes('pos',[0.1 0.56 0.37 0.38]);  hold(ax1,'on');
ax2 = axes('pos',[0.55 0.56 0.37 0.38]); hold(ax2,'on');
ax3 = axes('pos',[0.1 0.08 0.37 0.38]);  hold(ax3,'on');
ax4 = axes('pos',[0.55 0.08 0.37 0.38]);  hold(ax4,'on');

% limits
latlim = [min(lats_ship) max(lats_ship)]+[-1 1]*.2/111; % 200m outside ship circle
lonlim = [min(lons_ship) max(lons_ship)]+[-1 1]*.2/111; % 200m outside ship circle
londist = distance_km(mean(lats_ship),lonlim(1),mean(lats_ship),lonlim(2));




% Plot Velocities
scatter(ax1,lons_ship,lats_ship,100,sqrt(sum(v_ship.^2,1)),'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on;
plot(ax1,mean(lon_sta),mean(lat_sta),'pk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1);

title(ax1,'Velocity (m/s)','fontsize',18);
xlabel(ax1,'Longitude','fontsize',18);
ylabel(ax1,'Latitude','fontsize',18);

set(ax1,'fontsize',16,'linewidth',2,'box','on',...
        'xlim',lonlim,'ylim',latlim); grid(ax1,'on');
colorbar('peer',ax1);
% axis(ax1,'equal');

% Plot twt correction at each site
scatter(ax2,lons_ship,lats_ship,100,dtwtcorr_bs*1000,'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on;
plot(ax2,mean(lon_sta),mean(lat_sta),'pk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1);

title(ax2,'TWT corrections (ms)','fontsize',18);
xlabel(ax2,'Longitude','fontsize',18);
ylabel(ax2,'Latitude','fontsize',18);

set(ax2,'fontsize',16,'linewidth',2,'box','on',...
        'xlim',lonlim,'ylim',latlim); grid(ax2,'on');
colorbar('peer',ax2);
% axis(ax2,'equal');


% Plot residuals at each site
scatter(ax3,lons_ship,lats_ship,100,dtwt_bs*1000,'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on;
plot(ax3,mean(lon_sta),mean(lat_sta),'pk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1);

title(ax3,'Residuals (ms)','fontsize',18);
xlabel(ax3,'Longitude','fontsize',18);
ylabel(ax3,'Latitude','fontsize',18);

set(ax3,'fontsize',16,'linewidth',2,'box','on',...
        'xlim',lonlim,'ylim',latlim); grid(ax3,'on');
colorbar('peer',ax3);
% axis(ax3,'equal');


% Plot residuals as a function of azimuth
obs_num = 1:length(dtwt_bs);
scatter(ax4,azi_locs{1},dtwt_bs*1000,100,dtwtcorr_bs*1000,'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on

title(ax4,'Observation Residuals (ms)','fontsize',18);
xlabel(ax4,'Ship Azimuth (deg)','fontsize',18);
ylabel(ax4,'Residuals (ms)','fontsize',18);
set(ax4,'fontsize',16,'linewidth',2,'box','on'); grid(ax4,'on');
cb = colorbar('peer',ax4);
axis(ax1,'equal');
ylabel(cb,'TWT corrections (ms)');




% subplot(2,2,4);
% obs_num = 1:length(dtwt_bs);
% scatter(obs_num,dtwt_bs*1000,100,sqrt(sum(v_ship.^2,1)),'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on;
% grid on; box on; set(gca,'fontsize',16,'linewidth',2);
% title('Observation Residuals (ms)','fontsize',18);
% xlabel('Observation #','fontsize',18);
% ylabel('Residuals (ms)','fontsize',18);
% cb = colorbar;
% ylabel(cb,'Velocity (m/s)');