% Test how important the ellipsoid correction is...
clear

addpath('/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/OBSrange_v1_MATLAB/functions');

% Load data
% trudata = load('/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/paper/figs/figdata/PacificORCA_synthtest4/trudata_syn12.mat');
trudata = load('/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/paper/figs/figdata/PacificORCA_synthtest4_REVISION1_GPScorr/trudata_syn12_z5000m_fr10.mat');

% dlat = [7:95];
dlat = +1.7; % add to latitude to see how it changes things...
olon = trudata.drop_location(2); %-133.62; %0:360;
olat = trudata.drop_location(1)+dlat; %-7.54;
lon = trudata.survlon;
lat = trudata.survlat+dlat;
TAT = trudata.tat;
xOBS = trudata.obs_location_xyz(1)*1000;
yOBS = trudata.obs_location_xyz(2)*1000;
zOBS = -trudata.obs_location_xyz(3)*1000;
lonOBS = trudata.obs_location_laloz(2);
latOBS = trudata.obs_location_laloz(1)+dlat;
vp = trudata.vp_actual*1000;
azi = trudata.survaz;
dforward = trudata.dforward;
dstarboard = trudata.dstarboard;
v_surv_true = trudata.data.v_surv_true;

% Calculate GPS and transponder locations
survcog = atan2d(v_surv_true(:,2),v_surv_true(:,1));
[dx,dy] = GPS_transp_correction(dforward,dstarboard,survcog');
[ x_ship_GPS, y_ship_GPS, z_ship_GPS ] = lonlat2xy_nomap( olon, olat, lon, lat );
x_ship_TR = x_ship_GPS + dx;
y_ship_TR = y_ship_GPS + dy;
z_ship_TR = z_ship_GPS;

%% Distance correction
r_TR = sqrt( (xOBS-x_ship_TR).^2 + (yOBS-y_ship_TR).^2 );
r_GPS = sqrt( (xOBS-x_ship_GPS).^2 + (yOBS-y_ship_GPS).^2 );
dr = r_TR - r_GPS;

r3_TR = sqrt( (xOBS-x_ship_TR).^2 + (yOBS-y_ship_TR).^2 + (zOBS-z_ship_TR).^2 );
r3_GPS = sqrt( (xOBS-x_ship_GPS).^2 + (yOBS-y_ship_GPS).^2 + (zOBS-z_ship_GPS).^2);
dr3 = r3_TR - r3_GPS;

figure(1); clf;
set(gcf,'Position',[237.0000   84.0000  683.0000  615.0000]);
ax1 = subplot(2,1,1);
scatter(ax1,lon,lat,80,dr3,'o','filled','MarkerEdgeColor',[0 0 0]); hold on;
plot(lonOBS,latOBS,'sk','MarkerFaceColor',[0.5 0.5 0.5],'markersize',13);
text(max(lon)-.011,max(lat)-.002,{'mean : ',[num2str(mean(dr3),'%.1f'),' m']},'fontsize',14,'fontweight','bold');
axis equal;
ylim(ax1,[min(lat)-(max(abs(lat))-min(abs(lat)))*0.05 , max(lat)+(max(abs(lat))-min(abs(lat)))*0.05]);
xlim(ax1,[min(lon)-(max(abs(lon))-min(abs(lon)))*0.05 , max(lon)+(max(abs(lon))-min(abs(lon)))*0.05]);
set(gca,'fontsize',15,'linewidth',1.5,'box','on');
ylabel(ax1,'Latitude (\circ)','fontsize',15);
xlabel(ax1,'Longitude (\circ)','fontsize',15);
title(ax1,'Transponder-GPS offset relative to OBS','FontWeight','bold','fontsize',18);
pos = get(gca,'Position');
cb = colorbar(ax1);
set(gca,'Position',pos);
ylabel(cb,'r_{transp.} - r_{GPS} (m)','fontsize',15);
colormap(ax1,redbluecmap)
caxis(ax1,max(abs(dr3))*[-1 1]);

%% Residual travel-time correction 
TWT_TR = 2 .* r3_TR ./ vp + TAT;
TWT_GPS = 2 .* r3_GPS ./ vp + TAT;
dTWT = TWT_TR-TWT_GPS;
azi(azi<1) = azi(azi<1)+360;

subplot(2,1,2);
% scatter(azi,dTWT_ellipsoid*1000,130,dr_ellipsoid2,'o','filled','MarkerEdgeColor',[0 0 0]);
plot(azi,dTWT*1000,'ok','MarkerFaceColor',[0.5 0.5 0.5],'markersize',14);
text(300,-2,['mean : ',num2str(mean(dTWT*1000),'%.1f'),' ms'],'fontsize',14,'fontweight','bold');
set(gca,'fontsize',15,'linewidth',1.5,'box','on');
xlabel('Ship Azimuth (\circ)','fontsize',15);
ylabel('TWT_{transp.} - TWT_{GPS} (ms)','fontsize',15);
title(['Perturbation to travel times'],'FontWeight','bold','fontsize',18);
% ylim([-6 2]);
xlim([50 370]);
pos = get(gca,'Position');
% cb2 = colorbar;
% ylabel(cb2,'r_{ellip} - r_{sphere} (m)','fontsize',15);
% caxis([-4 4]);
% colormap(redbluecmap);

%%

save2pdf('GPS_corr.pdf',1,100);