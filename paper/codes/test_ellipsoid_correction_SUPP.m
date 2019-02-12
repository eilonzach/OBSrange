% Test how important the ellipsoid correction is...
clear

addpath('/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/OBSrange_v1_MATLAB/functions');

% Load data
trudata = load('/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/paper/figs/figdata/PacificORCA_synthtest4/trudata_syn12.mat');

% dlat = [7:95];
dlat = -89;%+1.7; % add to latitude to see how it changes things...
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

[ x_elli, y_elli , z_elli] = lonlat2xy_nomap( olon, olat, lon, lat );
[ x_sph, y_sph , z_sph] = lonlat2xy_nomap_sphere( olon, olat, lon, lat );

% [ lon_sph, lat_sph ] = xy2lonlat_nomap( olon, olat, x_sph, y_sph );
% dr_ellipsoid2 = distance(lon,lat,lon_sph,lat_sph);
% dr_ellipsoid2 = deg2km(dr_ellipsoid2)*1000;
%% Distance correction
dr_ellipsoid = sqrt((y_sph-y_elli).^2+(x_sph-x_elli).^2+(z_sph-z_elli).^2);

figure(1); clf;
set(gcf,'Position',[237.0000   84.0000  683.0000  615.0000]);
ax1 = subplot(2,1,1);
scatter(ax1,lon,lat,80,dr_ellipsoid,'o','filled','MarkerEdgeColor',[0 0 0]); hold on;
plot(lonOBS,latOBS,'sk','MarkerFaceColor',[0.5 0.5 0.5],'markersize',13);
axis equal;
ylim(ax1,[min(lat)-(max(abs(lat))-min(abs(lat)))*0.05 , max(lat)+(max(abs(lat))-min(abs(lat)))*0.05]);
xlim(ax1,[min(lon)-(max(abs(lon))-min(abs(lon)))*0.05 , max(lon)+(max(abs(lon))-min(abs(lon)))*0.05]);
set(gca,'fontsize',15,'linewidth',1.5,'box','on');
ylabel(ax1,'Latitude (\circ)','fontsize',15);
xlabel(ax1,'Longitude (\circ)','fontsize',15);
title(ax1,'Perturbation to ship location','FontWeight','bold','fontsize',18);
pos = get(gca,'Position');
cb = colorbar(ax1);
set(gca,'Position',pos);
ylabel(cb,'|\delta r_{ship}| (m)','fontsize',15);
colormap(ax1,parula)
caxis(ax1,[0 max(dr_ellipsoid)]);

%% Residual travel-time correction 
TWT_elli = 2 .* sqrt((x_elli-xOBS).^2+(y_elli-yOBS).^2+(z_elli-zOBS).^2) ./ vp + TAT;
TWT_sph = 2 .* sqrt((x_sph-xOBS).^2+(y_sph-yOBS).^2+(z_sph-zOBS).^2) ./ vp + TAT;
dTWT_ellipsoid = TWT_elli-TWT_sph;
azi(azi<1) = azi(azi<1)+360;

r_elli = sqrt((x_elli-xOBS).^2+(y_elli-yOBS).^2+(z_elli-zOBS).^2);
r_sph = sqrt((x_sph-xOBS).^2+(y_sph-yOBS).^2+(z_sph-zOBS).^2);
dr_ellipsoid2 = r_elli-r_sph;

subplot(2,1,2);
scatter(azi,dTWT_ellipsoid*1000,130,dr_ellipsoid2,'o','filled','MarkerEdgeColor',[0 0 0]);
set(gca,'fontsize',15,'linewidth',1.5,'box','on');
xlabel('Ship Azimuth (\circ)','fontsize',15);
ylabel('TWT_{ellip} - TWT_{sphere} (ms)','fontsize',15);
title('Perturbation to travel times','FontWeight','bold','fontsize',18);
ylim([-6 2]);
xlim([50 370]);
pos = get(gca,'Position');
cb2 = colorbar;
% set(gca,'Position',pos);
ylabel(cb2,'r_{ellip} - r_{sphere} (m)','fontsize',15);
caxis([-4 4]);
colormap(redbluecmap);

dTWT_p2p = (max(dTWT_ellipsoid)-min(dTWT_ellipsoid))*1000;

%%
disp(['Max/Min dr of ship due to ellipsoid correction: ',num2str(max(dr_ellipsoid)),'/',num2str(min(dr_ellipsoid)),' m'])
disp(['TWT residual peak-to-peak: ',num2str(dTWT_p2p),' ms'])

save2pdf('ellipsoid_corr.pdf',1,100);