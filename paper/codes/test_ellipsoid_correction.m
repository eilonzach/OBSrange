% Test how important the ellipsoid correction is...
clear

addpath('/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/OBSrange_v1_MATLAB/functions');

% Load data
trudata = load('/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/paper/figs/figdata/PacificORCA_synthtest4/trudata_syn12.mat');

% dlat = [7:95];
dlat = 0; % add to latitude to see how it changes things...
olon = trudata.drop_location(2); %-133.62; %0:360;
olat = trudata.drop_location(1); %-7.54;
lon = trudata.survlon;
lat = trudata.survlat;
TAT = trudata.tat;
xOBS = trudata.obs_location_xyz(1);
yOBS = trudata.obs_location_xyz(2);
zOBS = -trudata.obs_location_xyz(3)*1000;
vp = trudata.vp_actual*1000;
azi = trudata.survaz;

[ x_elli, y_elli , z_elli] = lonlat2xy_nomap( olon, olat+dlat, lon, lat+dlat );
[ x_sph, y_sph , z_sph] = lonlat2xy_nomap_sphere( olon, olat+dlat, lon, lat+dlat );

%% Distance correction
dr_ellipsoid = sqrt((y_sph-y_elli).^2+(x_sph-x_elli).^2);

figure(1); clf;
subplot(2,1,1);
scatter(lon,lat,80,dr_ellipsoid,'o','filled','MarkerEdgeColor',[0 0 0]); hold on;
axis equal;
ylim([min(lat)-(max(abs(lat))-min(abs(lat)))*0.05 , max(lat)+(max(abs(lat))-min(abs(lat)))*0.05]);
xlim([min(lon)-(max(abs(lon))-min(abs(lon)))*0.05 , max(lon)+(max(abs(lon))-min(abs(lon)))*0.05]);
set(gca,'fontsize',16,'linewidth',1.5,'box','on');
cb = colorbar;
ylabel(cb,'Ellipsoid \delta r (m)','FontWeight','bold','fontsize',18);
caxis([0 max(dr_ellipsoid)]);

%% Residual travel-time correction 
TWT_elli = 2 .* sqrt((x_elli-xOBS).^2+(y_elli-yOBS).^2+(zOBS).^2) ./ vp + TAT;
TWT_sph = 2 .* sqrt((x_sph-xOBS).^2+(y_sph-yOBS).^2+(zOBS).^2) ./ vp + TAT;
dTWT_ellipsoid = TWT_elli-TWT_sph;

subplot(2,1,2);
plot(azi,dTWT_ellipsoid*1000,'ok','MarkerFaceColor',[0.5 0.5 0.5],'markersize',15);
set(gca,'fontsize',16,'linewidth',1.5,'box','on');
xlabel('Azimuth');
ylabel('Ellipsoid \delta TWT (ms)','FontWeight','bold','fontsize',18);

dTWT_p2p = (max(dTWT_ellipsoid)-min(dTWT_ellipsoid))*1000;

%%
disp(['Max/Min dr of ship due to ellipsoid correction: ',num2str(max(dr_ellipsoid)),'/',num2str(min(dr_ellipsoid)),' m'])
disp(['TWT residual peak-to-peak: ',num2str(dTWT_p2p),' ms'])
