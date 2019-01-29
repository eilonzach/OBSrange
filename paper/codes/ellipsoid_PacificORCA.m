% Plot example of residuals with and without ellipsoid correction for
% Pacific ORCA
clear

addpath('/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/OBSrange_v1_MATLAB/functions');

staname = 'CC06';

% Load corrected station
datacorr =   load(['/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_noTAT_ellipse/OUT_OBSrange/OUT_nocorr/mats/',staname,'_data.mat']);
datanocorr = load(['/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_noTAT_sphere/OUT_OBSrange/OUT_nocorr/mats/',staname,'_data.mat']);


% Calculate ship distance azimuth from OBS
dx_ship = datacorr.datamat.x_ship;
dy_ship = datacorr.datamat.y_ship;
dist = sqrt(dx_ship.^2 + dy_ship.^2);
azi_ship = -atan2d(dy_ship , dx_ship) + 90;
azi_ship(azi_ship<0) = azi_ship(azi_ship<0) + 360;

%% Plot survey patterns
figure(1); clf;
set(gcf,'Position',[237.0000   84.0000  683.0000  615.0000]);

cmap = flip(brewermap(100,'RdBu'));

% No correction
subplot(2,2,1);
lon = datanocorr.datamat.lons_ship;
lat = datanocorr.datamat.lats_ship;
dtwt = datanocorr.datamat.dtwt_bs;
scatter(lon,lat,80,dtwt*1000,'o','filled','MarkerEdgeColor',[0 0 0]); hold on;
axis equal;
ylim([min(lat)-(max(abs(lat))-min(abs(lat)))*0.05 , max(lat)+(max(abs(lat))-min(abs(lat)))*0.05]);
xlim([min(lon)-(max(abs(lon))-min(abs(lon)))*0.05 , max(lon)+(max(abs(lon))-min(abs(lon)))*0.05]);
set(gca,'fontsize',15,'linewidth',1.5,'box','on');
ylabel('Latitude (\circ)','fontsize',15);
xlabel('Longitude (\circ)','fontsize',15);
title('No correction','FontWeight','bold','fontsize',18);
caxis([-6 6]);

% Corrected
subplot(2,2,2);
lon = datacorr.datamat.lons_ship;
lat = datacorr.datamat.lats_ship;
dtwt = datacorr.datamat.dtwt_bs;
scatter(lon,lat,80,dtwt*1000,'o','filled','MarkerEdgeColor',[0 0 0]); hold on;
axis equal;
ylim([min(lat)-(max(abs(lat))-min(abs(lat)))*0.05 , max(lat)+(max(abs(lat))-min(abs(lat)))*0.05]);
xlim([min(lon)-(max(abs(lon))-min(abs(lon)))*0.05 , max(lon)+(max(abs(lon))-min(abs(lon)))*0.05]);
set(gca,'fontsize',15,'linewidth',1.5,'box','on','YTickLabel',[]);
% ylabel('Latitude (\circ)','fontsize',15);
xlabel('Longitude (\circ)','fontsize',15);
title('Corrected','FontWeight','bold','fontsize',18);
pos = get(gca,'Position');
cb = colorbar;
set(gca,'Position',[pos(1)-0.06 pos(2:4)]);
ylabel(cb,'\delta TWT (ms)','fontsize',15);
set(cb,'LineWidth',1.5);
caxis([-6 6]);
colormap(cmap);

%% Plot residuals as a function of ship location
subplot(2,2,[3 4]);
azi_ship(azi_ship>150) = azi_ship(azi_ship>150)-360;
plot([min(azi_ship) max(azi_ship)],[0 0],'--k','linewidth',1.5); hold on;
h(1) = plot(azi_ship,datanocorr.datamat.dtwt_bs*1000,'ok','markerfacecolor',[0 0 0],'markersize',11); hold on;
h(2) = plot(azi_ship,datacorr.datamat.dtwt_bs*1000,'ok','markerfacecolor',[0.7 0.7 0.7],'markersize',11); hold on;
set(gca,'fontsize',15,'linewidth',1.5,'box','on');
xlabel('Ship Azimuth (\circ)','fontsize',15);
ylabel('\delta TWT (ms)','fontsize',15);
title('Travel time residuals','FontWeight','bold','fontsize',18);
axis tight;
ylim([-7 7]);
% xlim([50 370]);
legend(h,{'no correction','corrected'},'location','southeast');

save2pdf(['ellipsoid_PacificORCA_',staname,'.pdf'],1,100);