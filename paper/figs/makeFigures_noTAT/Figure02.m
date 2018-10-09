%% Function to produce Figure 02 
%  Figure 02 shows the survey pattern and residuals for an example station
function Figure02

ofile = '../Figure02';
ifsave = 1;

sta = 'EC03';

%% load 
load(['../figdata/',sta,'_data.mat']);


%% ---------------------   PLOTTING   ---------------------   
%% set up figure
f902 = figure(902); clf;
set(f902,'position',[903 20 1108 797],'color','white');

% axes
ax1 = axes(f902,'pos',[0.08 0.54 0.40 0.40]);  hold(ax1,'on');
ax2 = axes(f902,'pos',[0.08 0.08 0.40 0.40]);  hold(ax2,'on');
ax3 = axes(f902,'pos',[0.58 0.55 0.37 0.38]); hold(ax3,'on');
ax4 = axes(f902,'pos',[0.58 0.08 0.37 0.36]);  hold(ax4,'on');


% limits
latlim = [min(datamat.lats_ship) max(datamat.lats_ship)]+[-1 1]*.2/111; % 200m outside ship circle
lonlim = [min(datamat.lons_ship) max(datamat.lons_ship)]+[-1 1]*.2/111; % 200m outside ship circle
londist = vdist(mean(datamat.lats_ship),lonlim(1),mean(datamat.lats_ship),lonlim(2));
resid_max = prctile(abs(datamat.dtwt_bs*1000),97);

% Divergent colormap
% cmap = cmocean('balance');
cmap = flip(brewermap(100,'RdBu'));

%% Plot Velocities
plot(ax1,datamat.databad.lons,datamat.databad.lats,'xk','markersize',12,'linewidth',2.2);
plot(ax1,datamat.drop_lonlatz(1),datamat.drop_lonlatz(2),'sk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1);
plot(ax1,datamat.loc_lolaz(1),datamat.loc_lolaz(2),'pk','markerfacecolor',[1 0 0],'markersize',21,'linewidth',1);
scatter(ax1,datamat.lons_ship,datamat.lats_ship,120,sqrt(sum(datamat.v_ship.^2,1)),'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on;

colormap(ax1,summer)
title(ax1,'\textbf{Survey pattern}','fontsize',18,'interpreter','latex');
% xlabel(ax1,'Longitude','fontsize',18,'interpreter','latex');
ylabel(ax1,'Latitude','fontsize',18,'interpreter','latex');

set(ax1,'fontsize',16,'linewidth',2,'box','on',...
        'xlim',lonlim,'ylim',latlim); grid(ax1,'on');
hcb1 = colorbar('peer',ax1); set(hcb1,'linewidth',2)
ylabel(hcb1,'\textbf{Ship Velocity (m/s)}','fontsize',16,'interpreter','latex','rotation',270,'position',[3.5190 mean(get(hcb1,'limits')) 0]);
axis(ax1,'equal');
xlim(ax1,lonlim);
ylim(ax1,latlim);


%% Plot residuals at each site
plot(ax2,datamat.drop_lonlatz(1),datamat.drop_lonlatz(2),'sk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1);
plot(ax2,datamat.loc_lolaz(1),datamat.loc_lolaz(2),'pk','markerfacecolor',[1 0 0],'markersize',21,'linewidth',1);
scatter(ax2,datamat.lons_ship,datamat.lats_ship,100,datamat.dtwt_bs*1000,'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on;

% title(ax2,'\textbf{Travel time residuals (ms)}','fontsize',18,'interpreter','latex');
xlabel(ax2,'Longitude','fontsize',18,'interpreter','latex');
ylabel(ax2,'Latitude','fontsize',18,'interpreter','latex');

set(ax2,'fontsize',16,'linewidth',2,'box','on',...
        'xlim',lonlim,'ylim',latlim); grid(ax2,'on');
colormap(ax2,cmap);
hcb2 = colorbar('peer',ax2); set(hcb2,'linewidth',2)
ylabel(hcb2,'\textbf{Travel time residuals (ms)}','fontsize',16,'interpreter','latex','rotation',270,'position',[3.1190 mean(get(hcb2,'limits')) 0]);
caxis(ax2,[-resid_max resid_max])
% cb3.Ticks = linspace(-max(cb3.Ticks),max(cb3.Ticks),max(cb3.Ticks)+1);
axis(ax2,'equal');
xlim(ax2,lonlim);
ylim(ax2,latlim);


%% Plot residuals as a function of azimuth + Doppler correction
azi_locss = -atan2d(datamat.y_ship-datamat.loc_xyz(2), datamat.x_ship-datamat.loc_xyz(1)) + 90;
azi_locss(azi_locss<0) = azi_locss(azi_locss<0) + 360;
obs_num = 1:length(datamat.dtwt_bs);
scatter(ax3,azi_locss,datamat.dtwt_bs*1000,100,datamat.dtwtcorr_bs*1000,'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on

title(ax3,'\textbf{Travel time Residuals (ms)}','fontsize',18,'interpreter','latex');
xlabel(ax3,'Ship Azimuth (deg)','fontsize',18,'interpreter','latex');
ylabel(ax3,'Travel time residuals (ms)','fontsize',18,'interpreter','latex');
set(ax3,'fontsize',16,'linewidth',2,'box','on',...
    'xlim',[0 360],'xtick',[0:45:360],'ylim',4.5*[-1 1]); 
grid(ax3,'on');
colormap(ax3,cmap)
hcb3 = colorbar('peer',ax3); set(hcb3,'linewidth',2)
caxis(ax3,[-max(abs(datamat.dtwtcorr_bs)*1000) max(abs(datamat.dtwtcorr_bs)*1000)])
hcb3.Ticks = linspace(-max(hcb3.Ticks),max(hcb3.Ticks),11);
ylabel(hcb3,'\textbf{Doppler corrections (ms)}','interpreter','latex','fontsize',16,'interpreter','latex','rotation',270,'position',[3.1 mean(get(hcb3,'limits')) 0]);

%% Misfit histogram
Nbins = 15;
plot_hist(ax4,datamat.E_rms*1000,Nbins);
plot(ax4,rms(datamat.dtwt_bs)*1000,axlim(ax4,4)*.03,'pk','markerfacecolor',[1 0 0],'markersize',21,'linewidth',1);
set(ax4,'fontsize',16,'linewidth',2); box on;
title(ax4,'\textbf{Bootstrap misfits}','fontsize',18,'interpreter','latex');
xlabel(ax4,'RMS (ms)','fontsize',18,'interpreter','latex');

%% Figure numbers
figN_add('a)',ax1,-0.25,0.95,25);
figN_add('b)',ax2,-0.25,0.95,25);
figN_add('c)',ax3,-0.15,0.97,25);
figN_add('d)',ax4,-0.123,1.057,25);


%% SAVE
if ifsave
    save2pdf(ofile,f902)
end

end