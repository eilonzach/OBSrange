%% Function to produce Figure 01 
%  Figure 01 contains the results of a synthetic test for a single station,
%  comparing the true input model parameters to the inverted solution
function Figure01

ofile = '../Figure01';

systa = 'syn10'; % name of synthetic station

ifsave = 1; % Save pdf?
Nbins = 15; % Bins for histogram plots
col = colormap(parula); % colormap for f test plots
trucol = [0.2 0.9 0.6];
close



%% load 
addpath('~/Work/OBSrange/OBSrange_v1_MATLAB/functions/');
trudata = load(['../figdata/trudata_',systa,'.mat']);
load(['../figdata/',systa,'_data.mat']);

% % Data from Figure 5
% trudata = load(['../figdata/PacificORCA_synthtest4/trudata_syn12.mat']);
% load(['../figdata/PacificORCA_synthtest4/OUT_OBSrange/2_OUT_wcorr_xrec/mats/syn12_data.mat']);

%% ---------------------   PLOTTING   ---------------------   
%% set up figure
f901 = figure(901); clf;
set(f901,'position',[3 20 1208 797],'color','white');

% set up axes
dy = 0.05;
ax1 = axes(f901,'pos',[0.06 0.12-dy 0.245 0.36]); hold(ax1,'on');
ax2 = axes(f901,'pos',[0.37 0.12-dy 0.245 0.36]); hold(ax2,'on');
ax3 = axes(f901,'pos',[0.68 0.12-dy 0.305 0.36]); hold(ax3,'on');
ax4 = axes(f901,'pos',[0.07 0.57 0.25 0.33]); hold(ax4,'on');
ax5 = axes(f901,'pos',[0.38 0.57 0.25 0.33]);
ax6 = axes(f901,'pos',[0.68 0.57 0.25 0.33]);

%% HISTOGRAMS

% [ Ncount_TAT, cent_TAT ] = plot_hist(ax4,datamat.TAT_bs*1000,Nbins);
% title(ax4,'\textbf{TAT (ms)}','fontsize',18,'interpreter','latex');

[ Ncount_vw, cent_vw ] = plot_hist(ax5,datamat.V_w_bs,Nbins);
title(ax5,'\textbf{Water Velocity (m/s)}','fontsize',18,'interpreter','latex');

[ Ncount_drift, cent_drift ] = plot_hist(ax6,datamat.drift_bs,Nbins);
title(ax6,'\textbf{Drift (m)}','fontsize',18,'interpreter','latex');

% true values
% plot(ax4,1e3*trudata.tat*[1 1],axlim(ax4,[3,4]),'linewidth',3,'color',trucol)
plot(ax5,1e3*trudata.vp_actual*[1 1],axlim(ax5,[3,4]),'linewidth',3,'color',trucol)
plot(ax6,sqrt(2)*1e3*rms(trudata.obs_location_xyz(1:2))*[1 1],axlim(ax6,[3,4]),'linewidth',3,'color',trucol)

% drift directions
ax7 = polaraxes('pos',[0.625 0.34+0.45 0.16 0.16]); hold on
driftaz = atan2d(datamat.x_sta_bs,datamat.y_sta_bs);
plot_drift_polar( ax7,driftaz );
% draw drift azimuth
rlim = get(ax7,'Rlim'); 
truaz = atan2d(trudata.obs_location_xyz(1),trudata.obs_location_xyz(2));
polarplot(ax7,d2r(truaz+[0 0 7 0 -7]),rlim(2)*.97*[0 0.97 0.85 0.97 0.85],':','color',trucol,'linewidth',2)
polarplot(ax7,0,0,'^','linewidth',2,'markersize',12,'markerfacecolor',[ 0.1980 0.9500 0.7250],'markeredgecolor','k'); % put station back on

% remove ax6 tick underneath the polar plot
ytl=get(ax6,'YTicklabel');set(ax6,'YTickLabel',{ytl{1:end-1},''});




%% F-test plots
x_grid = datamat.Ftest_res.x_grid;
y_grid = datamat.Ftest_res.y_grid;
z_grid = datamat.Ftest_res.z_grid;
P = datamat.Ftest_res.Pstat;
[Xgrd,Ygrd,Zgrd] = meshgrid(x_grid,y_grid,z_grid);
[Pz_max, Iz_max] = max(max(max(P)));
[Py_max, Iy_max] = max(max(P(:,:,Iz_max)));
[Px_max, Ix_max] = max(P(:,Iy_max,Iz_max));


% bounds and ticks (zoom in within test space to fill plot)
x_bds = [min(x_grid),max(x_grid)] + [25 -25];
y_bds = [min(y_grid),max(y_grid)] + [25 -25];
z_bds = [min(z_grid),max(z_grid)];

% centre points:
x_mid = median(x_grid);
y_mid = median(y_grid);
z_mid = median(z_grid);

% bounds and ticks
xlim1 = [-20 20]; dx1 = 4;
ylim1 = [-20 20]; dy1 = 4;
xlim2 = [-40 40]; dx2 = 8;
ylim2 = [-40 40]; dy2 = 8;
zlim = [-40 40];  dz = 8;

% shade in background
fill(ax1,[-2000,2000,2000,-2000,-2000],[-2000,-2000,2000,2000,-2000],col(1,:))
fill(ax2,[-2000,2000,2000,-2000,-2000],[-8000,-8000,8000,8000,-8000],col(1,:))
fill(ax3,[-2000,2000,2000,-2000,-2000],[-8000,-8000,8000,8000,-8000],col(1,:))
ang = [0:0.1:2*pi];
%X-Y
contourf(ax1,Xgrd(:,:,Iz_max)-x_mid,Ygrd(:,:,Iz_max)-y_mid,P(:,:,Iz_max),'linestyle','none');
% shading(ax1,'flat');
contour(ax1,Xgrd(:,:,Iz_max)-x_mid,Ygrd(:,:,Iz_max)-y_mid,P(:,:,Iz_max),[[0.05 0.05],[0.32 0.32]],'-w','linewidth',2); hold on;
plot(ax1,Xgrd(Ix_max,Iy_max,Iz_max)-x_mid,Ygrd(Ix_max,Iy_max,Iz_max)-y_mid,'sk','markerfacecolor',[1 0 0],'markersize',12,'linewidth',1)
set(ax1,'fontsize',16,'linewidth',2,'box','on','layer','top',...
    'xlim',xlim1,'ylim',ylim1,'xtick',[xlim1(1):dx1:xlim1(2)],'ytick',[ylim1(1):dy1:ylim1(2)],'TickDir','out');
xlabel(ax1,'$\delta$X (m)','fontsize',18,'interpreter','latex');
ylabel(ax1,'$\delta$Y (m)','fontsize',18,'interpreter','latex');
title(ax1,'\textbf{X-Y}','fontsize',18,'interpreter','latex');
% average values
htx1 = text(ax1,xlim1(1)+0.04*diff(axlim(ax1,1:2)),ylim1(1)+0.11*diff(axlim(ax1,3:4)),...
    sprintf('$\\mathbf{\\bar{x}}$\\textbf{ = %.1f m}',x_mid),'color','white','interpreter','latex','fontsize',17);
hty1 = text(ax1,xlim1(1)+0.04*diff(axlim(ax1,1:2)),ylim1(1)+0.05*diff(axlim(ax1,3:4)),...
    sprintf('$\\mathbf{\\bar{y}}$\\textbf{ = %.1f m}',y_mid),'color','white','interpreter','latex','fontsize',17);
% axis equal;

%X-Z
contourf(ax2,squeeze(Xgrd(Iy_max,:,:))-x_mid,squeeze(Zgrd(Iy_max,:,:))-z_mid,squeeze(P(Iy_max,:,:)),'linestyle','none');
% shading(ax1,'flat');
contour(ax2,squeeze(Xgrd(Iy_max,:,:))-x_mid,squeeze(Zgrd(Iy_max,:,:))-z_mid,squeeze(P(Iy_max,:,:)),[[0.05 0.05],[0.32 0.32]],'-w','linewidth',2); hold on;
plot(ax2,Xgrd(Ix_max,Iy_max,Iz_max)-x_mid,Zgrd(Ix_max,Iy_max,Iz_max)-z_mid,'sk','markerfacecolor',[1 0 0],'markersize',12,'linewidth',1)
set(ax2,'fontsize',16,'linewidth',2,'box','on','layer','top',...
    'xlim',xlim2,'ylim',zlim,'xtick',[xlim2(1):dx2:xlim2(2)],'ytick',[zlim(1):dz:zlim(2)],'TickDir','out');
xlabel(ax2,'$\delta$X (m)','fontsize',18,'interpreter','latex');
ylabel(ax2,'$\delta$Z (m)','fontsize',18,'interpreter','latex');
title(ax2,'\textbf{X-Z}','fontsize',18,'interpreter','latex');
% average values
htx2 = text(ax2,xlim2(1)+0.04*diff(axlim(ax2,1:2)),zlim(1)+0.11*diff(axlim(ax2,3:4)),...
    sprintf('$\\mathbf{\\bar{x}}$\\textbf{ = %.1f m}',x_mid),'color','white','interpreter','latex','fontsize',17);
htz2 = text(ax2,xlim2(1)+0.04*diff(axlim(ax2,1:2)),zlim(1)+0.05*diff(axlim(ax2,3:4)),...
    sprintf('$\\mathbf{\\bar{z}}$\\textbf{ = %.1f m}',z_mid),'color','white','interpreter','latex','fontsize',17);
% axis equal;

%Y-Z
contourf(ax3,squeeze(Ygrd(:,Ix_max,:))-y_mid,squeeze(Zgrd(:,Ix_max,:))-z_mid,squeeze(P(:,Ix_max,:)),'linestyle','none');
% shading(ax1,'flat');
contour(ax3,squeeze(Ygrd(:,Ix_max,:))-y_mid,squeeze(Zgrd(:,Ix_max,:))-z_mid,squeeze(P(:,Ix_max,:)),[[0.05 0.05],[0.32 0.32]],'-w','linewidth',2); hold on;
plot(ax3,Ygrd(Ix_max,Iy_max,Iz_max)-y_mid,Zgrd(Ix_max,Iy_max,Iz_max)-z_mid,'sk','markerfacecolor',[1 0 0],'markersize',12,'linewidth',1)
set(ax3,'fontsize',16,'linewidth',2,'box','on','layer','top',...
    'xlim',ylim2,'ylim',zlim,'xtick',[ylim2(1):dy2:ylim2(2)],'ytick',[zlim(1):dz:zlim(2)],'TickDir','out');
xlabel(ax3,'$\delta$Y (m)','fontsize',18,'interpreter','latex');
ylabel(ax3,'$\delta$Z (m)','fontsize',18,'interpreter','latex');
title(ax3,'\textbf{Y-Z}','fontsize',18,'interpreter','latex');
% average values
hty3 = text(ax3,ylim2(1)+0.04*diff(axlim(ax3,1:2)),zlim(1)+0.11*diff(axlim(ax3,3:4)),...
    sprintf('$\\mathbf{\\bar{y}}$\\textbf{ = %.1f m}',y_mid),'color','white','interpreter','latex','fontsize',17);
htz3 = text(ax3,ylim2(1)+0.04*diff(axlim(ax3,1:2)),zlim(1)+0.05*diff(axlim(ax3,3:4)),...
    sprintf('$\\mathbf{\\bar{z}}$\\textbf{ = %.1f m}',z_mid),'color','white','interpreter','latex','fontsize',17);
% axis equal;
colorbar('peer',ax3)% plot the results of the F-tests for error in the final inverted location.

% Add on the true location for comparison
% plot the true location on the F-test plot in each dimension
plot(ax1,trudata.obs_location_xyz(1)*1e3 - x_mid,trudata.obs_location_xyz(2)*1e3 - y_mid,'p',...
    'markerfacecolor',trucol,'markeredgecolor','k','linewidth',1.5,'markersize',13)
plot(ax2,trudata.obs_location_xyz(1)*1e3 - x_mid,trudata.obs_location_xyz(3)*-1e3 - z_mid,'p',...
    'markerfacecolor',trucol,'markeredgecolor','k','linewidth',1.5,'markersize',13)
plot(ax3,trudata.obs_location_xyz(2)*1e3 - y_mid,trudata.obs_location_xyz(3)*-1e3 - z_mid,'p',...
    'markerfacecolor',trucol,'markeredgecolor','k','linewidth',1.5,'markersize',13)

%% Survey pattern & residuals
latlim = [min(datamat.y_ship) max(datamat.y_ship)]+[-1 1]*200; % 200m outside ship circle
lonlim = [min(datamat.x_ship) max(datamat.x_ship)]+[-1 1]*200; % 200m outside ship circle
resid_max = prctile(abs(datamat.dtwt_bs*1000),100)*1;

trudata.obs_location_xyz

scatter(ax4,datamat.x_ship,datamat.y_ship,100,datamat.dtwt_bs*1000,'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on;
plot(ax4,0,0,'sk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1); hold on;
plot(ax4,datamat.loc_xyz(1),datamat.loc_xyz(2),'sk','markerfacecolor',[1 0 0],'markersize',18,'linewidth',1);
plot(ax4,trudata.obs_location_xyz(1)*1000,trudata.obs_location_xyz(2)*1000,'p','markerfacecolor',trucol,'markeredgecolor','k','linewidth',1.5,'markersize',18);

% title(ax2,'\textbf{Travel time residuals (ms)}','fontsize',18,'interpreter','latex');
xlabel(ax4,'X (m)','fontsize',18,'interpreter','latex');
ylabel(ax4,'Y (m)','fontsize',18,'interpreter','latex');
title(ax4,'\textbf{Survey pattern}','fontsize',18,'interpreter','latex');

set(ax4,'fontsize',16,'linewidth',2,'box','on',...
        'xlim',lonlim,'ylim',latlim); grid(ax4,'on');
cmap = flip(brewermap(100,'RdBu'));
colormap(ax4,cmap);
hcb2 = colorbar('peer',ax4); set(hcb2,'linewidth',2)
ylabel(hcb2,'\textbf{Travel time residuals (ms)}','fontsize',16,'interpreter','latex','rotation',270,'position',[3.1190 mean(get(hcb2,'limits')) 0]);
caxis(ax4,[-resid_max resid_max])
% cb3.Ticks = linspace(-max(cb3.Ticks),max(cb3.Ticks),max(cb3.Ticks)+1);
axis(ax4,'equal');
xlim(ax4,lonlim);
ylim(ax4,latlim);
% set axis limits
% set(ax4,'xlim',mean(cent_TAT) + 1.3*([min(cent_TAT),max(cent_TAT)]-mean(cent_TAT)));

%%
% delete(ax4);
dx_shift = 0;%-0.16;
dx = 1;
dy_shift = 0;%-0.04;
ax5.Position = [ax5.Position(1)+dx_shift, ax5.Position(2)+dy_shift, ax5.Position(3)*dx, ax5.Position(4)];
ax6.Position = [ax6.Position(1)+dx_shift, ax6.Position(2)+dy_shift, ax6.Position(3)*dx, ax6.Position(4)];
ax7.Position = [ax7.Position(1)+dx_shift, ax7.Position(2)+dy_shift, ax7.Position(3)*dx, ax7.Position(4)];

%% Figure numbers
x = -0.1;
y = 1.1;
figN_add('d)',ax1,x,y,25);
figN_add('e)',ax2,x,y,25);
figN_add('f)',ax3,x,y,25);
figN_add('a)',ax4,x,y+0.06,25);
figN_add('b)',ax5,x,y,25);
figN_add('c)',ax6,x,y,25);

%% SAVE
if ifsave
save2pdf(ofile,f901)
end

end