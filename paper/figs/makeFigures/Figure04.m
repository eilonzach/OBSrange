%% Function to produce Figure 04 
%  Figure 04 shows the Ftest plots
function Figure04

ofile = '../Figure04';
ifsave = 1;

sta = 'EC03';

col = colormap(parula); % colormap for f test plots

%% load 
load(['../figdata/',sta,'_data.mat']);

%% ---------------------   PLOTTING   ---------------------   

%% F-plot plots
f904 = figure(904); clf;
set(f904,'position',[3 602 1208 396]);
ax1 = axes(f904,'pos',[0.06 0.15 0.25 0.76]); hold(ax1,'on');
ax2 = axes(f904,'pos',[0.37 0.15 0.25 0.76]); hold(ax2,'on');
ax3 = axes(f904,'pos',[0.68 0.15 0.30 0.76]); hold(ax3,'on');

x_grid = datamat.Ftest_res.x_grid;
y_grid = datamat.Ftest_res.y_grid;
z_grid = datamat.Ftest_res.z_grid;
P = datamat.Ftest_res.Pstat;
[Xgrd,Ygrd,Zgrd] = meshgrid(x_grid,y_grid,z_grid);
[Pz_max, Iz_max] = max(max(max(P)));
[Py_max, Iy_max] = max(max(P(:,:,Iz_max)));
[Px_max, Ix_max] = max(P(:,Iy_max,Iz_max));


% bounds and ticks
% centre points:
x_mid = median(x_grid);
y_mid = median(y_grid);
z_mid = median(z_grid);

% bounds and ticks
xlim1 = [-8 8]; dx1 = 2;
ylim1 = [-8 8]; dy1 = 2;
xlim2 = [-15 15]; dx2 = 5;
ylim2 = [-15 15]; dy2 = 5;
zlim = [-15 15];  dz = 5;


% shade in background
fill(ax1,[-2000,2000,2000,-2000,-2000],[-2000,-2000,2000,2000,-2000],col(1,:))
fill(ax2,[-2000,2000,2000,-2000,-2000],[-8000,-8000,8000,8000,-8000],col(1,:))
fill(ax3,[-2000,2000,2000,-2000,-2000],[-8000,-8000,8000,8000,-8000],col(1,:))

%% X-Y
contourf(ax1,Xgrd(:,:,Iz_max)-x_mid,Ygrd(:,:,Iz_max)-y_mid,P(:,:,Iz_max),'linestyle','none');
plot(ax1,datamat.x_sta_bs-x_mid,datamat.y_sta_bs-y_mid,'.k','MarkerSize',5)
% shading(ax1,'flat');
contour(ax1,Xgrd(:,:,Iz_max)-x_mid,Ygrd(:,:,Iz_max)-y_mid,P(:,:,Iz_max),[[0.05 0.05],[0.32 0.32]],'-w','linewidth',2); hold on;
plot(ax1,Xgrd(Ix_max,Iy_max,Iz_max)-x_mid,Ygrd(Ix_max,Iy_max,Iz_max)-y_mid,'pk','markerfacecolor',[1 0 0],'markersize',20,'linewidth',1)
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


%% X-Z
contourf(ax2,squeeze(Xgrd(Iy_max,:,:))-x_mid,squeeze(Zgrd(Iy_max,:,:))-z_mid,squeeze(P(Iy_max,:,:)),'linestyle','none');
plot(ax2,datamat.x_sta_bs-x_mid,datamat.z_sta_bs-z_mid,'.k','MarkerSize',5)
% shading(ax1,'flat');
contour(ax2,squeeze(Xgrd(Iy_max,:,:))-x_mid,squeeze(Zgrd(Iy_max,:,:))-z_mid,squeeze(P(Iy_max,:,:)),[[0.05 0.05],[0.32 0.32]],'-w','linewidth',2); hold on;
plot(ax2,Xgrd(Ix_max,Iy_max,Iz_max)-x_mid,Zgrd(Ix_max,Iy_max,Iz_max)-z_mid,'pk','markerfacecolor',[1 0 0],'markersize',20,'linewidth',1)
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

%% Y-Z
contourf(ax3,squeeze(Ygrd(:,Ix_max,:))-y_mid,squeeze(Zgrd(:,Ix_max,:))-z_mid,squeeze(P(:,Ix_max,:)),'linestyle','none');
plot(ax3,datamat.y_sta_bs-y_mid,datamat.z_sta_bs-z_mid,'.k','MarkerSize',5)
% shading(ax1,'flat');
contour(ax3,squeeze(Ygrd(:,Ix_max,:))-y_mid,squeeze(Zgrd(:,Ix_max,:))-z_mid,squeeze(P(:,Ix_max,:)),[[0.05 0.05],[0.32 0.32]],'-w','linewidth',2); hold on;
plot(ax3,Ygrd(Ix_max,Iy_max,Iz_max)-y_mid,Zgrd(Ix_max,Iy_max,Iz_max)-z_mid,'pk','markerfacecolor',[1 0 0],'markersize',20,'linewidth',1)
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

%% colourbar
hcb = colorbar('peer',ax3); 
set(hcb,'linewidth',2)



%% SAVE
if ifsave
    save2pdf(ofile,f904)
end

end
