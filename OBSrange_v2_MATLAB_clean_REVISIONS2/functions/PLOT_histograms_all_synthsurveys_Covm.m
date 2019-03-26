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
dx = 0.75;
dy = 0.9;

covm_95 = sqrt(diag(Cm_mat(:,:,1)))*2;

% make the histograms
[ Ncount_lat, cent_lat ] = plot_hist(ax1,data(1).misfit_xsta,Nbins);
Ncount_norm = Ncount_lat./sum(Ncount_lat);
plot(ax1,median(data(1).misfit_xsta)+[covm_95(1) covm_95(1)],[0 1.2*max(Ncount_norm)],':','color',[160, 57, 239]/255,'linewidth',3);
plot(ax1,median(data(1).misfit_xsta)-[covm_95(1) covm_95(1)],[0 1.2*max(Ncount_norm)],':','color',[160, 57, 239]/255,'linewidth',3);
title(ax1,['$\mathbf{\delta X \, (m)}$'],'fontsize',18,'interpreter','latex');
text(ax1,diff(ax1.XLim)*dx+ax1.XLim(1),diff(ax1.YLim)*dy+ax1.YLim(1),...
    sprintf('\\textbf{ %.1f m}',rms(data(1).misfit_xsta)),'color',[0 0 0],'interpreter','latex','fontsize',14);

[ Ncount_lon, cent_lon ] = plot_hist(ax2,data(1).misfit_ysta,Nbins);
Ncount_norm = Ncount_lon./sum(Ncount_lon);
plot(ax2,median(data(1).misfit_ysta)+[covm_95(2) covm_95(2)],[0 1.2*max(Ncount_norm)],':','color',[160, 57, 239]/255,'linewidth',3);
plot(ax2,median(data(1).misfit_ysta)-[covm_95(2) covm_95(2)],[0 1.2*max(Ncount_norm)],':','color',[160, 57, 239]/255,'linewidth',3);
title(ax2,['$\mathbf{\delta Y \, (m)}$ '],'fontsize',18,'interpreter','latex');
text(ax2,diff(ax2.XLim)*dx+ax2.XLim(1),diff(ax2.YLim)*dy+ax2.YLim(1),...
    sprintf('\\textbf{ %.1f m}',rms(data(1).misfit_ysta)),'color',[0 0 0],'interpreter','latex','fontsize',14);

[ Ncount_z, cent_z ] = plot_hist(ax3,data(1).misfit_zsta,Nbins);
Ncount_norm = Ncount_z./sum(Ncount_z);
plot(ax3,median(data(1).misfit_zsta)+[covm_95(3) covm_95(3)],[0 1.2*max(Ncount_norm)],':','color',[160, 57, 239]/255,'linewidth',3);
plot(ax3,median(data(1).misfit_zsta)-[covm_95(3) covm_95(3)],[0 1.2*max(Ncount_norm)],':','color',[160, 57, 239]/255,'linewidth',3);
title(ax3,['$\mathbf{\delta Z \, (m)}$ '],'fontsize',18,'interpreter','latex');
text(ax3,diff(ax3.XLim)*dx+ax3.XLim(1),diff(ax3.YLim)*dy+ax3.YLim(1),...
    sprintf('\\textbf{ %.1f m}',rms(data(1).misfit_zsta)),'color',[0 0 0],'interpreter','latex','fontsize',14);

[ Ncount_TAT, cent_TAT ] = plot_hist(ax4,data(1).misfit_TAT*1000,Nbins);
title(ax4,['$\mathbf{\delta TAT \, (ms) }$ '],'fontsize',18,'interpreter','latex');
text(ax4,diff(ax4.XLim)*dx+ax4.XLim(1),diff(ax4.YLim)*dy+ax4.YLim(1),...
    sprintf('\\textbf{ %.1f ms}',rms(data(1).misfit_TAT*1000)),'color',[0 0 0],'interpreter','latex','fontsize',14);

[ Ncount_vw, cent_vw ] = plot_hist(ax5,data(1).misfit_Vw,Nbins);
Ncount_norm = Ncount_vw./sum(Ncount_vw);
plot(ax5,median(data(1).misfit_Vw)+[covm_95(4) covm_95(4)],[0 1.2*max(Ncount_norm)],':','color',[160, 57, 239]/255,'linewidth',3);
plot(ax5,median(data(1).misfit_Vw)-[covm_95(4) covm_95(4)],[0 1.2*max(Ncount_norm)],':','color',[160, 57, 239]/255,'linewidth',3);
title(ax5,['$\mathbf{\delta V_w \, (m/s)}$ '],'fontsize',18,'interpreter','latex');
text(ax5,diff(ax5.XLim)*dx+ax5.XLim(1),diff(ax5.YLim)*dy+ax5.YLim(1),...
    sprintf('\\textbf{ %.1f m/s}',rms(data(1).misfit_Vw)),'color',[0 0 0],'interpreter','latex','fontsize',14);

[ Ncount_drift, cent_drift ] = plot_hist(ax6,data(1).misfit_r_xyz,Nbins);
title(ax6,['$\mathbf{\delta r_{xyz} \, (m)}$ '],'fontsize',18,'interpreter','latex');
text(ax6,diff(ax6.XLim)*dx+ax6.XLim(1),diff(ax6.YLim)*dy+ax6.YLim(1),...
    sprintf('\\textbf{ %.1f m}',rms(data(1).misfit_r_xyz)),'color',[0 0 0],'interpreter','latex','fontsize',14);

%% ticks for the histograms - in m about central value
% [cent_x,cent_y] = lonlat2xy_nomap(mean(cent_lon),mean(cent_lat),cent_lon,cent_lat);
% % ticks every two m
% tick_x = [2*floor(min(cent_x)/2):2:2*ceil(max(cent_x)/2)];
% tick_y = [2*floor(min(cent_y)/2):2:2*ceil(max(cent_y)/2)];
% [~,tick_lat] = xy2lonlat_nomap(mean(cent_lon),mean(cent_lat),0,tick_y);
% [tick_lon,~] = xy2lonlat_nomap(mean(cent_lon),mean(cent_lat),tick_x,0);
% [~,latlim] = xy2lonlat_nomap(mean(cent_lon),mean(cent_lat),0,1.2*[min(cent_y),max(cent_y)]);
% [lonlim,~] = xy2lonlat_nomap(mean(cent_lon),mean(cent_lat),1.2*[min(cent_x),max(cent_x)],0);
% set(ax1,'xtick',tick_lat,'xticklabel',tick_y,'xlim',latlim)
% set(ax2,'xtick',tick_lon,'xticklabel',tick_x,'xlim',lonlim)
% xlabel(ax1,sprintf('m from %.5f$^{\\circ}$',mean(cent_lat)),'interpreter','latex');
% xlabel(ax2,sprintf('m from %.5f$^{\\circ}$',mean(cent_lon)),'interpreter','latex');
% set(ax3,'xlim',mean(cent_z) + 1.3*([min(cent_z),max(cent_z)]-mean(cent_z)));
% set(ax4,'xlim',mean(cent_TAT) + 1.3*([min(cent_TAT),max(cent_TAT)]-mean(cent_TAT)));


% %% drift directions
% ax7 = polaraxes('pos',[0.621 0.35 0.16 0.16]); hold on
% driftaz = atan2d(x_sta,y_sta);
% plot_drift_polar( ax7,driftaz );

% %% prettify
% % remove ax6 tick underneath the polar plot
% ytl=get(ax6,'YTicklabel');set(ax6,'YTickLabel',{ytl{1:end-1},''});

% set up figure
f101 = figure(101); clf;
set(gcf,'position',[3 1 1208 697],'color','white');

% set up axes
ax1 = axes('pos',[0.08 0.58 0.25 0.35]);
ax2 = axes('pos',[0.38 0.58 0.25 0.35]);
ax3 = axes('pos',[0.68 0.58 0.25 0.35]);
ax4 = axes('pos',[0.08 0.12 0.25 0.35]);
ax5 = axes('pos',[0.38 0.12 0.25 0.35]);
% ax6 = axes('pos',[0.68 0.12 0.25 0.35]);

% make the histograms
[ Ncount_lat, cent_lat ] = plot_hist(ax1,data(1).E_rms*1000,Nbins);
title(ax1,['$\mathbf{RMS \, (ms)}$'],'fontsize',18,'interpreter','latex');
text(ax1,diff(ax1.XLim)*dx+ax1.XLim(1),diff(ax1.YLim)*dy+ax1.YLim(1),...
    sprintf('\\textbf{ %.1f ms}',mean(data(1).E_rms*1000)),'color',[0 0 0],'interpreter','latex','fontsize',14);

[ Ncount_lon, cent_lon ] = plot_hist(ax2,data(1).misfit_v_ship_all(1,:),Nbins);
title(ax2,['$\mathbf{\delta V_x\, (m/s) }$'],'fontsize',18,'interpreter','latex');
text(ax2,diff(ax2.XLim)*dx+ax2.XLim(1),diff(ax2.YLim)*dy+ax2.YLim(1),...
    sprintf('\\textbf{ %.1f m/s}',rms(data(1).misfit_v_ship_all(1,:))),'color',[0 0 0],'interpreter','latex','fontsize',14);

[ Ncount_z, cent_z ] = plot_hist(ax3,data(1).misfit_v_ship_all(2,:),Nbins);
title(ax3,['$\mathbf{\delta V_y\, (m/s) }$'],'fontsize',18,'interpreter','latex');
text(ax3,diff(ax3.XLim)*dx+ax3.XLim(1),diff(ax3.YLim)*dy+ax3.YLim(1),...
    sprintf('\\textbf{ %.1f m/s}',rms(data(1).misfit_v_ship_all(2,:))),'color',[0 0 0],'interpreter','latex','fontsize',14);

[ Ncount_TAT, cent_TAT ] = plot_hist(ax4,data(1).dtwt_all*1000,Nbins);
title(ax4,['$\mathbf{\delta TWT\, (ms) }$ '],'fontsize',18,'interpreter','latex');
text(ax4,diff(ax4.XLim)*dx+ax4.XLim(1),diff(ax4.YLim)*dy+ax4.YLim(1),...
    sprintf('\\textbf{ %.1f ms}',rms(data(1).dtwt_all*1000)),'color',[0 0 0],'interpreter','latex','fontsize',14);

[ Ncount_vw, cent_vw ] = plot_hist(ax5,data(1).misfit_dtwtcorr_all*1000,Nbins);
title(ax5,['$\mathbf{\delta TWT_{corr} \,(ms)} $'],'fontsize',18,'interpreter','latex');
text(ax5,diff(ax5.XLim)*dx+ax5.XLim(1),diff(ax5.YLim)*dy+ax5.YLim(1),...
    sprintf('\\textbf{ %.1f ms}',rms(data(1).misfit_dtwtcorr_all*1000)),'color',[0 0 0],'interpreter','latex','fontsize',14);

% [ Ncount_drift, cent_drift ] = plot_hist(ax6,data(1).misfit_r_xyz,Nbins);
% title(ax6,['\mathbf{\delta r_{xyz} (m)} ',num2str(rms(data(1).misfit_r_xyz))],'fontsize',18,'interpreter','latex');




