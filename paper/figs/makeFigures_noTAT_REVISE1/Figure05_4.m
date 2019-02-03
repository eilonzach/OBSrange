%% Function to produce Figure 05 
%  Figure 05 compares our algorithm to SIO
function Figure05

ofile = '../Figure05';
ifsave = 1;

%% load 
iSIO = [7 8]+1;
data_dirs = {
    '2_OUT_wcorr_xrec';
    '1_OUT_nocorr';
    '7_OUT_nocorr_noellipsoid';
%     '3_OUT_nocorr_TAT';
    '4_OUT_nocorr_Vp';
    '5_OUT_nocorr_Z';
    '6_OUT_nocorr_TAT_Vp_Z';
    '8_SIO_compare_nobads';
    '9_SIO_compare_wbads';
    };

synth_dirs = {
    '2_OUT_wcorr_xrec';
    '1_OUT_nocorr';
    '7_OUT_wcorr_xrec_noellipsoid';
    '10_OUT_wcorr_xrec_noGPScorr';
%     '3_OUT_wcorr_xrec_TAT';
    '4_OUT_wcorr_xrec_Vp';
    '5_OUT_wcorr_xrec_Z';
    '6_OUT_wcorr_xrec_TAT_Vp_Z';
    '8_SIO_compare_nobads';
    '9_SIO_compare_wbads';
    };

% xlabels = {
%     'OBSrange';
%     'No Doppler';
%     'No Ellipsoid';
%     'XYZ$V_p$';
%     'XYZ$\tau$';
%     'XY$\tau V_p$';
%     'XY';
%     'SIOgs';
%     'SIOgs no QC';
%     };

xlabels = {
    'OBSrange (1)';
    'No Doppler (2)';
    'No Ellipsoid (3)';
    'No GPS (4)';
%     'Fix-$\tau$';
    'Fix-$V_p$ (5)';
    'Fix-Z (6)';
    'XY-only (7)';
    'SIOgs (8)';
    'SIOgs no QC (9)';
    };

symbols = {
    'pk';
    'ok';
    'ok';
    'ok';
%     'ok';
    'ok';
    'ok';
    'ok';
    'pk';
    'ok';
    };

sizes = [
    20;
    14;
    14;
    14;
%     14;
    14;
    14;
    14;
    20;
    14 ];

%           X Y Z TAT Vp
issolve = [ 1 1 1  1  1;
            1 1 1  1  1;
            1 1 1  1  1;
            1 1 1  1  1;
%             1 1 1  0  1;
            1 1 1  1  0;
            1 1 0  1  1;
            1 1 0  0  0;
            1 1 0  0  0;
            1 1 0  0  0 ];

data_path = '../figdata/PacificORCA_EC03/OUT_OBSrange';

% synth_path = '../figdata/PacificORCA_synthtest4_noTAT/OUT_OBSrange';
% % Load synthetic
% trudata = load('../figdata/PacificORCA_synthtest4_noTAT/trudata_syn12.mat');

% REVISION1
synth_path = '../figdata/PacificORCA_synthtest4_REVISION1_GPScorr/OUT_OBSrange';
trudata = load('../figdata/PacificORCA_synthtest4_REVISION1_GPScorr/trudata_syn12_z5000m_fr10.mat');

% Load data
for ifil = 1:length(data_dirs)
    data_mat = dir(fullfile(data_path,data_dirs{ifil},'mats/*.mat'));
    data = load(fullfile(data_path,data_dirs{ifil},'mats',data_mat.name));
    if ~isempty(data.datamat.E_rms)
        RMS_data(ifil) = median(data.datamat.E_rms);
        RMS_95(:,ifil) = abs(median(data.datamat.E_rms) - prctile(data.datamat.E_rms,[2.5 97.5]));
    else
        RMS_data(ifil) = rms(data.datamat.dtwt_bs);
        RMS_95(:,ifil) = [nan nan];
    end
end

for ifil = 1:length(synth_dirs)
    synth_mat = dir(fullfile(synth_path,synth_dirs{ifil},'mats/*.mat'));
    synth = load(fullfile(synth_path,synth_dirs{ifil},'mats',synth_mat.name));
    misfit_xsta(ifil) = median(synth.datamat.x_sta_bs) - trudata.obs_location_xyz(1)*1000;
    misfit_ysta(ifil) = median(synth.datamat.y_sta_bs) - trudata.obs_location_xyz(2)*1000;
    misfit_zsta(ifil) = median(synth.datamat.z_sta_bs) - (-trudata.obs_location_xyz(3)*1000);
%     if ~isempty(data.datamat.E_rms)
%         misfit_r_xy(ifil) = median(synth.datamat.E_rms);
%     else
%         misfit_r_xy(ifil) = rms(synth.datamat.dtwt_bs);
%     end

    misfit_xsta_bs(:,ifil) = synth.datamat.x_sta_bs - trudata.obs_location_xyz(1)*1000;
    misfit_ysta_bs(:,ifil) = synth.datamat.y_sta_bs - trudata.obs_location_xyz(2)*1000;
    misfit_zsta_bs(:,ifil) = synth.datamat.z_sta_bs - (-trudata.obs_location_xyz(3)*1000);
    misfit_zsta(:,ifil) = rms(misfit_zsta_bs(:,ifil));
    RMS_data(ifil) = mean(synth.datamat.E_rms);
    dtwt(:,ifil) = synth.datamat.dtwt_bs;
%     RMS_data(ifil) = rms(dtwt(:,ifil));
    lon(ifil) = synth.datamat.loc_lolaz(1);
    lat(ifil) = synth.datamat.loc_lolaz(2);
    x_sta(ifil) = synth.datamat.loc_xyz(1);
    y_sta(ifil) = synth.datamat.loc_xyz(2);
    if ~(ifil == iSIO(1) || ifil == iSIO(2))
        misfit_TAT_bs(:,ifil) = synth.datamat.TAT_bs - trudata.tat;
        misfit_TAT(:,ifil) = rms(misfit_TAT_bs(:,ifil));
        misfit_Vp_bs(:,ifil) = synth.datamat.V_w_bs - trudata.vp_actual*1000;
        misfit_Vp(:,ifil) = rms(misfit_Vp_bs(:,ifil));
    else
        misfit_TAT(:,ifil) = rms(0.013 - trudata.tat);
        misfit_Vp(:,ifil) = rms(1500 - trudata.vp_actual*1000);
        misfit_Vp_bs(:,ifil) = misfit_Vp(:,ifil);
    end
    misfit_r_xy_bs(:,ifil) = sqrt( misfit_xsta_bs(:,ifil).^2 + misfit_ysta_bs(:,ifil).^2 );
%     misfit_r_xy_bs(:,ifil) = sqrt( (misfit_xsta_bs(:,ifil)./trudata.obs_location_xyz(1)*1000).^2 +...
%                                    (misfit_ysta_bs(:,ifil)./trudata.obs_location_xyz(2)*1000).^2 +...
%                                    (misfit_zsta_bs(:,ifil)./(-trudata.obs_location_xyz(3)*1000)).^2 +...
%                                    (misfit_TAT_bs(:,ifil)./trudata.tat).^2 +...
%                                    (misfit_Vp_bs(:,ifil)./trudata.vp_actual).^2);
    misfit_r_xy(ifil) = rms(misfit_r_xy_bs(:,ifil));    
    r_xy_95(:,ifil) = abs(median(misfit_r_xy_bs(:,ifil)) - prctile(misfit_r_xy_bs(:,ifil),[2.5 97.5]));
    
    if ifil == 1
        Ftest_res = synth.datamat.Ftest_res;
    end
    
    fprintf('%s : diff r_xy = %.2f m\n',xlabels{ifil},sqrt(misfit_xsta(ifil).^2+misfit_ysta(ifil).^2))
    fprintf('               = %.2f m (+/- %.2f)\n',mean(misfit_r_xy_bs(:,ifil)),std(misfit_r_xy_bs(:,ifil)))
    fprintf('   : diff z    = %.2f m (+/- %.2f)\n',mean(misfit_zsta_bs(:,ifil)),std(misfit_zsta_bs(:,ifil)))
    fprintf('   : diff Vp   = %.2f m (+/- %.2f)\n',mean(misfit_Vp_bs(:,ifil)),std(misfit_Vp_bs(:,ifil)))
end
r_xy_95(:,8) = [nan nan]';
r_xy_95(:,9) = [nan nan]';


%% ---------------------   PLOTTING   ---------------------   

%% Plot comparisons
f905 = figure(905); clf;
set(f905,'position',[118     1   891   704]);
ax1 = subplot(5,4,[1 2]+2); hold(ax1,'on'); box on;
ax2 = subplot(5,4,[5 6]+2); hold(ax2,'on'); box on;
ax3 = subplot(5,4,[9 10]+2); hold(ax3,'on'); box on;
ax4 = subplot(5,4,[13 14]+2); hold(ax4,'on'); box on;
ax5 = subplot(5,4,[17 18]+2); hold(ax5,'on'); box on;
ax6 = subplot(5,4,[3,4,7,8]-2); hold(ax6,'on'); box on;
ax7 = subplot(5,4,[15,16,19,20]-2); hold(ax7,'on'); box on;
dy = 1.5;
dy_space = -0.035; %-0.04;
dx_space = 0.055;
ax1.Position = [ax1.Position(1)+dx_space, ax1.Position(2)+dy_space, ax1.Position(3), ax1.Position(4)*dy];
ax2.Position = [ax2.Position(1)+dx_space, ax2.Position(2)+dy_space*2, ax2.Position(3), ax2.Position(4)*dy];
ax3.Position = [ax3.Position(1)+dx_space, ax3.Position(2)+dy_space*3, ax3.Position(3), ax3.Position(4)*dy];
ax4.Position = [ax4.Position(1)+dx_space, ax4.Position(2)+dy_space*4, ax4.Position(3), ax4.Position(4)*dy];
ax5.Position = [ax5.Position(1)+dx_space, ax5.Position(2)+dy_space*5, ax5.Position(3), ax5.Position(4)*dy];
dyy = 1.4;
dy = 1.25;
dxx_space = -0.01;
dyy_space = -0.08;
ax6.Position = [ax6.Position(1)+dxx_space, ax6.Position(2)+dyy_space, ax6.Position(3), ax6.Position(4)*dyy];
ax7.Position = [ax7.Position(1)+dxx_space, ax7.Position(2)+dyy_space+0.05, ax7.Position(3), ax7.Position(4)*dy];

%clr = parula(Nfils);
% clr = [brewermap(7,'Blues'); brewermap(3,'Greys'); brewermap(2,'Greens'); brewermap(1,'Reds'); brewermap(2,'Purples'); brewermap(2,'RdPu')];
clr = [repmat([0.7 0.7 0.7],length(synth_dirs)-length(iSIO),1); [0.8 0 0]; [0.8 0 0] ];
% for ifil = 1:length(data_dirs)
%     plot(ax1,[ifil ifil],[0.1 RMS_data(ifil)*1000],'-k','linewidth',1.5)
%     h(ifil) = plot(ax1,ifil,RMS_data(ifil)*1000,symbols{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
% %     errorbar(ax1,ifil,RMS_data(ifil)*1000,RMS_95(1,ifil)*1000,RMS_95(2,ifil)*1000,'.k','markerfacecolor',[0.5 0.5 0.5],'markersize',markersize,'linewidth',1.5,'CapSize',13); hold on;
%     set(ax1,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[]);
%     xticks(ax1,[1:9]);
%     ylabel(ax1,'$\delta$TWT (ms)','fontsize',18,'Interpreter','latex')
%     title(ax1,'$\textbf{RMS}$','fontsize',18,'Interpreter','latex');
%     xlim(ax1,[0 length(data_dirs)+1]);   
%     ylim(ax1,[1 max(RMS_data*1000)+10^(floor(log10(max(RMS_data*1000))))]);
% end

for ifil = 1:length(synth_dirs)
    markersize = sizes(ifil);
    
    %%%%%%% TWT %%%%%%%
    plot(ax1,[ifil ifil],[0.1 RMS_data(ifil)*1000],'-k','linewidth',1.5)
    if issolve(ifil,1)
        h(ifil) = plot(ax1,ifil,RMS_data(ifil)*1000,symbols{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    else
        h(ifil) = plot(ax1,ifil,RMS_data(ifil)*1000,symbols{ifil},'markerfacecolor','w','markeredgecolor',clr(ifil,:),'markersize',markersize,'linewidth',2); hold on;
    end
%     errorbar(ax1,ifil,RMS_data(ifil)*1000,RMS_95(1,ifil)*1000,RMS_95(2,ifil)*1000,'.k','markerfacecolor',[0.5 0.5 0.5],'markersize',markersize,'linewidth',1.5,'CapSize',13); hold on;
    set(ax1,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[]);
    xticks(ax1,[1:9]);
    yticks(ax1,[0.1 1 10 100]);
    ylabel(ax1,'$\delta$TWT (ms)','fontsize',18,'Interpreter','latex')
    title(ax1,'$\textbf{RMS Misfit}$','fontsize',18,'Interpreter','latex');
    xlim(ax1,[0 length(synth_dirs)+1]);   
    ylim(ax1,[1 max(RMS_data*1000)+10^(floor(log10(max(RMS_data*1000))))*2]);
    
    %%%%%%% r_xy %%%%%%%
    plot(ax2,[ifil ifil],[0.1 misfit_r_xy(ifil)],'-k','linewidth',1.5)
    if issolve(ifil,2)
        h(ifil) = plot(ax2,ifil,misfit_r_xy(ifil),symbols{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    else
        h(ifil) = plot(ax2,ifil,misfit_r_xy(ifil),symbols{ifil},'markerfacecolor','w','markeredgecolor',clr(ifil,:),'markersize',markersize,'linewidth',2); hold on;
    end
%     errorbar(ax2,ifil,misfit_r_xy(ifil),r_xy_95(1,ifil),r_xy_95(2,ifil),'.k','markerfacecolor',[0.5 0.5 0.5],'markersize',markersize,'linewidth',1.5,'CapSize',13); hold on;
    set(ax2,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[]);
    ylabel(ax2,'$\delta r_{xy}$ (m)','fontsize',18,'Interpreter','latex')
    yticks(ax2,[0.1 1 10 100]);
    xticks(ax2,[1:9]);
    xlim(ax2,[0 length(synth_dirs)+1]);   
    ylim(ax2,[1 max(misfit_r_xy)+10^(floor(log10(max(misfit_r_xy))))]);
    
    %%%%%%% Z %%%%%%%
    plot(ax3,[ifil ifil],[0.1 misfit_zsta(ifil)],'-k','linewidth',1.5)
    if issolve(ifil,3)
        h(ifil) = plot(ax3,ifil,misfit_zsta(ifil),symbols{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    else
        h(ifil) = plot(ax3,ifil,misfit_zsta(ifil),symbols{ifil},'markerfacecolor','w','markeredgecolor',clr(ifil,:),'markersize',markersize,'linewidth',2); hold on;
    end
%     errorbar(ax3,ifil,misfit_r_xy(ifil),r_xy_95(1,ifil),r_xy_95(2,ifil),'.k','markerfacecolor',[0.5 0.5 0.5],'markersize',markersize,'linewidth',1.5,'CapSize',13); hold on;
    set(ax3,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[]);
    ylabel(ax3,'${\delta Z}$ (m)','fontsize',18,'Interpreter','latex')
    yticks(ax3,[0.1 1 10 100]);
    xticks(ax3,[1:9]);
    xlim(ax3,[0 length(synth_dirs)+1]);   
    ylim(ax3,[1 max(misfit_zsta)+10^(floor(log10(max(misfit_zsta))))*5]);
    
    %%%%%%% TAT %%%%%%%
    plot(ax5,[ifil ifil],[0 misfit_TAT(ifil)*1000],'-k','linewidth',1.5)
    if issolve(ifil,4)
        h(ifil) = plot(ax5,ifil,misfit_TAT(ifil)*1000,symbols{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    else
        h(ifil) = plot(ax5,ifil,misfit_TAT(ifil)*1000,symbols{ifil},'markerfacecolor','w','markeredgecolor',clr(ifil,:),'markersize',markersize,'linewidth',2); hold on;
    end
%     errorbar(ax2,ifil,misfit_r_xy(ifil),r_xy_95(1,ifil),r_xy_95(2,ifil),'.k','markerfacecolor',[0.5 0.5 0.5],'markersize',markersize,'linewidth',1.5,'CapSize',13); hold on;
    set(ax5,'yscale','linear','linewidth',1.5,'fontsize',16,'XTickLabel',[]);
    ylabel(ax5,'{$\delta\tau$} (ms)','fontsize',18,'Interpreter','latex')
%         yticks(ax4,[0.1 1 10 100]);
    xticks(ax5,[1:9]);
    xlim(ax5,[0 length(synth_dirs)+1]);   
%         ylim(ax4,[1 max(misfit_TAT*1000)+10^(floor(log10(max(misfit_TAT*1000))))]);
    ylim(ax5,[0 1.3]);
    
    %%%%%%% Vp %%%%%%%
    plot(ax4,[ifil ifil],[0.1 misfit_Vp(ifil)],'-k','linewidth',1.5)
    if issolve(ifil,5)
        h(ifil) = plot(ax4,ifil,misfit_Vp(ifil),symbols{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    else
        h(ifil) = plot(ax4,ifil,misfit_Vp(ifil),symbols{ifil},'markerfacecolor','w','markeredgecolor',clr(ifil,:),'markersize',markersize,'linewidth',2); hold on;
    end
%     errorbar(ax5,ifil,misfit_r_xy(ifil),r_xy_95(1,ifil),r_xy_95(2,ifil),'.k','markerfacecolor',[0.5 0.5 0.5],'markersize',markersize,'linewidth',1.5,'CapSize',13); hold on;
    set(ax4,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[]);
    ylabel(ax4,'${\delta V_{P}}$ (m/s)','fontsize',18,'Interpreter','latex')
%         yticks(ax5,[0.1 1 10 100]);
    xticks(ax4,[1:9]);
    xlim(ax4,[0 length(synth_dirs)+1]);   
    ylim(ax4,[1 max(misfit_Vp)+10^(floor(log10(max(misfit_Vp))))]);
    htxt = text(ax4,ifil,0.8,xlabels(ifil),'HorizontalAlignment','right','Interpreter','latex','FontSize',16);
    set(htxt,'Rotation',45);
%     xticklabels(ax2,xlabels);
%     xtickangle(ax2,45);
end
%% F-test plots
trucol = [0.2 0.9 0.6];
col = colormap(parula); % colormap for f test plots
x_grid = Ftest_res.x_grid;
y_grid = Ftest_res.y_grid;
z_grid = Ftest_res.z_grid;
P = Ftest_res.Pstat;
[Xgrd,Ygrd,Zgrd] = meshgrid(x_grid,y_grid,z_grid);
[Pz_max, Iz_max] = max(max(max(P)));
[Py_max, Iy_max] = max(max(P(:,:,Iz_max)));
[Px_max, Ix_max] = max(P(:,Iy_max,Iz_max));


% bounds and ticks (zoom in within test space to fill plot)
x_bds = [min(x_grid),max(x_grid)] + [25 -25];
y_bds = [min(y_grid),max(y_grid)] + [25 -25];
z_bds = [min(z_grid),max(z_grid)];

% centre points:
x_mid = 0; %median(x_grid);
y_mid = 0; %median(y_grid);
z_mid = 0; %median(z_grid);

% bounds and ticks
xlim1 = [-20 20]; dx1 = 4;
ylim1 = [-20 20]; dy1 = 4;
xlim2 = [-40 40]; dx2 = 8;
ylim2 = [-40 40]; dy2 = 8;
zlim = [-40 40];  dz = 8;

% shade in background
% fill(ax6,[-2000,2000,2000,-2000,-2000],[-2000,-2000,2000,2000,-2000],col(1,:))
ang = [0:0.1:2*pi];
%X-Y
% contourf(ax6,Xgrd(:,:,Iz_max)-x_mid,Ygrd(:,:,Iz_max)-y_mid,P(:,:,Iz_max),'linestyle','none');
% shading(ax1,'flat');
contour(ax6,Xgrd(:,:,Iz_max)-x_mid,Ygrd(:,:,Iz_max)-y_mid,P(:,:,Iz_max),[[0.05 0.05],[0.32 0.32]],'-k','linewidth',3); hold on;
h6(1) = plot(ax6,Xgrd(Ix_max,Iy_max,Iz_max)-x_mid,Ygrd(Ix_max,Iy_max,Iz_max)-y_mid,'pk','markerfacecolor',clr(1,:),'markersize',21,'linewidth',1);
set(ax6,'fontsize',16,'linewidth',1.5,'box','on','layer','top',...
    'xlim',xlim1,'ylim',ylim1,'xtick',[xlim1(1):dx1:xlim1(2)],'ytick',[ylim1(1):dy1:ylim1(2)],'TickDir','out');
xlabel(ax6,'$\delta$X (m)','fontsize',18,'interpreter','latex');
ylabel(ax6,'$\delta$Y (m)','fontsize',18,'interpreter','latex');
% title(ax6,'\textbf{X-Y}','fontsize',18,'interpreter','latex');
% average values
htx1 = text(ax6,xlim1(1)+0.04*diff(axlim(ax1,1:2)),ylim1(1)+0.02*diff(axlim(ax1,3:4)),...
    sprintf('$\\mathbf{\\bar{x}}$\\textbf{ = %.1f m}',x_mid),'color','white','interpreter','latex','fontsize',17);
hty1 = text(ax6,xlim1(1)+0.04*diff(axlim(ax1,1:2)),ylim1(1)+0.01*diff(axlim(ax1,3:4)),...
    sprintf('$\\mathbf{\\bar{y}}$\\textbf{ = %.1f m}',y_mid),'color','white','interpreter','latex','fontsize',17);
% axis equal;

% Add on the true location for comparison
% plot the true location on the F-test plot in each dimension
h6(3) = plot(ax6,trudata.obs_location_xyz(1)*1e3 - x_mid,trudata.obs_location_xyz(2)*1e3 - y_mid,'p',...
    'markerfacecolor',trucol,'markeredgecolor','k','linewidth',1,'markersize',20);

% Plot SIO result
h6(2) = plot(ax6,x_sta(iSIO(1)) - x_mid,y_sta(iSIO(1)) - y_mid,'pk','MarkerFaceColor',clr(iSIO(1),:),'markersize',21,'linewidth',1);

% expand bounds and ticks if needed to include SIO station
    set(ax6,'xlim',[min([x_grid,x_sta(iSIO(1))-5]),max([x_grid,x_sta(iSIO(1))+5])],...
            'ylim',[min([y_grid,y_sta(iSIO(1))-5]),max([y_grid,y_sta(iSIO(1))+5])])
    xlim(ax6,[160 220]);
    ylim(ax6,[-420 -360]);
    dtick = [1,2,4,5,10,20,40,50,100];
    dtickx = dtick(mindex(abs((axlim(ax6,2)-axlim(ax6,1))/7 - dtick)));
    dticky = dtick(mindex(abs((axlim(ax6,4)-axlim(ax6,3))/7 - dtick)));
    set(ax6,'xtick',[-1000:dtickx:1000],'xticklabel',[-1000:dtickx:1000],...
            'ytick',[-1000:dticky:1000],'yticklabel',[-1000:dticky:1000]);
    % kill text
    delete([htx1,hty1])
    % plot range
    x12 = [x_sta(iSIO(1)),trudata.obs_location_xyz(1)*1e3];
    y12 = [y_sta(iSIO(1)),trudata.obs_location_xyz(2)*1e3];
    %dxy
    plot(ax6,x12,y12,'-','Color',clr(8,:),'linewidth',2);
    text(ax6,mean(x12),mean(y12),num2str(sqrt(diff(x12).^2 + diff(y12).^2),'%.1f m'),...
        'fontsize',14,'color',clr(8,:),'fontweight','bold',...
        'rotation',atand(vertexag(ax6)*diff(y12)/diff(x12)),...
        'horizontalalignment','center','verticalalignment','top')
    
    % plot range
    x12 = [x_sta(1),trudata.obs_location_xyz(1)*1e3];
    y12 = [y_sta(1),trudata.obs_location_xyz(2)*1e3];
    %dxy
    plot(ax6,x12,y12,'-','Color',clr(1,:),'linewidth',2);
    text(ax6,mean(x12-2),mean(y12-2),num2str(sqrt(diff(x12).^2 + diff(y12).^2),'%.1f m'),...
        'fontsize',14,'color',[0.4 0.4 0.4],'fontweight','bold',...
        'rotation',atand(vertexag(ax6)*diff(y12)/diff(x12)),...
        'horizontalalignment','center','verticalalignment','top')
    
    legend(ax6,h6,{xlabels{1},xlabels{iSIO(1)},'True'},'fontsize',15,'location','southwest','box','off');

%% RESIDUALS
plot(ax7,[0 360],[0 0],'--k','linewidth',1.5);
h1(1) = plot(ax7,trudata.survaz(~isnan(trudata.data.tot_dt)),dtwt(:,1)*1000,'ok','markersize',10,'MarkerFaceColor',clr(1,:),'linewidth',1); hold on;
h1(2) = plot(ax7,trudata.survaz(~isnan(trudata.data.tot_dt)),dtwt(:,iSIO(1))*1000,'ok','markersize',10,'MarkerFaceColor',clr(iSIO(1),:),'linewidth',1); hold on;
% plot(ax7,trudata.survaz(~isnan(trudata.data.tot_dt)),dtwt(:,1)*1000-dtwt(:,3)*1000,'ok','markersize',10,'MarkerFaceColor',[0 0 1],'linewidth',1); hold on;
% rms(dtwt(:,1)*1000-dtwt(:,3)*1000)
% max(dtwt(:,1)*1000-dtwt(:,3)*1000)-min(dtwt(:,1)*1000-dtwt(:,3)*1000)
xlim(ax7,[0 360]);
ylim(ax7,[-50 20]);
xlabel(ax7,'Ship azimuth ($^\circ$)','fontsize',18,'interpreter','latex');
ylabel(ax7,'$\delta$TWT (ms)','fontsize',18,'interpreter','latex');
set(ax7,'fontsize',16,'box','on','linewidth',1.5)
legend(ax7,h1,{xlabels{1},xlabels{iSIO(1)}},'fontsize',15,'location','southeast','box','off');


% l = legend(h,lgd,'position',[0.8 0.05 0.1852 0.95],'interpreter','none','fontsize',14,'box','off');

%% Figure numbers
x = 0.015;
y = 0.82;
text(ax1,x,y,...
'\textbf{c)}','color',[0 0 0],'interpreter','latex','fontsize',20,'Units','normalized');
text(ax2,x,y,...
'\textbf{d)}','color',[0 0 0],'interpreter','latex','fontsize',20,'Units','normalized');
text(ax3,x,y,...
'\textbf{e)}','color',[0 0 0],'interpreter','latex','fontsize',20,'Units','normalized');
text(ax4,x,y,...
'\textbf{f)}','color',[0 0 0],'interpreter','latex','fontsize',20,'Units','normalized');
text(ax5,x,y,...
'\textbf{g)}','color',[0 0 0],'interpreter','latex','fontsize',20,'Units','normalized');

text(ax6,x,y+0.12,...
'\textbf{a)}','color',[0 0 0],'interpreter','latex','fontsize',20,'Units','normalized');
text(ax7,x,y+0.11,...
'\textbf{b)}','color',[0 0 0],'interpreter','latex','fontsize',20,'Units','normalized');

%%
delete(ax5);
%% Plot histograms
% figure(1);
% 
% for ii = 1:2
%     plot_hist(subplot(1,2,ii),misfit_r_xy_bs(:,ii),10)
%     xlabel('\delta r_{xy}');
% end

%% SAVE
if ifsave
    save2pdf(ofile,f905)
end

%% Statistics for Doppler improvements
RMS_improve = (RMS_data(2)-RMS_data(1)) / RMS_data(2)*100;
misfit_r_xy_improve = (misfit_r_xy(2)-misfit_r_xy(1)) / misfit_r_xy(2)*100;
fprintf('\n** DOPPLER IMPROVEMENTS **\n')
fprintf(' dTWT improve: %.3f ms (%.3f%%)\n',(RMS_data(2)-RMS_data(1))*1000,RMS_improve)
fprintf('dr_xy improve: %.3f m (%.3f%%)\n',(misfit_r_xy(2)-misfit_r_xy(1)),misfit_r_xy_improve)


end
