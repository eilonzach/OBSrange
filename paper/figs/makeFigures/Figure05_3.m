%% Function to produce Figure 05 
%  Figure 05 compares our algorithm to SIO
function Figure05

ofile = '../Figure05';
ifsave = 0;

%% load 
data_dirs = {
    '2_OUT_wcorr_xrec';
    '1_OUT_nocorr';
    '7_OUT_nocorr_noellipsoid';
    '3_OUT_nocorr_TAT';
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
    '3_OUT_wcorr_xrec_TAT';
    '4_OUT_wcorr_xrec_Vp';
    '5_OUT_wcorr_xrec_Z';
    '6_OUT_wcorr_xrec_TAT_Vp_Z';
    '8_SIO_compare_nobads';
    '9_SIO_compare_wbads';
    };

xlabels = {
    'OBSrange';
    'No Doppler';
    'No Ellipsoid';
    'X, Y, Z, $V_p$';
    'X, Y, Z, $\tau$';
    'X, Y, $\tau$, $V_p$';
    'X, Y';
    'SIO';
    'SIO (+ bads)';
    };

data_path = '../figdata/PacificORCA_EC03/OUT_OBSrange';
synth_path = '../figdata/PacificORCA_synthtest3/OUT_OBSrange';

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

% Load synthetic
trudata = load('../figdata/PacificORCA_synthtest3/trudata_syn12.mat');
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
    lon(ifil) = synth.datamat.loc_lolaz(1);
    lat(ifil) = synth.datamat.loc_lolaz(2);
    if ifil ~= 8 && ifil ~= 9
        misfit_TAT_bs(:,ifil) = synth.datamat.TAT_bs - trudata.tat;
        misfit_TAT(:,ifil) = rms(misfit_TAT_bs(:,ifil));
        misfit_Vp_bs(:,ifil) = synth.datamat.V_w_bs - trudata.vp_actual*1000;
        misfit_Vp(:,ifil) = rms(misfit_Vp_bs(:,ifil));
    end
    misfit_r_xy_bs(:,ifil) = sqrt( misfit_xsta_bs(:,ifil).^2 + misfit_ysta_bs(:,ifil).^2 );
%     misfit_r_xy_bs(:,ifil) = sqrt( (misfit_xsta_bs(:,ifil)./trudata.obs_location_xyz(1)*1000).^2 +...
%                                    (misfit_ysta_bs(:,ifil)./trudata.obs_location_xyz(2)*1000).^2 +...
%                                    (misfit_zsta_bs(:,ifil)./(-trudata.obs_location_xyz(3)*1000)).^2 +...
%                                    (misfit_TAT_bs(:,ifil)./trudata.tat).^2 +...
%                                    (misfit_Vp_bs(:,ifil)./trudata.vp_actual).^2);
    misfit_r_xy(ifil) = rms(misfit_r_xy_bs(:,ifil));    
    r_xy_95(:,ifil) = abs(median(misfit_r_xy_bs(:,ifil)) - prctile(misfit_r_xy_bs(:,ifil),[2.5 97.5]));
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
dy = 1.25;
dy_space = 1.6;
dx_space = 0.055;
ax1.Position = [ax1.Position(1)+dx_space, ax1.Position(2), ax1.Position(3), ax1.Position(4)*dy];
ax2.Position = [ax2.Position(1)+dx_space, ax2.Position(2), ax2.Position(3), ax2.Position(4)*dy];
ax3.Position = [ax3.Position(1)+dx_space, ax3.Position(2), ax3.Position(3), ax3.Position(4)*dy];
ax4.Position = [ax4.Position(1)+dx_space, ax4.Position(2), ax4.Position(3), ax4.Position(4)*dy];
ax5.Position = [ax5.Position(1)+dx_space, ax5.Position(2), ax5.Position(3), ax5.Position(4)*dy];
dyy = 1.4;
dxx_space = -0.01;
dyy_space = -0.08;
ax6.Position = [ax6.Position(1)+dxx_space, ax6.Position(2)+dyy_space, ax6.Position(3), ax6.Position(4)*dyy];
ax7.Position = [ax7.Position(1)+dxx_space, ax7.Position(2)+dyy_space+0.05, ax7.Position(3), ax7.Position(4)*dy];

markersize = 14;
%clr = parula(Nfils);
% clr = [brewermap(7,'Blues'); brewermap(3,'Greys'); brewermap(2,'Greens'); brewermap(1,'Reds'); brewermap(2,'Purples'); brewermap(2,'RdPu')];
clr = [repmat([0.7 0.7 0.7],7,1); [0.8 0 0]; [0.8 0 0] ];
for ifil = 1:length(data_dirs)
    plot(ax1,[ifil ifil],[0.1 RMS_data(ifil)*1000],'-k','linewidth',1.5)
    h(ifil) = plot(ax1,ifil,RMS_data(ifil)*1000,'ok','markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
%     errorbar(ax1,ifil,RMS_data(ifil)*1000,RMS_95(1,ifil)*1000,RMS_95(2,ifil)*1000,'.k','markerfacecolor',[0.5 0.5 0.5],'markersize',markersize,'linewidth',1.5,'CapSize',13); hold on;
    set(ax1,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[]);
    xticks(ax1,[1:9]);
    ylabel(ax1,'$\delta$TWT (ms)','fontsize',18,'Interpreter','latex')
    title(ax1,'$\textbf{RMS}$','fontsize',18,'Interpreter','latex');
    xlim(ax1,[0 length(data_dirs)+1]);   
    ylim(ax1,[1 max(RMS_data*1000)+10^(floor(log10(max(RMS_data*1000))))]);
end

for ifil = 1:length(synth_dirs)
    plot(ax2,[ifil ifil],[0.1 misfit_r_xy(ifil)],'-k','linewidth',1.5)
    h(ifil) = plot(ax2,ifil,misfit_r_xy(ifil),'ok','markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
%     errorbar(ax2,ifil,misfit_r_xy(ifil),r_xy_95(1,ifil),r_xy_95(2,ifil),'.k','markerfacecolor',[0.5 0.5 0.5],'markersize',markersize,'linewidth',1.5,'CapSize',13); hold on;
    set(ax2,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[]);
    ylabel(ax2,'$\delta r_{xy}$ (m)','fontsize',18,'Interpreter','latex')
    yticks(ax2,[0.1 1 10 100]);
    xticks(ax2,[1:9]);
    xlim(ax2,[0 length(synth_dirs)+1]);   
    ylim(ax2,[1 max(misfit_r_xy)+10^(floor(log10(max(misfit_r_xy))))]);
    
    plot(ax3,[ifil ifil],[0.1 misfit_zsta(ifil)],'-k','linewidth',1.5)
    h(ifil) = plot(ax3,ifil,misfit_zsta(ifil),'ok','markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
%     errorbar(ax3,ifil,misfit_r_xy(ifil),r_xy_95(1,ifil),r_xy_95(2,ifil),'.k','markerfacecolor',[0.5 0.5 0.5],'markersize',markersize,'linewidth',1.5,'CapSize',13); hold on;
    set(ax3,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[]);
    ylabel(ax3,'${\delta Z}$ (m)','fontsize',18,'Interpreter','latex')
    yticks(ax3,[0.1 1 10 100]);
    xticks(ax3,[1:9]);
    xlim(ax3,[0 length(synth_dirs)+1]);   
    ylim(ax3,[10 max(misfit_zsta)+10^(floor(log10(max(misfit_zsta))))]);
    
    if ifil ~= 8 && ifil ~= 9
        plot(ax4,[ifil ifil],[0 misfit_TAT(ifil)*1000],'-k','linewidth',1.5)
        h(ifil) = plot(ax4,ifil,misfit_TAT(ifil)*1000,'ok','markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    %     errorbar(ax2,ifil,misfit_r_xy(ifil),r_xy_95(1,ifil),r_xy_95(2,ifil),'.k','markerfacecolor',[0.5 0.5 0.5],'markersize',markersize,'linewidth',1.5,'CapSize',13); hold on;
        set(ax4,'yscale','linear','linewidth',1.5,'fontsize',16,'XTickLabel',[]);
        ylabel(ax4,'{$\delta\tau$} (ms)','fontsize',18,'Interpreter','latex')
%         yticks(ax4,[0.1 1 10 100]);
        xticks(ax4,[1:9]);
        xlim(ax4,[0 length(synth_dirs)+1]);   
%         ylim(ax4,[1 max(misfit_TAT*1000)+10^(floor(log10(max(misfit_TAT*1000))))]);
        ylim(ax4,[0 1.2]);

        plot(ax5,[ifil ifil],[0.1 misfit_Vp(ifil)],'-k','linewidth',1.5)
        h(ifil) = plot(ax5,ifil,misfit_Vp(ifil),'ok','markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    %     errorbar(ax5,ifil,misfit_r_xy(ifil),r_xy_95(1,ifil),r_xy_95(2,ifil),'.k','markerfacecolor',[0.5 0.5 0.5],'markersize',markersize,'linewidth',1.5,'CapSize',13); hold on;
        set(ax5,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[]);
        ylabel(ax5,'${\delta V_{P}}$ (m/s)','fontsize',18,'Interpreter','latex')
%         yticks(ax5,[0.1 1 10 100]);
        xticks(ax5,[1:9]);
        xlim(ax5,[0 length(synth_dirs)+1]);   
        ylim(ax5,[1 max(misfit_Vp)+10^(floor(log10(max(misfit_Vp))))]);
    end
    htxt = text(ax5,ifil,0.8,xlabels(ifil),'HorizontalAlignment','right','Interpreter','latex','FontSize',16);
    set(htxt,'Rotation',45);
%     xticklabels(ax2,xlabels);
%     xtickangle(ax2,45);
end

% SURVEY PATTERN
latlim = [min(trudata.survlat) max(trudata.survlat)]+[-1 1]*.2/111; % 200m outside ship circle
lonlim = [min(trudata.survlon) max(trudata.survlon)]+[-1 1]*.2/111; % 200m outside ship circle
plot(ax6,trudata.survlon(~isnan(trudata.data.tot_dt)),trudata.survlat(~isnan(trudata.data.tot_dt)),'ok','markersize',10,'MarkerFaceColor',[0 0 0],'linewidth',1); hold on;
plot(ax6,trudata.drop_location(2),trudata.drop_location(1),'sk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1);
plot(ax6,trudata.obs_location_laloz(2),trudata.obs_location_laloz(1),'pk','MarkerFaceColor',[1 1 0],'markersize',21,'linewidth',1);
plot(ax6,lon(8),lat(8),'p','MarkerEdgeColor',clr(8,:),'markersize',21,'linewidth',2);
plot(ax6,lon(1),lat(1),'p','MarkerEdgeColor',clr(1,:),'markersize',21,'linewidth',2);
axis(ax6,'equal');
xlabel(ax6,'Longitude','fontsize',18,'interpreter','latex');
ylabel(ax6,'Latitude','fontsize',18,'interpreter','latex');
grid(ax6,'on');
set(ax6,'fontsize',16,'box','on','linewidth',1.5,...
        'xlim',lonlim,'ylim',latlim);

% RESIDUALS
plot(ax7,[0 360],[0 0],'--k','linewidth',1.5);
h1(1) = plot(ax7,trudata.survaz(~isnan(trudata.data.tot_dt)),dtwt(:,1)*1000,'ok','markersize',10,'MarkerFaceColor',clr(1,:),'linewidth',1); hold on;
h1(2) = plot(ax7,trudata.survaz(~isnan(trudata.data.tot_dt)),dtwt(:,8)*1000,'ok','markersize',10,'MarkerFaceColor',clr(8,:),'linewidth',1); hold on;
xlim(ax7,[0 360]);
ylim(ax7,[-50 20]);
xlabel(ax7,'Ship azimuth ($^\circ$)','fontsize',18,'interpreter','latex');
ylabel(ax7,'$\delta$TWT (ms)','fontsize',18,'interpreter','latex');
set(ax7,'fontsize',16,'box','on','linewidth',1.5)
legend(ax7,h1,{xlabels{1},xlabels{8}},'fontsize',16,'location','southeast','box','off');


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

text(ax6,x,y+0.13,...
'\textbf{a)}','color',[0 0 0],'interpreter','latex','fontsize',20,'Units','normalized');
text(ax7,x,y+0.13,...
'\textbf{b)}','color',[0 0 0],'interpreter','latex','fontsize',20,'Units','normalized');


%% Plot histograms
figure(1);

for ii = 1:2
    plot_hist(subplot(1,2,ii),misfit_r_xy_bs(:,ii),10)
    xlabel('\delta r_{xy}');
end

%% SAVE
if ifsave
    save2pdf(ofile,f905)
end

end
