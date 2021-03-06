%% Function to produce Figure 06 
%  Figure 06 compare survey pattern geometries
function Figure06

ofile = '../Figure06_mean';
ifsave = 1;

projname = 'mats_SynthBoot_summary';

x_rms = 100;
y_rms = 100;
rxy_rms = sqrt(x_rms^2+y_rms^2);
z_rms = 50;
TAT_rms = 3; % ms
Vp_rms = 10;

%% Load *.mat files
files = dir(['../figdata/',projname,'/*.mat']);

Nfils = length(files);
for ifil = 1:Nfils
    load(['../figdata/',projname,'/',files(ifil).name]);
    
    data_all(ifil).x_ship = data_summary.x_ship;
    data_all(ifil).y_ship = data_summary.y_ship;
    data_all(ifil).v_surv_true = data_summary.v_surv_true;
    data_all(ifil).vmag = data_summary.vmag;
    data_all(ifil).survx = data_summary.survx;
    data_all(ifil).survy = data_summary.survy;
    data_all(ifil).survey = data_summary.survey;
    data_all(ifil).radius = data_summary.radius;
    
    
    misfit_xsta(ifil,1) = data_summary.misfit_xsta;
    misfit_xsta_std(ifil,1) = data_summary.misfit_xsta_std;
    misfit_xsta_mean(ifil,1) = abs(data_summary.misfit_xsta_mean);
    misfit_ysta(ifil,1) = data_summary.misfit_ysta;
    misfit_ysta_std(ifil,1) = data_summary.misfit_ysta_std;
    misfit_ysta_mean(ifil,1) = abs(data_summary.misfit_ysta_mean);
    misfit_zsta(ifil,1) = data_summary.misfit_zsta;
    misfit_zsta_std(ifil,1) = data_summary.misfit_zsta_std;
    misfit_zsta_mean(ifil,1) = abs(data_summary.misfit_zsta_mean);
    misfit_zsta_med(ifil,1) = abs(data_summary.misfit_zsta_med);
    misfit_r_xy(ifil,1) = data_summary.misfit_r_xy;
    misfit_r_xy_std(ifil,1) = data_summary.misfit_r_xy_std;
    misfit_r_xy_mean(ifil,1) = sqrt(data_summary.misfit_xsta_mean.^2 + data_summary.misfit_ysta_mean.^2); %abs(data_summary.misfit_r_xy_mean);
    misfit_r_xy_med(ifil,1) = abs(data_summary.misfit_r_xy_med);
    misfit_r_xyz(ifil,1) = data_summary.misfit_r_xyz;
    misfit_r_xyz_std(ifil,1) = data_summary.misfit_r_xyz_std;
    misfit_TAT(ifil,1) = data_summary.misfit_TAT;
    misfit_TAT_std(ifil,1) = data_summary.misfit_TAT_std;
    misfit_TAT_mean(ifil,1) = abs(data_summary.misfit_TAT_mean);
    misfit_TAT_med(ifil,1) = abs(data_summary.misfit_TAT_med);
    misfit_Vw(ifil,1) = data_summary.misfit_Vw;
    misfit_Vw_std(ifil,1) = data_summary.misfit_Vw_std;
    misfit_Vw_mean(ifil,1) = abs(data_summary.misfit_Vw_mean);
    misfit_Vw_med(ifil,1) = abs(data_summary.misfit_Vw_med);
    E_rms(ifil,1) = data_summary.E_rms;
    E_rms_std(ifil,1) = data_summary.E_rms_std;
    misfit_v_ship_all(1:2,ifil) = data_summary.misfit_v_ship_all;
    misfit_v_ship_all_std(1:2,ifil) = data_summary.misfit_v_ship_all_std;
    misfit_dtwtcorr_all(ifil,1) = data_summary.misfit_dtwtcorr_all;
    misfit_dtwtcorr_all_std(ifil,1) = data_summary.misfit_dtwtcorr_all_std;
    dtwt_all(ifil,1) = data_summary.dtwt_all;
    dtwt_all_std(ifil,1) = data_summary.dtwt_all_std;
    
    lgd{ifil} = [data_all(ifil).survey,' ',num2str(data_all(ifil).radius),' Nm']; %files(ifil).name(1:end-4);
    if any(regexp(lgd{ifil},'PACMAN'))
        symbol{ifil} = 'ok';
    elseif any(regexp(lgd{ifil},'cross'))
        symbol{ifil} = 'sk';
    elseif any(regexp(lgd{ifil},'diamond'))
        symbol{ifil} = 'dk';
    elseif any(regexp(lgd{ifil},'line'))
        symbol{ifil} = 'pk';
    elseif any(regexp(lgd{ifil},'tri'))
        symbol{ifil} = '^k';
    elseif any(regexp(lgd{ifil},'circle'))
        symbol{ifil} = 'ok';
    end
    
end

%% Plots
% Model Parameters
f906 = figure(906); clf;
set(gcf,'position',[1     1   918   697])
ax1 = subplot(4,4,[1 3]); hold(ax1,'on'); box on;
ax2 = subplot(4,4,[5 7]); hold(ax2,'on'); box on;
ax3 = subplot(4,4,[9 11]); hold(ax3,'on'); box on;
ax4 = subplot(4,4,[13 15]); hold(ax4,'on'); box on;
dyy = 1.3;
dxx = 1.05;
dx = -0.05;
dy = -0.05;
yshift = 0.01;
ax1.Position = [ax1.Position(1)+dx, ax1.Position(2)+dy+yshift*4, ax1.Position(3)/2*dxx, ax1.Position(4)*dyy];
ax2.Position = [ax2.Position(1)+dx, ax2.Position(2)+dy+yshift*3, ax2.Position(3)/2*dxx, ax2.Position(4)*dyy];
ax3.Position = [ax3.Position(1)+dx, ax3.Position(2)+dy+yshift*2, ax3.Position(3)/2*dxx, ax3.Position(4)*dyy];
ax4.Position = [ax4.Position(1)+dx, ax4.Position(2)+dy+yshift*1, ax4.Position(3)/2*dxx, ax4.Position(4)*dyy];

ax5 = subplot(4,4,[4]); hold(ax5,'on'); box on;
ax6 = subplot(4,4,[8]); hold(ax6,'on'); box on;
ax7 = subplot(4,4,[12]); hold(ax7,'on'); box on;
ax8 = subplot(4,4,[16]); hold(ax8,'on'); box on;
dxcol2 = ax1.Position(1)+ax1.Position(3);
dxspace = -0.03;
ax5.Position = [ax1.Position(1)+dxcol2+dxspace, ax1.Position(2), ax1.Position(3), ax1.Position(4)];
ax6.Position = [ax2.Position(1)+dxcol2+dxspace, ax2.Position(2), ax2.Position(3), ax2.Position(4)];
ax7.Position = [ax3.Position(1)+dxcol2+dxspace, ax3.Position(2), ax3.Position(3), ax3.Position(4)];
ax8.Position = [ax4.Position(1)+dxcol2+dxspace, ax4.Position(2), ax4.Position(3), ax4.Position(4)];

markersize = 12;
%clr = parula(Nfils);
purp = brewermap(3,'Purples');
Nlines = 2;
purp2 = brewermap(4,'Purples');
clr = [brewermap(7,'Blues'); brewermap(3,'Greys'); brewermap(2,'Greens'); brewermap(1,'Reds'); purp(end-1:end,:); brewermap(2,'RdPu')];
iline = 0;
for ifil = 1:Nfils
    %% MEAN
    h(ifil) = plot(ax1,ifil,misfit_r_xy_mean(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    if any(regexp(lgd{ifil},'line'))
        iline = iline + 1;
        plot(ax1,ifil,misfit_xsta_mean(ifil),symbol{ifil},'markeredgecolor',purp2(end-(Nlines-iline),:),'markersize',markersize,'linewidth',2); hold on;
    end
    set(ax1,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3);
    ylabel(ax1,'$\mathbf{\delta r_{xy}\, (m)}$','fontsize',18,'Interpreter','latex')
    xlim(ax1,[0 Nfils+1]);   
    title(ax1,'\textbf{Average}','interpreter','latex','fontsize',18);
    ylim(ax1,[0.01 max(misfit_r_xy_med)+10^(floor(log10(max(misfit_r_xy_med))))*7]);
    yticks(ax1,[0.001 0.01 0.1 1 10 100]);

    plot(ax2,ifil,misfit_zsta_mean(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax2,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3);
    ylabel(ax2,'$\mathbf{\delta Z\, (m)}$','fontsize',18,'Interpreter','latex')
    xlim(ax2,[0 Nfils+1]);
    ylim(ax2,[0.05 max(misfit_zsta_mean)+10^(floor(log10(max(misfit_zsta_mean))))]);
    yticks(ax2,[0.1 1 10 100]);
    
    plot(ax3,ifil,misfit_TAT_mean(ifil)*1000,symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax3,'yscale','log','linewidth',1.5,'fontsize',16,'xticklabel',[],'yminortick','on','xtick',[],'TickLength',[0.01, 0.001]*3);
    ylabel(ax3,'{$\delta$\boldmath$\tau$ (\textbf{ms})}','fontsize',18,'Interpreter','latex')
    xlim(ax3,[0 Nfils+1]);
    ylim(ax3,[0.003 0.5]);
%     yticks(ax3,[0:0.1:0.4]);
    yticks(ax3,[0.001, 0.01, 0.1, 1]);

    plot(ax4,ifil,misfit_Vw_mean(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax4,'yscale','log','linewidth',1.5,'fontsize',16,'xticklabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3);
    ylabel(ax4,'$\mathbf{\delta V_{P} \, (m/s)}$','fontsize',18,'Interpreter','latex')
    xlim(ax4,[0 Nfils+1]);
%     ylim(ax4,[1 max(misfit_Vw_mean)+10^(floor(log10(max(misfit_Vw_mean))))]);
    yticks(ax4,[0.001 0.01 0.1 1 10 100]);

    %% RMS
    if ifil==1
        plot(ax5,[0 Nfils+1],[rxy_rms rxy_rms],'-','color',[0.7 0 0],'linewidth',3); hold on;
    end
    h(ifil) = plot(ax5,ifil,misfit_r_xy(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    if any(regexp(lgd{ifil},'line'))
        plot(ax5,ifil,misfit_xsta(ifil),symbol{ifil},'markeredgecolor',purp2(end-(Nlines-iline),:),'markersize',markersize,'linewidth',2); hold on;
    end
    set(ax5,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3);
%     ylabel(ax5,'$\mathbf{\delta r_{xy}\, (m)}$','fontsize',18,'Interpreter','latex')
    xlim(ax5,[0 Nfils+1]);   
    ylim(ax5,[1 max(misfit_r_xy)+10^(floor(log10(max(misfit_r_xy))))*10]);
    title(ax5,'\textbf{RMS}','interpreter','latex','fontsize',18);
    yticks(ax5,[0.001 0.01 0.1 1 10 100 1000]);
    
    if ifil==1
        plot(ax6,[0 Nfils+1],[z_rms z_rms],'-','color',[0.7 0 0],'linewidth',3); hold on;
    end
    plot(ax6,ifil,misfit_zsta(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax6,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3);
%     ylabel(ax6,'$\mathbf{\delta Z\, (m)}$','fontsize',18,'Interpreter','latex')
    xlim(ax6,[0 Nfils+1]);
    ylim(ax6,[4 max(misfit_zsta)+10^(floor(log10(max(misfit_zsta))))*2]);
    yticks(ax6,[0.001 0.01 0.1 1 10 100 1000]);
    
    if ifil==1
        plot(ax7,[0 Nfils+1],[TAT_rms TAT_rms],'-','color',[0.7 0 0],'linewidth',3); hold on;
    end
    plot(ax7,ifil,misfit_TAT(ifil)*1000,symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax7,'yscale','log','linewidth',1.5,'fontsize',16,'xticklabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3);
%     ylabel(ax7,'{$\delta$\boldmath$\tau$ (\textbf{ms})}','fontsize',18,'Interpreter','latex')
    xlim(ax7,[0 Nfils+1]);
    ylim(ax7,[1 10]);
%     yticks(ax7,[0.001 0.01 0.1 1 10 100]);
    
    if ifil==1
        plot(ax8,[0 Nfils+1],[Vp_rms Vp_rms],'-','color',[0.7 0 0],'linewidth',3); hold on;
    end
    plot(ax8,ifil,misfit_Vw(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax8,'yscale','log','linewidth',1.5,'fontsize',16,'xticklabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3);
%     ylabel(ax8,'$\mathbf{\delta V_{P} \, (m/s)}$','fontsize',18,'Interpreter','latex')
    xlim(ax8,[0 Nfils+1]);
    ylim(ax8,[1 max(misfit_Vw)+10^(floor(log10(max(misfit_Vw))))]);
    yticks(ax8,[0.001 0.01 0.1 1 10 100]);
    
    %% Make shape legend
    dy = 0.958/(Nfils+0.25);
    ax(ifil) = axes('pos',[0.75 1-dy*ifil 0.05 0.05]);
    plot(ax(ifil),data_all(ifil).survx,data_all(ifil).survy,'-k','linewidth',1.5); hold on;
    if any(regexp(lgd{ifil},'line2'))
        plot(ax(ifil),0,-0.5,'ok','linewidth',2,'MarkerFaceColor',[0 0 0],'markersize',4);
    elseif any(regexp(lgd{ifil},'line 1'))
        plot(ax(ifil),0,0,'ok','linewidth',2,'MarkerFaceColor',[0 0 0],'markersize',4);
    end
    axis equal;
    xlim([-4 4]);
    ax(ifil).Visible = 'off';
end
l = legend(h,lgd,'position',[0.8 0.05 0.1852 0.95],'interpreter','none','fontsize',14,'box','off');

%%
% Plot all surveys
fig4 = figure(4); clf;
% set(gcf,'position',[53 292 1157 413]);
% clr = parula(Nfils);
for ifil = 1:Nfils
%     plot(data_all(ifil).x_ship,data_all(ifil).y_ship,symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',13,'linewidth',1); hold on;
    plot(data_all(ifil).x_ship,data_all(ifil).y_ship,'-','color',clr(ifil,:),'linewidth',3); hold on;
    grid off; box on; set(gca,'fontsize',16,'linewidth',2);
%     title([lgd{ifil}],'fontsize',18,'Interpreter','none');
    xlabel('X (m)','fontsize',18);
    ylabel('Y (m)','fontsize',18);
    axis equal;
end


%% SAVE
if ifsave
    save2pdf(ofile,f906)
end

end