%% Function to produce Figure 06 
%  Figure 06 compare survey pattern geometries
function Figure06

ofile = '../Figure06';
ifsave = 1;

x_rms = 100;
y_rms = 100;
rxy_rms = sqrt(x_rms^2+y_rms^2);
z_rms = 50;
TAT_rms = 3; % ms
Vp_rms = 10;

%% Load *.mat files
files = dir('../figdata/mats_SynthBoot_summary/*.mat');

Nfils = length(files);
for ifil = 1:Nfils
    load(['../figdata/mats_SynthBoot_summary/',files(ifil).name]);
    
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
    misfit_r_xy(ifil,1) = data_summary.misfit_r_xy;
    misfit_r_xy_std(ifil,1) = data_summary.misfit_r_xy_std;
    misfit_r_xy_mean(ifil,1) = abs(data_summary.misfit_r_xy_mean);
    misfit_r_xyz(ifil,1) = data_summary.misfit_r_xyz;
    misfit_r_xyz_std(ifil,1) = data_summary.misfit_r_xyz_std;
    misfit_TAT(ifil,1) = data_summary.misfit_TAT;
    misfit_TAT_std(ifil,1) = data_summary.misfit_TAT_std;
    misfit_TAT_mean(ifil,1) = abs(data_summary.misfit_TAT_mean);
    misfit_Vw(ifil,1) = data_summary.misfit_Vw;
    misfit_Vw_std(ifil,1) = data_summary.misfit_Vw_std;
    misfit_Vw_mean(ifil,1) = abs(data_summary.misfit_Vw_mean);
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
dy = 1.2;
ax1.Position = [ax1.Position(1), ax1.Position(2), ax1.Position(3), ax1.Position(4)*dy];
ax2.Position = [ax2.Position(1), ax2.Position(2), ax2.Position(3), ax2.Position(4)*dy];
ax3.Position = [ax3.Position(1), ax3.Position(2), ax3.Position(3), ax3.Position(4)*dy];
ax4.Position = [ax4.Position(1), ax4.Position(2), ax4.Position(3), ax4.Position(4)*dy];

markersize = 14;
%clr = parula(Nfils);

purp = brewermap(3,'Purples');
Nlines = 2;
purp2 = brewermap(4,'Purples');
clr = [brewermap(7,'Blues'); brewermap(3,'Greys'); brewermap(2,'Greens'); brewermap(1,'Reds'); purp(end-1:end,:); brewermap(2,'RdPu')];
iline = 0;
for ifil = 1:Nfils
    if ifil==1
        plot(ax1,[0 Nfils+1],[rxy_rms rxy_rms],'-','color',[0.7 0 0],'linewidth',1.5); hold on;
    end
    h(ifil) = plot(ax1,ifil,misfit_r_xy(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    if any(regexp(lgd{ifil},'line'))
        iline = iline+1
        plot(ax1,ifil,misfit_xsta(ifil),symbol{ifil},'markeredgecolor',purp2(end-(Nlines-iline),:),'markersize',markersize,'linewidth',2); hold on;
    end
    set(ax1,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3);
%     ylabel(ax5,'$\mathbf{\delta r_{xy}\, (m)}$','fontsize',18,'Interpreter','latex')
    xlim(ax1,[0 Nfils+1]);   
    ylim(ax1,[1 max(misfit_r_xy)+10^(floor(log10(max(misfit_r_xy))))*10]);
    title(ax1,'\textbf{RMS}','interpreter','latex','fontsize',18);
    yticks(ax1,[0.001 0.01 0.1 1 10 100 1000]);
    
    if ifil==1
        plot(ax2,[0 Nfils+1],[z_rms z_rms],'-','color',[0.7 0 0],'linewidth',1.5); hold on;
    end
    plot(ax2,ifil,misfit_zsta(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax2,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3);
%     ylabel(ax6,'$\mathbf{\delta Z\, (m)}$','fontsize',18,'Interpreter','latex')
    xlim(ax2,[0 Nfils+1]);
    ylim(ax2,[4 max(misfit_zsta)+10^(floor(log10(max(misfit_zsta))))*2]);
    yticks(ax2,[0.001 0.01 0.1 1 10 100 1000]);
    
    if ifil==1
        plot(ax3,[0 Nfils+1],[TAT_rms TAT_rms],'-','color',[0.7 0 0],'linewidth',1.5); hold on;
    end
    plot(ax3,ifil,misfit_TAT(ifil)*1000,symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax3,'yscale','log','linewidth',1.5,'fontsize',16,'xticklabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3);
%     ylabel(ax7,'{$\delta$\boldmath$\tau$ (\textbf{ms})}','fontsize',18,'Interpreter','latex')
    xlim(ax3,[0 Nfils+1]);
    ylim(ax3,[2.9 3.25]);
%     yticks(ax7,[0.001 0.01 0.1 1 10 100]);
    
    if ifil==1
        plot(ax4,[0 Nfils+1],[Vp_rms Vp_rms],'-','color',[0.7 0 0],'linewidth',1.5); hold on;
    end
    plot(ax4,ifil,misfit_Vw(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax4,'yscale','log','linewidth',1.5,'fontsize',16,'xticklabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3);
%     ylabel(ax8,'$\mathbf{\delta V_{P} \, (m/s)}$','fontsize',18,'Interpreter','latex')
    xlim(ax4,[0 Nfils+1]);
    ylim(ax4,[1 max(misfit_Vw)+10^(floor(log10(max(misfit_Vw))))]);
    yticks(ax4,[0.001 0.01 0.1 1 10 100]);
    
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