%% Function to produce Figure 06 & Supplemental figure 02
%  Figure 06 compare survey pattern geometries
%  Figure S02 fitting function to drxy-radius plot
%  This version calculates for different depths
function Figure06
addpath('/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/OBSrange_v1_MATLAB/functions/');
ofile = '../Figure06';
ofile_supp = '../FigureS2';
ifsave = 1;

% projname = 'mats_SynthBoot_summary_noTAT_REVISION1';
projname = 'mats_SynthBoot_summary_noTAT_REVISION1_fixtat';

depth = 2000; %5000, 2000, 500 % station depth (m)
% surveynames = {'PACMAN','cross','diamond','line','tri','circle'};
surveynames = {'PACMAN','cross','diamond','line','tri','circle','hourglass','cardinal'};

x_rms = 100;
y_rms = 100;
rxy_rms = sqrt(x_rms^2+y_rms^2);
z_rms = 50;
TAT_rms = 3; % ms
Vp_rms = 10;
FONTSIZE = 14;

markersize = 12;
%clr = parula(Nfils);
blues = brewermap(9,'Blues');
purp = brewermap(3,'Purples');
Nlines = 2;
purp2 = brewermap(4,'Purples');
grays = brewermap(4,'Greys');
greens = brewermap(4,'Greens');
pacmanclr = blues(end-6:end,:);
if length(surveynames) == 6
%           PACMAN        circle              cross2 cross        diamond    line2    line      tri_center  tri_edge
    clr = [pacmanclr; pacmanclr([3,5,7],:); [0 0 0; 0.5 0.5 0.5]; [0 0 0]; [0 0 0; 0.5 0.5 0.5]; [0 0 0; 0.5 0.5 0.5]];
%     Ilgd = [5 9 11:17];
    Ilgd = [5 9 13 12 11 15 14 16 17];
    dy_lgd = 0.5;
    dy_lgd_shift = 0;
elseif length(surveynames) == 8
%           PACMAN   cardinal      circle              cross2 cross        diamond hourglass  line2    line      tri_center  tri_edge
    clr = [pacmanclr; [0 0 0]; pacmanclr([3,5,7],:); [0 0 0; 0.5 0.5 0.5]; [0 0 0]; [0 0 0]; [0 0 0; 0.5 0.5 0.5]; [0 0 0; 0.5 0.5 0.5]];
%     Ilgd = [5 8 10 12:19];
    Ilgd = [5 10 13 12 17 16 18 19 14 15 8];
    dy_lgd = 0.6;
    dy_lgd_shift = -0.1;
end
%% Load *.mat files
files = dir(['../figdata/',projname,'/*.mat']);

Nfils = length(files);
ifil = 0;
for ifils = 1:Nfils
    if ~any(regexp(files(ifils).name,['z',num2str(depth),'m']))
        continue
    end
    if cellfun('isempty',regexp(files(ifils).name,surveynames))
        disp(files(ifils).name);
        continue
    end
    load(['../figdata/',projname,'/',files(ifils).name]);
    ifil = ifil+1;
    
    data_all(ifil).x_ship = data_summary.x_ship;
    data_all(ifil).y_ship = data_summary.y_ship;
    data_all(ifil).v_surv_true = data_summary.v_surv_true;
    data_all(ifil).vmag = data_summary.vmag;
    data_all(ifil).survx = data_summary.survx;
    data_all(ifil).survy = data_summary.survy;
    data_all(ifil).survt = data_summary.survt;
    data_all(ifil).survey = data_summary.survey;
    data_all(ifil).radius = data_summary.radius;
    data_all(ifil).survey_length = sum(sqrt(diff(data_all(ifil).survx).^2+diff(data_all(ifil).survy).^2))*1000;
    
    % If cardinal, load diamond shape
    if strcmp('cardinal',data_all(ifil).survey)
        diamond = load(['../figdata/',projname,'/SynthBoot_diamond_rad1.00_z5000m_fr10_OUT.mat']);
        data_all(ifil).survx = data_all(ifil).x_ship/1000;
        data_all(ifil).survy = data_all(ifil).y_ship/1000;
        data_all(ifil).survt = diamond.data_summary.survt;
        data_all(ifil).survey_length = sum(sqrt(diff(diamond.data_summary.survx).^2+diff(diamond.data_summary.survy).^2))*1000;        
    end
    
    misfit_xsta(ifil,1) = data_summary.misfit_xsta;
    misfit_xsta_std(ifil,1) = data_summary.misfit_xsta_std;
    misfit_xsta_mean(ifil,1) = data_summary.misfit_xsta_mean;
    misfit_ysta(ifil,1) = data_summary.misfit_ysta;
    misfit_ysta_std(ifil,1) = data_summary.misfit_ysta_std;
    misfit_ysta_mean(ifil,1) = data_summary.misfit_ysta_mean;
    misfit_zsta(ifil,1) = data_summary.misfit_zsta;
    misfit_zsta_std(ifil,1) = data_summary.misfit_zsta_std;
    misfit_zsta_mean(ifil,1) = data_summary.misfit_zsta_mean;
    misfit_zsta_absmean(ifil,1) = data_summary.misfit_zsta_absmean;
    misfit_zsta_med(ifil,1) = abs(data_summary.misfit_zsta_med);
    misfit_r_xy(ifil,1) = data_summary.misfit_r_xy;
    misfit_r_xy_std(ifil,1) = data_summary.misfit_r_xy_std;
%     misfit_r_xy_mean(ifil,1) = sqrt(data_summary.misfit_xsta_mean.^2 + data_summary.misfit_ysta_mean.^2); %abs(data_summary.misfit_r_xy_mean);
    misfit_r_xy_mean(ifil,1) = abs(data_summary.misfit_r_xy_mean);
    misfit_r_xy_med(ifil,1) = abs(data_summary.misfit_r_xy_med);
    misfit_r_xyz(ifil,1) = data_summary.misfit_r_xyz;
    misfit_r_xyz_std(ifil,1) = data_summary.misfit_r_xyz_std;
    misfit_TAT(ifil,1) = data_summary.misfit_TAT;
    misfit_TAT_std(ifil,1) = data_summary.misfit_TAT_std;
    misfit_TAT_mean(ifil,1) = data_summary.misfit_TAT_mean;
    misfit_TAT_absmean(ifil,1) = data_summary.misfit_TAT_absmean;
    misfit_TAT_med(ifil,1) = abs(data_summary.misfit_TAT_med);
    misfit_Vw(ifil,1) = data_summary.misfit_Vw;
    misfit_Vw_std(ifil,1) = data_summary.misfit_Vw_std;
    misfit_Vw_mean(ifil,1) = data_summary.misfit_Vw_mean;
    misfit_Vw_absmean(ifil,1) = data_summary.misfit_Vw_absmean;
    misfit_Vw_med(ifil,1) = abs(data_summary.misfit_Vw_med);
    E_rms(ifil,1) = data_summary.E_rms;
    E_rms_std(ifil,1) = data_summary.E_rms_std;
    misfit_v_ship_all(1:2,ifil) = data_summary.misfit_v_ship_all;
    misfit_v_ship_all_std(1:2,ifil) = data_summary.misfit_v_ship_all_std;
    misfit_dtwtcorr_all(ifil,1) = data_summary.misfit_dtwtcorr_all;
    misfit_dtwtcorr_all_std(ifil,1) = data_summary.misfit_dtwtcorr_all_std;
    dtwt_all(ifil,1) = data_summary.dtwt_all;
    dtwt_all_std(ifil,1) = data_summary.dtwt_all_std;
    
    if any(regexp(files(ifils).name,'PACMAN'))
        x = data_all(ifil).survx;
        y = data_all(ifil).survy;
        t = data_all(ifil).survt;
        r = sqrt((x-mean(x)).^2+(y-mean(y)).^2)*1000;
        r_mean(ifil) = mean(r);
        r_max(ifil) = max(r)/1852;
        dist_tot(ifil) = data_all(ifil).survey_length/1000;
        dt(ifil) = max(t)/60;
        misfit(ifil) = misfit_r_xy(ifil);
    %     lam(ifil) = r_max(ifil) * misfit_r_xy(ifil);
        lam(ifil) = dist_tot(ifil) * misfit_r_xy(ifil);
    else
        lam(ifil) = 9999;
    end
    
%     lgd{ifil} = [data_all(ifil).survey,' ',num2str(data_all(ifil).radius),' Nm']; %files(ifil).name(1:end-4);
    lgd{ifil} = [data_all(ifil).survey];
    if any(regexp(lgd{ifil},'PACMAN'))
        symbol{ifil} = 'ok';
    elseif any(regexp(lgd{ifil},'cross'))
        symbol{ifil} = '+k';
    elseif any(regexp(lgd{ifil},'diamond'))
        symbol{ifil} = 'dk';
    elseif any(regexp(lgd{ifil},'line'))
        symbol{ifil} = 'pk';
    elseif any(regexp(lgd{ifil},'tri'))
        symbol{ifil} = '^k';
    elseif any(regexp(lgd{ifil},'circle'))
        symbol{ifil} = 'ok';
    elseif any(regexp(lgd{ifil},'hourglass'))
        symbol{ifil} = 'sk';
    elseif any(regexp(lgd{ifil},'cardinal'))
        symbol{ifil} = '*k';
    else
        continue
    end
    
    fprintf('%s %.2f Nm : zmisfit %.2f\n',lgd{ifil},data_all(ifil).radius,misfit_zsta(ifil,1));
    
end

%% Show mean errors
disp([misfit_xsta_mean, misfit_ysta_mean, misfit_zsta_mean, misfit_TAT_mean, misfit_Vw_mean]);

cmap = cmocean('balance');
figure(1); set(gcf,'position',[132     1   840   697]);
axb = subplot(1,2,1);
% Plot mean errors in percent of RMS error
imagesc(axb,100*[misfit_xsta_mean/100, misfit_ysta_mean/100, misfit_zsta_mean./50, misfit_TAT_mean./0.003, misfit_Vw_mean./10])
colormap(axb,cmap);
colorbar;
caxis(axb,[-100 100]);
yticks(axb,[1:length(files)]);
xticks(axb,[1:5]);
title('Mean error (% of RMS)','fontsize',18);

axa = subplot(1,2,2);
% Plot mean errors
imagesc(axa,[misfit_xsta_mean, misfit_ysta_mean, misfit_zsta_mean, misfit_TAT_mean*1000, misfit_Vw_mean])
colormap(axa,parula);
colorbar;
yticks(axa,[1:length(files)]);
xticks(axa,[1:5]);
title('Mean error (true units)','fontsize',18);

%% Plots
% Model Parameters
f906 = figure(906); clf;
set(gcf,'position',[253     6   666   697],'color','w')
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

iline = 0;
iilgd = 0;
Nfils = length(data_all);
[~,imin] = min(lam);
for ifil = 1:Nfils
    % MEAN
    if ifil==1
        plot(ax1,[0 Nfils+1],[rxy_rms rxy_rms],'-','color',[0.7 0 0],'linewidth',3); hold on;
    end
    h(ifil) = plot(ax1,ifil,misfit_r_xy_mean(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    if any(regexp(lgd{ifil},'line'))
        iline = iline + 1;
        plot(ax1,ifil,misfit_xsta_mean(ifil),symbol{ifil},'markeredgecolor',purp2(end-(Nlines-iline),:),'markersize',markersize,'linewidth',2); hold on;
    end
    set(ax1,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3,'layer','top');
    ylabel(ax1,'$\mathbf{\delta r_{xy}\, (m)}$','fontsize',18,'Interpreter','latex')
    xlim(ax1,[0 Nfils+1]);   
    title(ax1,'\textbf{MAE}','interpreter','latex','fontsize',18);
    ylim(ax1,[0.7 max(misfit_r_xy)+10^(floor(log10(max(misfit_r_xy))))*4]);
    yticks(ax1,[0.001 0.01 0.1 1 10 100 1000]);
    
    if ifil==1
        plot(ax2,[0 Nfils+1],[z_rms z_rms],'-','color',[0.7 0 0],'linewidth',3); hold on;
    end
    plot(ax2,ifil,misfit_zsta_absmean(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax2,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3,'layer','top');
    ylabel(ax2,'$\mathbf{\delta Z\, (m)}$','fontsize',18,'Interpreter','latex')
    xlim(ax2,[0 Nfils+1]);
    ylim(ax2,[3 max(misfit_zsta_absmean)+10^(floor(log10(max(misfit_zsta_absmean))))*5]);
    yticks(ax2,[0.1 1 10 100]);
    
    if ifil==1
        semilogy(ax3,[0 Nfils+1],[TAT_rms TAT_rms],'-','color',[0.7 0 0],'linewidth',3); hold on;
    end
    plot(ax3,ifil,misfit_TAT_absmean(ifil)*1000,symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax3,'yscale','log','linewidth',1.5,'fontsize',16,'xticklabel',[],'yminortick','on','xtick',[],'TickLength',[0.01, 0.001]*3,'layer','top');
    ylabel(ax3,'{$\delta$\boldmath$\tau$ (\textbf{ms})}','fontsize',18,'Interpreter','latex')
    xlim(ax3,[0 Nfils+1]);
%     yticks(ax3,[0.001 0.01 0.1 1 10 100 1000]);
    ylim(ax3,[1 10]);
%     yticks(ax3,[0:0.1:0.4]);
    
    if ifil==1
        plot(ax4,[0 Nfils+1],[Vp_rms Vp_rms],'-','color',[0.7 0 0],'linewidth',3); hold on;
    end
    plot(ax4,ifil,misfit_Vw_absmean(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
%     set(ax4,'yscale','log','linewidth',1.5,'fontsize',16,'xticklabel',[],'xtick',[],'TickLength',[0.01, 0.001]*3,'layer','top');
%     ylabel(ax4,'$\mathbf{\delta V_{P} \, (m/s)}$','fontsize',18,'Interpreter','latex')
    xlim(ax4,[0 Nfils+1]);
    ylim(ax4,[0.8 max(misfit_Vw_absmean)+10^(floor(log10(max(misfit_Vw_absmean))))*2]);
%     yticks(ax4,[0.001 0.01 0.1 1 10 100]);
end
%% Place axes
delete(ax1); delete(ax2); delete(ax3); %delete(ax4); %delete(ax8);
dx_shift = -0.33;
dx = 1.75;
dy_shift = 0.012; %0.0;%-0.08;
dy_static = -0.032;
dy = 1;%1.4;
ax5.Position = [ax5.Position(1)+dx_shift, ax5.Position(2)+dy_shift*0+dy_static,   ax5.Position(3)*dx, ax5.Position(4)*dy];
ax6.Position = [ax6.Position(1)+dx_shift, ax6.Position(2)+dy_shift*1+dy_static, ax6.Position(3)*dx, ax6.Position(4)*dy];
ax7.Position = [ax7.Position(1)+dx_shift, ax7.Position(2)+dy_shift*2+dy_static, ax7.Position(3)*dx, ax7.Position(4)*dy];

dy_shift = -0.005;%-0.08;
dx = 0.7;
dx_shift = -0.33;
ax8.Position = [ax8.Position(1)+dx_shift, ax8.Position(2)+dy_shift, ax8.Position(3)*dx, ax8.Position(4)*dy];

ax4.Position = ax8.Position;
dx_shift = 0.24;
ax4.Position = [ax4.Position(1)+dx_shift, ax4.Position(2), ax4.Position(3), ax4.Position(4)];
for ifil = 1:Nfils
    %% RMS
    lw = 2;
    if ifil==1
%         plot(ax5,[0 Nfils+1],[rxy_rms rxy_rms],'-','color',[0.7 0 0],'linewidth',3); hold on;
    end
    if any(regexp(lgd{ifil},'line'))
        h(ifil) = plot(ax5,data_all(ifil).survey_length/1000+[-0.5 +0.5],misfit_r_xy(ifil)*[1 1],'-','color',clr(ifil,:),'linewidth',lw); hold on;
    elseif ~any(regexp(lgd{ifil},'PACMAN')) && ~any(regexp(lgd{ifil},'line'))
        h(ifil) = plot(ax5,data_all(ifil).survey_length/1000,misfit_r_xy(ifil),symbol{ifil},'color',clr(ifil,:),'markersize',markersize,'linewidth',lw); hold on;
    else
        h(ifil) = plot(ax5,data_all(ifil).survey_length/1000,misfit_r_xy(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    end
    if ifil == imin
        plot(ax5,data_all(ifil).survey_length/1000,misfit_r_xy(ifil),symbol{ifil},'color','r','markersize',markersize,'linewidth',2); hold on;
    end
    if any(regexp(lgd{ifil},'line'))
%         plot(ax5,data_all(ifil).survey_length/1000,misfit_xsta(ifil),symbol{ifil},'markeredgecolor',purp2(end-(Nlines-iline),:),'markersize',markersize,'linewidth',2); hold on;
        plot(ax5,data_all(ifil).survey_length/1000+[-0.5 +0.5],misfit_xsta(ifil)*[1 1],'-.','color',clr(ifil,:),'linewidth',2); hold on;
    end
    set(ax5,'yscale','log','linewidth',1.5,'fontsize',FONTSIZE,'XAxisLocation','top','TickLength',[0.01, 0.001]*2,'layer','top','xminortick','on','TickDir','out');
%     set(ax5,'yscale','linear','linewidth',1.5,'fontsize',FONTSIZE,'XAxisLocation','top','TickLength',[0.01, 0.001]*3,'layer','top','xminortick','on');    
%     ylim(ax5,ax1.YLim);
%     ylim(ax5,[0.7 max(misfit_r_xy)+10^(floor(log10(max(misfit_r_xy))))*4]);
    ylim(ax5,[min(misfit_r_xy)-0.3 100]);
    xlim(ax5,[0 20]);
    title(ax5,'\textbf{RMS}','interpreter','latex','fontsize',FONTSIZE);
    yticks(ax5,[0.001 0.01 0.1 1 10 100 1000]);
    ylabel(ax5,'$\mathbf{\delta r_{xy}\, (m)}$','fontsize',FONTSIZE,'Interpreter','latex')
    xlabel(ax5,'\bf{Total survey distance (km)}','fontsize',FONTSIZE,'Interpreter','latex');
    if misfit_r_xy(ifil) > ax5.YLim(2)                 % Difference
        qax5=axes('position',get(ax5,'position'));  % overlay second axis on top first
        Xdata = data_all(ifil).survey_length/1000;
        Ydata = ax5.YLim(2);
        Ypl = Ydata-15;
        plot_arrows( ax5, Xdata, Ydata )
        if any(regexp(lgd{ifil},'line'))
            plot(qax5,data_all(ifil).survey_length/1000+[-0.5 +0.5],Ypl*[1 1],'-','color',clr(ifil,:),'linewidth',lw); hold on;
        elseif ~any(regexp(lgd{ifil},'PACMAN')) && ~any(regexp(lgd{ifil},'line'))
            plot(qax5,data_all(ifil).survey_length/1000,Ypl,symbol{ifil},'color',clr(ifil,:),'markersize',markersize,'linewidth',lw); hold on;
        else
            plot(qax5,data_all(ifil).survey_length/1000,Ypl,symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
        end
        set(qax5,'Visible','off','XLim',ax5.XLim,'YLim',ax5.YLim);
    end
    
    if ifil==1
%         plot(ax6,[0 Nfils+1],[z_rms z_rms],'-','color',[0.7 0 0],'linewidth',3); hold on;
    end
    if any(regexp(lgd{ifil},'line'))
        plot(ax6,data_all(ifil).survey_length/1000+[-0.5 +0.5],misfit_zsta(ifil)*[1 1],'-','color',clr(ifil,:),'linewidth',lw); hold on;
    elseif ~any(regexp(lgd{ifil},'PACMAN')) && ~any(regexp(lgd{ifil},'line'))
        plot(ax6,data_all(ifil).survey_length/1000,misfit_zsta(ifil),symbol{ifil},'color',clr(ifil,:),'markersize',markersize,'linewidth',lw); hold on;
    else
        plot(ax6,data_all(ifil).survey_length/1000,misfit_zsta(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    end
    if ifil == imin
        plot(ax6,data_all(ifil).survey_length/1000,misfit_zsta(ifil),symbol{ifil},'color','r','markersize',markersize,'linewidth',2); hold on;
    end
    set(ax6,'yscale','log','linewidth',1.5,'fontsize',FONTSIZE,'XTickLabel',[],'TickLength',[0.01, 0.001]*2,'layer','top','xminortick','on','TickDir','out');
%     set(ax6,'yscale','linear','linewidth',1.5,'fontsize',FONTSIZE,'XTickLabel',[],'TickLength',[0.01, 0.001]*3,'layer','top','xminortick','on');
%     ylim(ax6,ax2.YLim);
%     ylim(ax6,[min(misfit_zsta)-0.2 max(misfit_zsta)+10^(floor(log10(max(misfit_zsta))))*5]);
    ylim(ax6,[min(misfit_zsta)-0.4 100]);
    xlim(ax6,[0 20]);
    yticks(ax6,[0.001 0.01 0.1 1 10 100 1000]);
    ylabel(ax6,'$\mathbf{\delta Z\, (m)}$','fontsize',FONTSIZE,'Interpreter','latex')
    if misfit_zsta(ifil) > ax6.YLim(2)                      % Difference
        qax6=axes('position',get(ax6,'position'));  % overlay second axis on top first
        Xdata = data_all(ifil).survey_length/1000;
        Ydata = ax6.YLim(2);
        Ypl = Ydata-15;
        plot_arrows( ax6, Xdata, Ydata )
        if any(regexp(lgd{ifil},'line'))
            plot(qax6,data_all(ifil).survey_length/1000+[-0.5 +0.5],Ypl*[1 1],'-','color',clr(ifil,:),'linewidth',lw); hold on;
        elseif ~any(regexp(lgd{ifil},'PACMAN')) && ~any(regexp(lgd{ifil},'line'))
            plot(qax6,data_all(ifil).survey_length/1000,Ypl,symbol{ifil},'color',clr(ifil,:),'markersize',markersize,'linewidth',lw); hold on;
        else
            plot(qax6,data_all(ifil).survey_length/1000,Ypl,symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
        end
        set(qax6,'Visible','off','XLim',ax6.XLim,'YLim',ax6.YLim);
    end
    
    if ifil==1
%         plot(ax7,[0 Nfils+1],[Vp_rms Vp_rms],'-','color',[0.7 0 0],'linewidth',3); hold on;
    end
    if any(regexp(lgd{ifil},'line'))
        plot(ax7,data_all(ifil).survey_length/1000+[-0.5 +0.5],misfit_Vw(ifil)*[1 1],'-','color',clr(ifil,:),'linewidth',lw); hold on;
    elseif ~any(regexp(lgd{ifil},'PACMAN')) && ~any(regexp(lgd{ifil},'line'))
        plot(ax7,data_all(ifil).survey_length/1000,misfit_Vw(ifil),symbol{ifil},'color',clr(ifil,:),'markersize',markersize,'linewidth',lw); hold on;
    else
        plot(ax7,data_all(ifil).survey_length/1000,misfit_Vw(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    end
    if ifil == imin
        plot(ax7,data_all(ifil).survey_length/1000,misfit_Vw(ifil),symbol{ifil},'color','r','markersize',markersize,'linewidth',2); hold on;
    end
    set(ax7,'yscale','log','linewidth',1.5,'fontsize',FONTSIZE,'XTickLabel',[],'TickLength',[0.01, 0.001]*2,'layer','top','xminortick','on','TickDir','out');
%     set(ax7,'yscale','linear','linewidth',1.5,'fontsize',FONTSIZE,'XTickLabel',[],'TickLength',[0.01, 0.001]*3,'layer','top','xminortick','on');
%     ylim(ax7,ax4.YLim);
%     ylim(ax7,[0.8 max(misfit_Vw)+10^(floor(log10(max(misfit_Vw))))*2]);
    ylim(ax7,[min(misfit_Vw)-0.2 100]);
    xlim(ax7,[0 20]);
    yticks(ax7,[0.001 0.01 0.1 1 10 100]);
    ylabel(ax7,'$\mathbf{\delta V_{P} \, (m/s)}$','fontsize',FONTSIZE,'Interpreter','latex')
    if misfit_Vw(ifil) > ax7.YLim(2)
        qax7=axes('position',get(ax7,'position'));  % overlay second axis on top first
        Xdata = data_all(ifil).survey_length/1000;
        Ydata = ax7.YLim(2);
        Ypl = Ydata-15;
        plot_arrows( ax7, Xdata, Ydata )        
        if any(regexp(lgd{ifil},'line'))
            plot(qax7,data_all(ifil).survey_length/1000+[-0.5 +0.5],Ypl*[1 1],'-','color',clr(ifil,:),'linewidth',lw); hold on;
        elseif ~any(regexp(lgd{ifil},'PACMAN')) && ~any(regexp(lgd{ifil},'line'))
            plot(qax7,data_all(ifil).survey_length/1000,Ypl,symbol{ifil},'color',clr(ifil,:),'markersize',markersize,'linewidth',lw); hold on;
        else
            plot(qax7,data_all(ifil).survey_length/1000,Ypl,symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
        end
        set(qax7,'Visible','off','XLim',ax7.XLim,'YLim',ax7.YLim);
    end
    %% Calculate diminishing return parameter
    if any(regexp(lgd{ifil},'PACMAN'))
%     x = data_all(ifil).survx;
%     y = data_all(ifil).survy;
%     t = data_all(ifil).survt;
%     r = sqrt((x-mean(x)).^2+(y-mean(y)).^2)*1000;
%     r_mean(ifil) = mean(r);
%     r_max(ifil) = max(r)/1852;
%     dist_tot(ifil) = data_all(ifil).survey_length/1000;
%     dt(ifil) = max(t)/60;
%     misfit(ifil) = misfit_r_xy(ifil);
% %     lam(ifil) = r_max(ifil) * misfit_r_xy(ifil);
%     lam(ifil) = dist_tot(ifil) * misfit_r_xy(ifil);
    
    plot(ax8,data_all(ifil).radius,lam(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    if ifil == imin
        plot(ax8,data_all(ifil).radius,lam(ifil),symbol{ifil},'color','r','markersize',markersize,'linewidth',2); hold on;
    end
    if any(regexp(lgd{ifil},'line'))
        lamx = r_max(ifil) * misfit_xsta(ifil);
        plot(ax8,ifil,lamx,symbol{ifil},'markeredgecolor',purp2(end-(Nlines-iline),:),'markersize',markersize,'linewidth',2); hold on;
    end
    set(ax8,'yscale','linear','linewidth',1.5,'fontsize',FONTSIZE,'xminortick','on','yminortick','on','TickLength',[0.01, 0.001]*3,'layer','top','TickDir','out');
%     ylabel(ax8,'{$\delta$\boldmath$\lambda$ (\textbf{m^2})}','fontsize',18,'Interpreter','latex')
    xlim(ax8,[min(r_max)-0.1 max(r_max)+0.1]);
%     ylim(ax8,[35 70]);
%     yticks(ax8,[0.001 0.01 0.1 1 10 100 1e3 1e4 1e5 1e6 1e7]);
    ylabel(ax8,'{\boldmath$\lambda$ (\textbf{m}$\cdot$\textbf{km})}','fontsize',FONTSIZE,'Interpreter','latex')
    xlabel(ax8,'\bf{Radius (Nm)}','fontsize',FONTSIZE,'Interpreter','latex')
    end
end

%% Make shape legend
% if ~isempty(intersect(ifil,Ilgd))
iilgd = 0;
for ifil = Ilgd
    iilgd = iilgd+1;
    Nfils = length(Ilgd);
    dy = dy_lgd/(Nfils+0.25);
    ax(ifil) = axes('pos',[0.85 0.93-dy*iilgd 0.08 0.08]);
    if strcmp('cardinal',data_all(ifil).survey)
        plot(ax(ifil),data_all(ifil).survx,data_all(ifil).survy,'ok','linewidth',1.5); hold on;
    %     fill(ax(ifil),data_all(ifil).survx,data_all(ifil).survy,[0.5 0.5 0.5]); hold on;
    else
        plot(ax(ifil),data_all(ifil).survx,data_all(ifil).survy,'-k','linewidth',1.5); hold on;
    end
    if any(regexp(lgd{ifil},'line2'))
        plot(ax(ifil),0,-0.5,'ok','linewidth',2,'MarkerFaceColor',[0 0 0],'markersize',3);
    elseif ~any(regexp(lgd{ifil},'line2')) %any(regexp(lgd{ifil},'line 1'))
        plot(ax(ifil),0,0,'ok','linewidth',2,'MarkerFaceColor',[0 0 0],'markersize',3);
    end
    axis equal;
    xlim([-4 4]);
    ax(ifil).Visible = 'off';
    if ifil == 5            
        text(ax(ifil),-0.15,1.05,...
            '\textbf{Geometry}','color',[0 0 0],'interpreter','latex','fontsize',14,'Units','normalized');
    end
end
% l = legend(h,lgd,'position',[0.75 0.05 0.1852 0.95],'interpreter','none','fontsize',14,'box','off');
l = legend(h(Ilgd),{lgd{Ilgd}},'position',[0.65 0.45+dy_lgd_shift 0.1852 dy_lgd],'interpreter','none','fontsize',14,'box','off');

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

%% fit line
% g = fittype('b*exp(-c*x)');
x = r_max';
% x = dist_tot';
y = misfit';

% g = fittype('a*(x)^(b)');
% f0 = fit(x,y,g,'StartPoint',[1,1]);

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,-Inf,-100],...
               'Upper',[Inf,Inf,100],...
               'StartPoint',[1 1,1]);
g = fittype('a*((x)^(b)+c)','options',fo);
f0 = fit(x,y,g);

xx = linspace(min(x),max(x),50);
xx_rad = linspace(min(r_max),max(r_max),50);
f12 = figure(12); clf;
for ifil = 1:7
    h(ifil) = plot(x(ifil),y(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
end
% plot(x,y,'o','linewidth',2,'markersize',10); hold on;
plot(xx,f0(xx),'-k','linewidth',2);
xlabel('Radius (Nm)','fontsize',14);
ylabel('$\mathbf{\delta r_{xy}\, (m)}$','fontsize',14,'Interpreter','latex')
xlabel('\bf{Radius (Nm)}','fontsize',14,'Interpreter','latex')
set(gca,'linewidth',1.5,'fontsize',14,'xminortick','on','yminortick','on');
% ylim([0 65]);

dx = diff(xx);
dx = dx(1);
df_dx = gradient(f0(xx),dx);
txt = sprintf('$y = %.2f(x^{%.2f} + %.2f)$',f0.a,f0.b,f0.c);
text(0.8,max(f0(xx))-1,txt,'interpreter','latex','fontsize',18);
% subplot(2,1,2);
% plot(xx,df_dx,'-b');

cla(ax4);
plot(ax4,xx,df_dx,'-','color',blues(7,:),'linewidth',2); hold on;
plot(ax4,[data_all(imin).radius data_all(imin).radius],[-90 0],'--r','linewidth',2);
set(ax4,'yscale','linear','linewidth',1.5,'fontsize',FONTSIZE,'TickLength',[0.01, 0.001]*3,...
    'layer','top','xminortick','on','yminortick','on','yaxislocation','right','TickDir','out');
% ylabel(ax4,'$\mathbf{m/Nm}$','fontsize',16,'Interpreter','latex')
ylabel(ax4,'$\mathbf{\nabla r_{xy}\, (m/Nm)}$','fontsize',FONTSIZE,'Interpreter','latex')
xlabel(ax4,'\bf{Radius (Nm)}','fontsize',FONTSIZE,'Interpreter','latex')
xlim(ax4,[min(r_max)-0.1 max(r_max)+0.1]);
ylim(ax4,[-90 0]);
% ylim(ax4,[0.8 max(misfit_Vw_absmean)+10^(floor(log10(max(misfit_Vw_absmean))))*2]);
% yticks(ax4,[0.001 0.01 0.1 1 10 100]);

%% Figure numbers
x = 0.92;
y = 0.9;
text(ax5,x,y,...
'\textbf{a)}','color',[0 0 0],'interpreter','latex','fontsize',18,'Units','normalized');
text(ax6,x,y,...
'\textbf{b)}','color',[0 0 0],'interpreter','latex','fontsize',18,'Units','normalized');
text(ax7,x,y,...
'\textbf{c)}','color',[0 0 0],'interpreter','latex','fontsize',18,'Units','normalized');

if depth == 500
    text(ax8,0.07,0.9,...
    '\textbf{d)}','color',[0 0 0],'interpreter','latex','fontsize',18,'Units','normalized');

    text(ax4,0.8,0.9,...
    '\textbf{e)}','color',[0 0 0],'interpreter','latex','fontsize',18,'Units','normalized');
else
    text(ax8,0.8,0.9,...
    '\textbf{d)}','color',[0 0 0],'interpreter','latex','fontsize',18,'Units','normalized');

    text(ax4,0.07,0.9,...
    '\textbf{e)}','color',[0 0 0],'interpreter','latex','fontsize',18,'Units','normalized');
end

%% SAVE
if ifsave
    save2pdf([ofile,'_z',num2str(depth),'m'],f906)
    save2pdf([ofile_supp,'_z',num2str(depth),'m'],f12);
end

end