% Compare synthetic survey model residuals
%

clear; close all;

%% INPUTS - MAKE SURE THESE ARE 
% path to project
projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA/'; 
% path to survey data from the project directory
datapath = '/Users/russell/Lamont/PROJ_OBSrange/synth_tests_paper/synth_surveys/'; 
% path to output directory from project directory(will be created if it does not yet exist)
outdir = [projpath,'/OUT_OBSrange_synthsurveys/OUT_wcorr_1pts/']; 

% prepend functions directory to MATLAB path
fullMAINpath = mfilename('fullpath');
functionspath = [fullMAINpath(1:regexp(fullMAINpath,mfilename)-1),'functions'];
addpath(functionspath);

%% Load *.mat files
files = dir([outdir,'/mats/*.mat']);

Nfils = length(files);
for ifil = 1:Nfils
    load([outdir,'/mats/',files(ifil).name]);
    
    lat_drop = data(1).drop(1);
    lon_drop = data(1).drop(2);
    z_drop = data(1).drop(3)*-1000;
    lats_ship = data(1).survlats;
    lons_ship = data(1).survlons;
    olon = lon_drop;
    olat = lat_drop;
    [ x_ship, y_ship ] = lonlat2xy_nomap( olon, olat, lons_ship, lats_ship );
    data_all(ifil).x_ship = x_ship;
    data_all(ifil).y_ship = y_ship;
    data_all(ifil).v_surv_true = data(1).v_surv_true;
    data_all(ifil).vmag = sqrt(sum(data(1).v_surv_true.^2,2));
    data_all(ifil).survx = data(1).survx;
    data_all(ifil).survy = data(1).survy;
    data_all(ifil).survey = data(1).survey;
    data_all(ifil).radius = data(1).radius;
    
    
    misfit_xsta(ifil,1) = rms(data(1).misfit_xsta);
    misfit_xsta_std(ifil,1) = std(data(1).misfit_xsta);
    misfit_ysta(ifil,1) = rms(data(1).misfit_ysta);
    misfit_ysta_std(ifil,1) = std(data(1).misfit_ysta);
    misfit_zsta(ifil,1) = rms(data(1).misfit_zsta);
    misfit_zsta_std(ifil,1) = std(data(1).misfit_zsta);
    misfit_r_xy(ifil,1) = rms(data(1).misfit_r_xy);
    misfit_r_xy_std(ifil,1) = std(data(1).misfit_r_xy);
    misfit_r_xyz(ifil,1) = rms(data(1).misfit_r_xyz);
    misfit_r_xyz_std(ifil,1) = std(data(1).misfit_r_xyz);
    misfit_TAT(ifil,1) = rms(data(1).misfit_TAT);
    misfit_TAT_std(ifil,1) = std(data(1).misfit_TAT);
    misfit_Vw(ifil,1) = rms(data(1).misfit_Vw);
    misfit_Vw_std(ifil,1) = std(data(1).misfit_Vw);
    E_rms(ifil,1) = mean(data(1).E_rms);
    E_rms_std(ifil,1) = std(data(1).E_rms);
    misfit_v_ship_all(1:2,ifil) = mean(data(1).misfit_v_ship_all,2);
    misfit_v_ship_all_std(1:2,ifil) = std(data(1).misfit_v_ship_all,0,2);
    misfit_dtwtcorr_all(ifil,1) = rms(data(1).misfit_dtwtcorr_all);
    misfit_dtwtcorr_all_std(ifil,1) = std(data(1).misfit_dtwtcorr_all);
    dtwt_all(ifil,1) = rms(data(1).dtwt_all);
    dtwt_all_std(ifil,1) = std(data(1).dtwt_all);
    
    lgd{ifil} = [data_all(ifil).survey,' ',num2str(data_all(ifil).radius),' nm']; %files(ifil).name(1:end-4);
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
    end
    
end
%% Plots
% Model Parameters
fig1 = figure(1); clf;
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
clr = [brewermap(7,'Blues'); brewermap(2,'Greens'); brewermap(1,'Reds'); brewermap(2,'Purples'); brewermap(2,'RdPu')];
for ifil = 1:Nfils
    h(ifil) = plot(ax1,ifil,misfit_r_xy(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax1,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[]);
    ylabel(ax1,'$\mathbf{\delta r_{xy}\, (m)}$','fontsize',18,'Interpreter','latex')
    xlim(ax1,[0 Nfils+1]);   

    plot(ax2,ifil,misfit_zsta(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax2,'yscale','log','linewidth',1.5,'fontsize',16,'XTickLabel',[]);
    ylabel(ax2,'$\mathbf{\delta Z\, (m)}$','fontsize',18,'Interpreter','latex')
    xlim(ax2,[0 Nfils+1]);
    
    plot(ax3,ifil,misfit_TAT(ifil)*1000,symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax3,'yscale','log','linewidth',1.5,'fontsize',16,'xticklabel',[]);
    ylabel(ax3,'$\mathbf{\delta TAT\, (ms)}$','fontsize',18,'Interpreter','latex')
    xlim(ax3,[0 Nfils+1]);
    ylim(ax3,[2.9 3.9]);

    plot(ax4,ifil,misfit_Vw(ifil),symbol{ifil},'markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    set(ax4,'yscale','log','linewidth',1.5,'fontsize',16);
    ylabel(ax4,'$\mathbf{\delta V_{H_2O} \, (m/s)}$','fontsize',18,'Interpreter','latex')
    xlim(ax4,[0 Nfils+1]);
    
    % Make shape legend
    dy = 0.95/(Nfils+0.25);
    ax(ifil) = axes('pos',[0.75 1-dy*ifil 0.05 0.05]); hold(ax1,'on');
    plot(ax(ifil),data_all(ifil).survx,data_all(ifil).survy,'-k','linewidth',1.5);
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


%% Save output
if ~exist(outdir)
    mkdir(outdir);
end

% Save plots
if ~exist([outdir,'/plots/'])
    mkdir([outdir,'/plots/']);
end
%     save2pdf([outdir,'/plots/',data.sta,'_1_OBSlocation.pdf'],f1,500)
%     save2pdf([outdir,'/plots/',data.sta,'_2_misfit.pdf'],f2,500)
%     save2pdf([outdir,'/plots/',data.sta,'_3_VelCorrs.pdf'],f3,500)
save2pdf([outdir,'/plots/','1_Comparisons.pdf'],fig1,500)
save2pdf([outdir,'/plots/','2_SurveyPatterns.pdf'],fig4,500)
