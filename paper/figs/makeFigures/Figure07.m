%% Function to produce Figure 07
%  Figure 07 Resolution and covariance matrices
function Figure07

ofile = '../Figure07';
ifsave = 0;

files = {
        'SynthBoot_PACMAN_rad1.00_OUT.mat';
        'SynthBoot_line_rad1.00_OUT.mat';
        'SynthBoot_circle_rad1.00_OUT.mat';
         };

%% Load *.mat files
% files = dir('../figdata/mats_SynthBoot_summary/*.mat');

Nfils = length(files);
for ifil = 1:Nfils
    load(['../figdata/mats_SynthBoot_summary/',files{ifil}]);
    
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
    misfit_ysta(ifil,1) = data_summary.misfit_ysta;
    misfit_ysta_std(ifil,1) = data_summary.misfit_ysta_std;
    misfit_zsta(ifil,1) = data_summary.misfit_zsta;
    misfit_zsta_std(ifil,1) = data_summary.misfit_zsta_std;
    misfit_r_xy(ifil,1) = data_summary.misfit_r_xy;
    misfit_r_xy_std(ifil,1) = data_summary.misfit_r_xy_std;
    misfit_r_xyz(ifil,1) = data_summary.misfit_r_xyz;
    misfit_r_xyz_std(ifil,1) = data_summary.misfit_r_xyz_std;
    misfit_TAT(ifil,1) = data_summary.misfit_TAT;
    misfit_TAT_std(ifil,1) = data_summary.misfit_TAT_std;
    misfit_Vw(ifil,1) = data_summary.misfit_Vw;
    misfit_Vw_std(ifil,1) = data_summary.misfit_Vw_std;
    E_rms(ifil,1) = data_summary.E_rms;
    E_rms_std(ifil,1) = data_summary.E_rms_std;
    misfit_v_ship_all(1:2,ifil) = data_summary.misfit_v_ship_all;
    misfit_v_ship_all_std(1:2,ifil) = data_summary.misfit_v_ship_all_std;
    misfit_dtwtcorr_all(ifil,1) = data_summary.misfit_dtwtcorr_all;
    misfit_dtwtcorr_all_std(ifil,1) = data_summary.misfit_dtwtcorr_all_std;
    dtwt_all(ifil,1) = data_summary.dtwt_all;
    dtwt_all_std(ifil,1) = data_summary.dtwt_all_std;
    R_mat{ifil} = data_summary.R_mat;
    Cm_mat{ifil} = data_summary.Cm_mat;
    
    [ExpSigma,Cm_mat{ifil}] = cov2corr(Cm_mat{ifil});
    
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
% % Model Parameters
% f906 = figure(906); clf;
% set(gcf,'position',[1     1   918   697])
% ax1 = subplot(4,4,[1 3]); hold(ax1,'on'); box on;
% ax2 = subplot(4,4,[5 7]); hold(ax2,'on'); box on;
% ax3 = subplot(4,4,[9 11]); hold(ax3,'on'); box on;
% ax4 = subplot(4,4,[13 15]); hold(ax4,'on'); box on;
% dy = 1.2;
% ax1.Position = [ax1.Position(1), ax1.Position(2), ax1.Position(3), ax1.Position(4)*dy];
% ax2.Position = [ax2.Position(1), ax2.Position(2), ax2.Position(3), ax2.Position(4)*dy];
% ax3.Position = [ax3.Position(1), ax3.Position(2), ax3.Position(3), ax3.Position(4)*dy];
% ax4.Position = [ax4.Position(1), ax4.Position(2), ax4.Position(3), ax4.Position(4)*dy];

%% Model resolution and Covariance
    f907 = figure(907); clf;
    set(gcf,'position',[91     1   575   697]);
    cmap = cmocean('balance');
    
for ifil = 1:length(files)    
    ax(ifil*2-1) = subplot(length(files),2,ifil*2-1);
    ax1Pos = ax(ifil*2-1).Position;
    imagesc(ax(ifil*2-1),R_mat{ifil}(:,:,1)); hold on;
    for i = 1:5
        plot([.5,5.5],[i-.5,i-.5],'k-','linewidth',2);
        plot([i-.5,i-.5],[.5,5.5],'k-','linewidth',2);
    end
    R_spread = sum(sum((R_mat{ifil}(:,:,1)-eye(size(R_mat{ifil}(:,:,1)))).^2));
    axis square;
    axis tight;
    set(ax(ifil*2-1),'fontsize',16,'linewidth',2, ...
        'XTickLabel',{'$\mathbf{X}$','$\mathbf{Y}$','$\mathbf{Z}$','{\boldmath$\tau$}','$\mathbf{V_{P}}$'},'TickLabelInterpreter','latex', ...
        'YTickLabel',{'$\mathbf{X}$','$\mathbf{Y}$','$\mathbf{Z}$','{\boldmath$\tau$}','$\mathbf{V_{P}}$'});
    title(ax(ifil*2-1),'\textbf{Model Resolution}','fontsize',18,'Interpreter','latex');
    cb1 = colorbar(ax(ifil*2-1));
    colormap(cmap)
%     caxis([-max(max(abs(mean(RR,3)))) max(max(abs(mean(RR,3))))]);
    caxis([-1 1]);
    
    ax(ifil*2) = subplot(length(files),2,ifil*2);
    ax2Pos = ax(ifil*2).Position;
%     imagesc(ax2,sign(Cm_mat{ifil}(:,:,1)).*log10(abs(Cm_mat{ifil}(:,:,1)))); hold on;
    imagesc(ax(ifil*2),Cm_mat{ifil}(:,:,1)); hold on;
    for i = 1:5
        plot([.5,5.5],[i-.5,i-.5],'k-','linewidth',2);
        plot([i-.5,i-.5],[.5,5.5],'k-','linewidth',2);
    end
    axis square;
    axis tight;
    set(ax(ifil*2),'fontsize',16,'linewidth',2, ...
        'XTickLabel',{'$\mathbf{X}$','$\mathbf{Y}$','$\mathbf{Z}$','{\boldmath$\tau$}','$\mathbf{V_{P}}$'},'TickLabelInterpreter','latex', ...
        'YTickLabel',{'$\mathbf{X}$','$\mathbf{Y}$','$\mathbf{Z}$','{\boldmath$\tau$}','$\mathbf{V_{P}}$'});
    title(ax(ifil*2),'\textbf{Model Covariance}','fontsize',18,'Interpreter','latex');
%     cb2 = colorbar(ax2);
%     ylabel(cb2,'$\log_{10}(C_m)$','fontsize',16,'Interpreter','latex');
    colormap(cmap)
    caxis([-max(max(abs(Cm_mat{ifil}(:,:,1)))) max(max(abs(Cm_mat{ifil}(:,:,1))))])
    
    
    set(ax(ifil*2-1),'Position',[ax1Pos(1)-0.04 ax1Pos(2:4)]);
    set(ax(ifil*2),'Position',ax2Pos);
    dx = 0.8;
    dy = 0.1;
    text(ax(ifil*2-1),diff(ax(ifil*2-1).XLim)*dx+ax(ifil*2-1).XLim(1),diff(ax(ifil*2-1).YLim)*dy+ax(ifil*2-1).YLim(1),...
    sprintf('\\textbf{ %.4f}',R_spread),'color',[0 0 0],'interpreter','latex','fontsize',14);

end
%% SAVE
if ifsave
    save2pdf(ofile,f906)
end

end