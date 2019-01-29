%% Function to produce Figure 07
%  Figure 07 Resolution and covariance matrices
function Figure07

ofile = '../Figure07';
ifsave = 1;

projname = 'mats_SynthBoot_summary_noTAT';

files = {
        'SynthBoot_PACMAN_rad1.00_OUT.mat';
        'SynthBoot_circle_rad1.00_OUT.mat';
        'SynthBoot_line_rad1.00_OUT.mat';
         };
     
titles = {
    'PACMAN 1 Nm';
    'Circle 1 Nm';
    'Line 1 Nm';
    };

%% Load *.mat files
% files = dir('../figdata/mats_SynthBoot_summary/*.mat');

Nfils = length(files);
for ifil = 1:Nfils
    load(['../figdata/',projname,'/',files{ifil}]);
    
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
    
    C_mat{ifil} = cov2cor(Cm_mat{ifil});
end

%% Model resolution and Covariance
    f907 = figure(907); clf;
    set(gcf,'position',[248   322   607   365]);
    cmap = cmocean('balance');
%     cmap_R = cmocean('tempo');
    cmap_R = flip(gray);
    
for ifil = 1:length(files)  
    %%%%%%%%% RESOLUTION
    ax(ifil*2-1) = subplot(2,length(files),ifil);
    ax1Pos = ax(ifil*2-1).Position;
    imagesc(ax(ifil*2-1),abs(R_mat{ifil}(:,:,1))); hold on;
    for i = 1:4
        plot([.5,4.5],[i-.5,i-.5],'k-','linewidth',1.5);
        plot([i-.5,i-.5],[.5,4.5],'k-','linewidth',1.5);
    end
    R_spread = sum(sum((R_mat{ifil}(:,:,1)-eye(size(R_mat{ifil}(:,:,1)))).^2));
%     R_spread = norm(R_mat{ifil}(:,:,1)-eye(size(R_mat{ifil}(:,:,1))),2).^2;
    axis square;
    axis tight;
    xticks([1:4]);
    if ifil == 1
        set(ax(ifil*2-1),'fontsize',16,'linewidth',3,'TickLength',[0 0], ...
            'XTickLabel',[],'TickLabelInterpreter','latex', ...
            'YTickLabel',{'$\mathbf{X}$','$\mathbf{Y}$','$\mathbf{Z}$','$\mathbf{V_{P}}$'});
        
        ylabel(ax(ifil*2-1),'\textbf{Resolution}','fontsize',16,'Interpreter','latex');
    else
        set(ax(ifil*2-1),'fontsize',16,'linewidth',3,'TickLength',[0 0], ...
            'XTickLabel',[],'TickLabelInterpreter','latex', ...
            'YTickLabel',[]);
    end
    t = title(ax(ifil*2-1),sprintf('\\textbf{%s}',titles{ifil}),'fontsize',16,'Interpreter','latex');
    if ifil == length(files)
        cb1 = colorbar(ax(ifil*2-1),'LineWidth',2);
%         cb1.Position = [cb1.Position(1) cb1.Position(2)-0.3 cb1.Position(3) cb1.Position(4)];
    end
    colormap(ax(ifil*2-1),cmap_R)
    caxis(ax(ifil*2-1),[0  1]);
    
    %%%%%%%%%% COVARIANCE
    ax(ifil*2) = subplot(2,length(files),ifil+length(files));
    ax2Pos = ax(ifil*2).Position;
%     imagesc(ax2,sign(Cm_mat{ifil}(:,:,1)).*log10(abs(Cm_mat{ifil}(:,:,1)))); hold on;
    imagesc(ax(ifil*2),C_mat{ifil}(:,:,1)); hold on;
    for i = 1:4
        plot([.5,4.5],[i-.5,i-.5],'k-','linewidth',1.5);
        plot([i-.5,i-.5],[.5,4.5],'k-','linewidth',1.5);
    end
    axis square;
    axis tight;
    xticks([1:4]);
    if ifil == 1
        set(ax(ifil*2),'fontsize',16,'linewidth',3,'TickLength',[0 0], ...
            'XTickLabel',{'$\mathbf{X}$','$\mathbf{Y}$','$\mathbf{Z}$','$\mathbf{V_{P}}$'},'TickLabelInterpreter','latex', ...
            'YTickLabel',{'$\mathbf{X}$','$\mathbf{Y}$','$\mathbf{Z}$','$\mathbf{V_{P}}$'});
        
        ylabel(ax(ifil*2),'\textbf{Correlation}','fontsize',16,'Interpreter','latex');
    else
        set(ax(ifil*2),'fontsize',16,'linewidth',3,'TickLength',[0 0], ...
            'XTickLabel',{'$\mathbf{X}$','$\mathbf{Y}$','$\mathbf{Z}$','$\mathbf{V_{P}}$'},'TickLabelInterpreter','latex', ...
            'YTickLabel',[]);
    end
    if ifil == length(files)
        cb2 = colorbar(ax(ifil*2),'linewidth',2);
    end
    colormap(ax(ifil*2),cmap)
    caxis(ax(ifil*2),[-1 1]);
    
    
    %% positioning...
    
    % Windows
    yfac = 1.2;
    dx = -0.04; %+0.02;
    dx2 = -0.01;
    dy = -0.08;
    set(ax(ifil*2-1),'Position',[ax1Pos(1)+dx+dx2*ifil ax1Pos(2)+dy ax1Pos(3)*yfac ax1Pos(4)*yfac]);
    dy = -0.05;
    set(ax(ifil*2),'Position',[ax2Pos(1)+dx+dx2*ifil ax2Pos(2)+dy ax2Pos(3)*yfac ax2Pos(4)*yfac]);
    
    
    % Title positioning
    t.Position = [t.Position(1) t.Position(2)-0.2 t.Position(3)];
    
    % Text for "spread"
    dx = 0.65;
    dy = 0.1;
%     text(ax(ifil*2-1),diff(ax(ifil*2-1).XLim)*dx+ax(ifil*2-1).XLim(1),diff(ax(ifil*2-1).YLim)*dy+ax(ifil*2-1).YLim(1),...
%     sprintf('\\textbf{ %.4f}',R_spread),'color',[0 0 0],'interpreter','latex','fontsize',14);
    axpos = ax(ifil*2-1).Position;
    annotation('textbox',[axpos(1)+0.17 axpos(2)+0.35 0.07 0.04],...
    'string',{sprintf('%.4f',R_spread)},'color',[0 0 0],'BackgroundColor',cmap_R(1,:),...
    'EdgeColor','none','fontsize',14,'fontweight','bold',...
    'VerticalAlignment','middle','HorizontalAlignment','center');
    
end
cb1.Position = [cb1.Position(1)+0.01 cb1.Position(2) cb1.Position(3) cb1.Position(4)];
cb2.Position = [cb2.Position(1)+0.01 cb2.Position(2) cb2.Position(3) cb2.Position(4)];

%% SAVE
if ifsave
    save2pdf(ofile,f907)
end

end