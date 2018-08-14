%% Model resolution and Covariance
    f103 = figure(103); clf;
    set(gcf,'position',[91   260   829   438]);
    cmap = cmocean('balance');
    cmap_R = flip(gray);
    C_mat = cov2cor(Cm_mat(:,:,1));
    
    ax1 = subplot(1,2,1);
    ax1Pos = ax1.Position;
    imagesc(ax1,R_mat(:,:,1)); hold on;
    for i = 1:5
        plot([.5,5.5],[i-.5,i-.5],'k-','linewidth',2);
        plot([i-.5,i-.5],[.5,5.5],'k-','linewidth',2);
    end
    R_spread = sum(sum((R_mat(:,:,1)-eye(size(R_mat(:,:,1)))).^2));
    axis square;
    axis tight;
    set(ax1,'fontsize',16,'linewidth',2, ...
        'XTickLabel',{'$\mathbf{X}$','$\mathbf{Y}$','$\mathbf{Z}$','{\boldmath$\tau$}','$\mathbf{V_{P}}$'},'TickLabelInterpreter','latex', ...
        'YTickLabel',{'$\mathbf{X}$','$\mathbf{Y}$','$\mathbf{Z}$','{\boldmath$\tau$}','$\mathbf{V_{P}}$'});
    title(ax1,'\textbf{Model Resolution}','fontsize',18,'Interpreter','latex');
    cb1 = colorbar(ax1);
    colormap(ax1,cmap_R)
%     caxis([-max(max(abs(mean(RR,3)))) max(max(abs(mean(RR,3))))]);
    caxis([0 1]);
    
    ax2 = subplot(1,2,2);
    ax2Pos = ax2.Position;
    imagesc(ax2,C_mat(:,:,1)); hold on;
    for i = 1:5
        plot([.5,5.5],[i-.5,i-.5],'k-','linewidth',2);
        plot([i-.5,i-.5],[.5,5.5],'k-','linewidth',2);
    end
    axis square;
    axis tight;
    set(ax2,'fontsize',16,'linewidth',2, ...
        'XTickLabel',{'$\mathbf{X}$','$\mathbf{Y}$','$\mathbf{Z}$','{\boldmath$\tau$}','$\mathbf{V_{P}}$'},'TickLabelInterpreter','latex', ...
        'YTickLabel',{'$\mathbf{X}$','$\mathbf{Y}$','$\mathbf{Z}$','{\boldmath$\tau$}','$\mathbf{V_{P}}$'});
    title(ax2,'\textbf{Model Covariance}','fontsize',18,'Interpreter','latex');
    cb2 = colorbar(ax2);
    colormap(ax2,cmap)
    caxis([-1 1])
    
%     tk = logspace(0,8,20);
%     cmap = parula(19);
%     ctable = [tk(1:end-1)' cmap*255 tk(2:end)' cmap*255];
% %     cptcmap(ctable, 'mapping', 'direct');
%     colormap(ctable)
% %     cptcbar(ax2, ctable, 'eastoutside', false);
    
    set(ax1,'Position',[ax1Pos(1)-0.04 ax1Pos(2:4)]);
    set(ax2,'Position',ax2Pos);
    dx = 0.8;
    dy = 0.1;
    text(ax1,diff(ax1.XLim)*dx+ax1.XLim(1),diff(ax1.YLim)*dy+ax1.YLim(1),...
    sprintf('\\textbf{ %.4f}',R_spread),'color',[0 0 0],'interpreter','latex','fontsize',14);