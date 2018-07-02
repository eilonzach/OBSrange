clear; close all;

%% Parameters
% projpath = '/Users/zeilon/Documents/MATLAB/OBS_range_dist/projects/Worcester_OBS-Relocations'; % path to project
projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA/'; % path to project
outdir_OBSrange_nocorr = [projpath,'/OUT_OBSrange_synthsurveys/OUT_nocorr/']; % Path to results WITHOUT velocity corrections
outdir_OBSrange_wcorr = [projpath,'/OUT_OBSrange_synthsurveys/OUT_wcorr_1pts/']; % Path to results WITH velocity corrections

% outdir_OBSrange = './OUT_OBSrange/OUT_Vshipcorrs/'; % output directory
% if ~exist(outdir_OBSrange)
%     mkdir(outdir_OBSrange);
% end

% prepend functions directory to MATLAB path
fullMAINpath = mfilename('fullpath');
functionspath = [fullMAINpath(1:regexp(fullMAINpath,mfilename)-1),'functions'];
addpath(functionspath);

files = dir([outdir_OBSrange_nocorr,'/mats/*.mat']);
Nfils = length(files);
% Load data without corrections
for ifil = 1:Nfils
    filenames_nocorr{ifil} = [outdir_OBSrange_nocorr,'/mats/',files(ifil).name];
end
[ data_all_nocorr,misfit_nocorr,lgd_nocorr,symbol_nocorr ] = agg_synthsurveys( filenames_nocorr,Nfils );
% Load data with corrections
for ifil = 1:Nfils
    filenames_wcorr{ifil} = [outdir_OBSrange_wcorr,'/mats/',files(ifil).name];
end
[ data_all_wcorr,misfit_wcorr,lgd_wcorr,symbol_wcorr ] = agg_synthsurveys( filenames_wcorr,Nfils );
        
%% Plot
% set up figure
f100 = figure(100); clf;
set(gcf,'position',[3 1 1208 697],'color','white');

% set up axes
ax1 = axes('pos',[0.08 0.58 0.25 0.35]); axis equal;
plot(ax1,[0 100],[0 100],'-k','linewidth',2); hold on;
ax2 = axes('pos',[0.38 0.58 0.25 0.35]); axis equal;
plot(ax2,[0 1000],[0 1000],'-k','linewidth',2); hold on;
ax3 = axes('pos',[0.68 0.58 0.25 0.35]); axis equal;
plot(ax3,[0 1000],[0 1000],'-k','linewidth',2); hold on;
ax4 = axes('pos',[0.08 0.12 0.25 0.35]); axis equal;
plot(ax4,[0 1000],[0 1000],'-k','linewidth',2); hold on;
ax5 = axes('pos',[0.38 0.12 0.25 0.35]); axis equal;
plot(ax5,[0 1000],[0 1000],'-k','linewidth',2); hold on;
ax6 = axes('pos',[0.68 0.12 0.25 0.35]); axis equal;
plot(ax6,[0 1000],[0 1000],'-k','linewidth',2); hold on;
clr = [brewermap(7,'Blues'); brewermap(2,'Greens'); brewermap(1,'Reds'); brewermap(2,'Purples'); brewermap(2,'RdPu')];
dx = 0.1;
dy = 0.85;
markersize = 14;
for ifil = 1:Nfils

    % make the plots
    % plot(dtwt_all_nocorr*1000,dtwt_all_wcorr*1000,'ok','markerfacecolor',[0 0 0],'markersize',5); hold on;
    plot(ax1,misfit_nocorr.dtwt_all(ifil)*1000,misfit_wcorr.dtwt_all(ifil)*1000,symbol_nocorr{ifil},'MarkerFaceColor',clr(ifil,:),'Markersize',markersize,'linewidth',1); hold on;
    % plot([-100:100],m*[-100:100]+b,'-r','linewidth',2);
    set(ax1,'fontsize',16,'linewidth',2,'box','on')
    xlabel(ax1,'No Correction','Interpreter','latex');
    ylabel(ax1,'Corrected','Interpreter','latex');
    title(ax1,'$\mathbf{\delta TWT \, (ms)}$','Interpreter','latex');
    xlim(ax1,[0 max(abs(misfit_nocorr.dtwt_all))]*1000);
    ylim(ax1,[0 max(abs(misfit_nocorr.dtwt_all))]*1000);
    colormap(clr);
    ax1.XTick = ax1.YTick;
    if ifil == 1
        [R_dtwt,P_dtwt] = corrcoef(misfit_nocorr.dtwt_all,misfit_wcorr.dtwt_all);
        fitvars = polyfit(misfit_nocorr.dtwt_all*1000, misfit_wcorr.dtwt_all*1000, 1);
        m = fitvars(1);
        b = fitvars(2);
        text(ax1,diff(ax1.XLim)*dx+ax1.XLim(1),diff(ax1.YLim)*dy+ax1.YLim(1),...
        sprintf('\\textbf{ %.2f}',R_dtwt(2,1)),'color',[0 0 0],'interpreter','latex','fontsize',14);
    end

    % plot(dtwt_all_nocorr*1000,dtwt_all_wcorr*1000,'ok','markerfacecolor',[0 0 0],'markersize',5); hold on;
    plot(ax2,misfit_nocorr.r_xy(ifil),misfit_wcorr.r_xy(ifil),symbol_nocorr{ifil},'MarkerFaceColor',clr(ifil,:),'Markersize',markersize,'linewidth',1); hold on;
    % plot([-100:100],m*[-100:100]+b,'-r','linewidth',2);
    set(ax2,'fontsize',16,'linewidth',2,'box','on')
    xlabel(ax2,'No Correction','Interpreter','latex');
    ylabel(ax2,'Corrected','Interpreter','latex');
    title(ax2,'$\mathbf{\delta r_{xy} \, (m)}$','Interpreter','latex');
    xlim(ax2,[0 max(abs(misfit_nocorr.r_xy))]);
    ylim(ax2,[0 max(abs(misfit_nocorr.r_xy))]);
    colormap(clr);
    ax2.XTick = ax2.YTick;
    if ifil == 1
        [R_rxy,P_rxy] = corrcoef(misfit_nocorr.r_xy,misfit_wcorr.r_xy);
        fitvars = polyfit(misfit_nocorr.r_xy, misfit_wcorr.r_xy, 1);
        m = fitvars(1);
        b = fitvars(2);
        text(ax2,diff(ax2.XLim)*dx+ax2.XLim(1),diff(ax2.YLim)*dy+ax2.YLim(1),...
        sprintf('\\textbf{ %.2f}',R_rxy(2,1)),'color',[0 0 0],'interpreter','latex','fontsize',14);
    end
    
    % plot(dtwt_all_nocorr*1000,dtwt_all_wcorr*1000,'ok','markerfacecolor',[0 0 0],'markersize',5); hold on;
    plot(ax3,misfit_nocorr.zsta(ifil),misfit_wcorr.zsta(ifil),symbol_nocorr{ifil},'MarkerFaceColor',clr(ifil,:),'Markersize',markersize,'linewidth',1); hold on;
    % plot([-100:100],m*[-100:100]+b,'-r','linewidth',2);
    set(ax3,'fontsize',16,'linewidth',2,'box','on')
    xlabel(ax3,'No Correction','Interpreter','latex');
    ylabel(ax3,'Corrected','Interpreter','latex');
    title(ax3,'$\mathbf{\delta Z \, (m)}$','Interpreter','latex');
    xlim(ax3,[0 max(abs(misfit_nocorr.zsta))]);
    ylim(ax3,[0 max(abs(misfit_nocorr.zsta))]);
    colormap(clr);
    ax3.XTick = ax3.YTick;
    if ifil == 1
        [R_zsta,P_zsta] = corrcoef(misfit_nocorr.zsta, misfit_wcorr.zsta);
        fitvars = polyfit(misfit_nocorr.zsta, misfit_wcorr.zsta, 1);
        m = fitvars(1);
        b = fitvars(2);
        text(ax3,diff(ax3.XLim)*dx+ax3.XLim(1),diff(ax3.YLim)*dy+ax3.YLim(1),...
        sprintf('\\textbf{ %.2f}',R_zsta(2,1)),'color',[0 0 0],'interpreter','latex','fontsize',14);
    end

    % plot(dtwt_all_nocorr*1000,dtwt_all_wcorr*1000,'ok','markerfacecolor',[0 0 0],'markersize',5); hold on;
    plot(ax4,misfit_nocorr.Vw(ifil),misfit_wcorr.Vw(ifil),symbol_nocorr{ifil},'MarkerFaceColor',clr(ifil,:),'Markersize',markersize,'linewidth',1); hold on;
    % plot([-100:100],m*[-100:100]+b,'-r','linewidth',2);
    set(ax4,'fontsize',16,'linewidth',2,'box','on')
    xlabel(ax4,'No Correction','Interpreter','latex');
    ylabel(ax4,'Corrected','Interpreter','latex');
    title(ax4,'$\mathbf{\delta V_{H_2O} \, (m/s)}$','Interpreter','latex');
    xlim(ax4,[0 max(abs(misfit_nocorr.Vw))]);
    ylim(ax4,[0 max(abs(misfit_nocorr.Vw))]);
    colormap(clr);
    ax4.XTick = ax4.YTick;
    if ifil == 1
        [R_Vw,P_Vw] = corrcoef(misfit_nocorr.Vw, misfit_wcorr.Vw);
        fitvars = polyfit(misfit_nocorr.Vw, misfit_wcorr.Vw, 1);
        m = fitvars(1);
        b = fitvars(2);
        text(ax4,diff(ax4.XLim)*dx+ax4.XLim(1),diff(ax4.YLim)*dy+ax4.YLim(1),...
        sprintf('\\textbf{ %.2f}',R_Vw(2,1)),'color',[0 0 0],'interpreter','latex','fontsize',14);
    end

    % plot(dtwt_all_nocorr*1000,dtwt_all_wcorr*1000,'ok','markerfacecolor',[0 0 0],'markersize',5); hold on;
    plot(ax5,misfit_nocorr.TAT(ifil)*1000,misfit_wcorr.TAT(ifil)*1000,symbol_nocorr{ifil},'MarkerFaceColor',clr(ifil,:),'Markersize',markersize,'linewidth',1); hold on;
    % plot([-100:100],m*[-100:100]+b,'-r','linewidth',2);
    set(ax5,'fontsize',16,'linewidth',2,'box','on')
    xlabel(ax5,'No Correction','Interpreter','latex');
    ylabel(ax5,'Corrected','Interpreter','latex');
    title(ax5,'$\mathbf{\delta TAT \, (ms)}$','Interpreter','latex');
    xlim(ax5,[0 max(abs(misfit_nocorr.TAT))]*1000);
    ylim(ax5,[0 max(abs(misfit_nocorr.TAT))]*1000);
    colormap(clr);
    ax5.XTick = ax5.YTick;
    if ifil == 1
        [R_TAT,P_TAT] = corrcoef(misfit_nocorr.TAT*1000, misfit_wcorr.TAT*1000);
        fitvars = polyfit(misfit_nocorr.TAT*1000, misfit_wcorr.TAT*1000, 1);
        m = fitvars(1);
        b = fitvars(2);
        text(ax5,diff(ax5.XLim)*dx+ax5.XLim(1),diff(ax5.YLim)*dy+ax5.YLim(1),...
        sprintf('\\textbf{ %.2f}',R_TAT(2,1)),'color',[0 0 0],'interpreter','latex','fontsize',14);
    end

%     [R_dtwt,P_dtwt] = corrcoef(misfit_nocorr.dtwt_all,misfit_wcorr.dtwt_all);
%     % plot(dtwt_all_nocorr*1000,dtwt_all_wcorr*1000,'ok','markerfacecolor',[0 0 0],'markersize',5); hold on;
%     plot(ax6,misfit_nocorr.dtwt_all(ifil)*1000,misfit_wcorr.dtwt_all(ifil)*1000,symbol_nocorr{ifil},'MarkerFaceColor',clr(ifil,:),'Markersize',markersize,'linewidth',1); hold on;
%     % plot([-100:100],m*[-100:100]+b,'-r','linewidth',2);
%     set(ax6,'fontsize',16,'linewidth',2,'box','on')
%     xlabel(ax6,'No Correction','Interpreter','latex');
%     ylabel(ax6,'Corrected','Interpreter','latex');
%     title(ax6,'\textbf{TWT Residuals (ms)}','Interpreter','latex');
%     xlim(ax6,[0 max(abs(misfit_nocorr.dtwt_all))]*1000);
%     ylim(ax6,[0 max(abs(misfit_nocorr.dtwt_all))]*1000);
%     colormap(clr);
%     ax6.XTick = ax6.YTick;
%     if ifil == 1
%         fitvars = polyfit(misfit_nocorr.dtwt_all*1000, misfit_wcorr.dtwt_all*1000, 1);
%         m = fitvars(1);
%         b = fitvars(2);
%         text(ax6,diff(ax6.XLim)*dx+ax6.XLim(1),diff(ax6.YLim)*dy+ax6.YLim(1),...
%         sprintf('\\textbf{ %.2f}',m),'color',[0 0 0],'interpreter','latex','fontsize',14);
%     end
end

save2pdf([outdir_OBSrange_wcorr,'analyze_Vshipcorrs.pdf'],f100,500)
%% Plot residuals and Erms
% figure(1); clf;
% set(gcf,'position',[102   348   845   357]);
% dt = 2/1000;
% 
% % Divergent colormap
% % cmap = cmocean('balance');
% dx = 0.1;
% dy = 0.9;
% 
% % ax1 = subplot(1,2,1);
% % [R_dtwt,P_dtwt] = corrcoef(misfit_nocorr.agg.dtwt_all,misfit_wcorr.agg.dtwt_all);
% % plot([-100 100],[-100 100],'-k','linewidth',2); hold on;
% % plot([0 0],[-100 100],'--','color',[0.5 0.5 0.5],'linewidth',2);
% % plot([-100 100],[0 0],'--','color',[0.5 0.5 0.5],'linewidth',2);
% % % plot(dtwt_all_nocorr*1000,dtwt_all_wcorr*1000,'ok','markerfacecolor',[0 0 0],'markersize',5); hold on;
% % scatter(misfit_nocorr.agg.dtwt_all*1000,misfit_wcorr.agg.dtwt_all*1000,60,misfit_wcorr.agg.dtwtcorr_all*1000,'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on;
% % fitvars = polyfit(misfit_nocorr.agg.dtwt_all*1000, misfit_wcorr.agg.dtwt_all*1000, 1);
% % m = fitvars(1);
% % b = fitvars(2);
% % plot([-100:100],m*[-100:100]+c,'-r','linewidth',2);
% % set(gca,'fontsize',16,'linewidth',2,'box','on')
% % axis equal;
% % xlabel('No Correction','Interpreter','latex');
% % ylabel('Corrected','Interpreter','latex');
% % title('\textbf{TWT Residuals (ms)}','Interpreter','latex');
% % xlim([-max(abs(misfit_nocorr.agg.dtwt_all))-dt max(abs(misfit_nocorr.agg.dtwt_all))+dt]*1000);
% % ylim([-max(abs(misfit_nocorr.agg.dtwt_all))-dt max(abs(misfit_nocorr.agg.dtwt_all))+dt]*1000);
% % cb = colorbar('peer',ax1);
% % ylabel(cb,'TWT corrections (ms)','fontsize',16,'Interpreter','latex');
% % colormap(cmap);
% % caxis([-max(abs(misfit_wcorr.agg.dtwtcorr_all)) max(abs(misfit_wcorr.agg.dtwtcorr_all))]*1000);
% % ax1.XTick = ax1.YTick;
% % text(ax1,diff(ax1.XLim)*dx+ax1.XLim(1),diff(ax1.YLim)*dy+ax1.YLim(1),...
% %     sprintf('\\textbf{ %.2f ms}',m),'color',[0 0 0],'interpreter','latex','fontsize',14);
% 
% ax1 = subplot(1,2,1);
% [R_dtwt,P_dtwt] = corrcoef(misfit_nocorr.dtwt_all,misfit_wcorr.dtwt_all);
% plot([-100 100],[-100 100],'-k','linewidth',2); hold on;
% % plot(dtwt_all_nocorr*1000,dtwt_all_wcorr*1000,'ok','markerfacecolor',[0 0 0],'markersize',5); hold on;
% h = scatter(misfit_nocorr.dtwt_all*1000,misfit_wcorr.dtwt_all*1000,60,misfit_wcorr.dtwtcorr_all*1000,'filled','markeredgecolor',[0 0 0],'linewidth',1); hold on;
% fitvars = polyfit(misfit_nocorr.dtwt_all*1000, misfit_wcorr.dtwt_all*1000, 1);
% m = fitvars(1);
% b = fitvars(2);
% % plot([-100:100],m*[-100:100]+b,'-r','linewidth',2);
% set(gca,'fontsize',16,'linewidth',2,'box','on')
% axis equal;
% xlabel('No Correction','Interpreter','latex');
% ylabel('Corrected','Interpreter','latex');
% title('\textbf{TWT Residuals (ms)}','Interpreter','latex');
% xlim([0 max(abs(misfit_nocorr.dtwt_all))+dt]*1000);
% ylim([0 max(abs(misfit_nocorr.dtwt_all))+dt]*1000);
% % cb = colorbar('peer',ax1);
% % ylabel(cb,'TWT corrections (ms)','fontsize',16,'Interpreter','latex');
% colormap(clr);
% % caxis([-max(abs(misfit_wcorr.dtwtcorr_all)) max(abs(misfit_wcorr.dtwtcorr_all))]*1000);
% ax1.XTick = ax1.YTick;
% text(ax1,diff(ax1.XLim)*dx+ax1.XLim(1),diff(ax1.YLim)*dy+ax1.YLim(1),...
%     sprintf('\\textbf{ %.2f}',m),'color',[0 0 0],'interpreter','latex','fontsize',14);
% 
% ax2 = subplot(1,2,2);
% [R_Erms,P_Erms] = corrcoef(misfit_nocorr.E_rms,misfit_wcorr.E_rms);
% plot(misfit_nocorr.E_rms*1000,misfit_wcorr.E_rms*1000,'ok','markerfacecolor',[0.5 0.5 0.5],'markersize',8); hold on;
% plot([-100 100],[-100 100],'-k','linewidth',2);
% set(gca,'fontsize',16,'linewidth',2,'box','on')
% axis equal;
% xlabel('No Correction','Interpreter','latex');
% ylabel('Corrected','Interpreter','latex');
% title('\textbf{RMS (ms)}','Interpreter','latex');
% xlim([0 max(misfit_nocorr.E_rms)+dt]*1000);
% ylim([0 max(misfit_nocorr.E_rms)+dt]*1000);
% ax2.XTick = ax2.YTick;
% 
% % resize axis 1
% ax1Pos = get(ax1,'position');
% ax2Pos = get(ax2,'position');
% ax1Pos(3:4) = [ax2Pos(3:4)];
% ax1Pos(1) = ax1Pos(1)-0.04;
% set(ax1,'position',ax1Pos);
% set(ax2,'position',[ax2Pos(1)+0.04,ax2Pos(2:4)]);
% 
% save2pdf([outdir_OBSrange_wcorr,'analyze_Vshipcorrs.pdf'],1,500)

%% print to csv file:
% fid = fopen([projname,'_allstas.csv'],'w');
% fprintf(fid,'#, STA,  LAT(deg.),  LON(deg.),  depth(m), driftX, driftY, TAT, V_water, Erms  \n');
% for is = 1:length(deployarray.slats)
%     if deployarray.selevs(is)==0, depstr = ' '; else, depstr = sprintf('%4.0f',-1000*deployarray.selevs(is)); end
%     fprintf(fid,'%2.0f, %s, %2.0f%s%06.3f%sS, %4.0f%s%06.3f%sW, %9.5f, %9.5f,   %s,    %6.1f,      %6.1f \n',...
%         deployarray.order(is),deployarray.stas{is},...
%         abs(fix(deployarray.slats(is))),'  ',60*abs(rem(deployarray.slats(is),1)),39,...
%         abs(fix(deployarray.slons(is))),'  ',60*abs(rem(deployarray.slons(is),1)),39,...
%         [deployarray.slats(is),deployarray.slons(is)],...
%         depstr,deployarray.nm_from_last(is),deployarray.nm_to_next(is));
% end
% fclose(fid);


