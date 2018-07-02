clear; close all;

%% Parameters
% projpath = '/Users/zeilon/Documents/MATLAB/OBS_range_dist/projects/Worcester_OBS-Relocations'; % path to project
projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA/'; % path to project
outdir_OBSrange_nocorr = [projpath,'/OUT_OBSrange/OUT_nocorr/']; % Path to results WITHOUT velocity corrections
outdir_OBSrange_wcorr = [projpath,'/OUT_OBSrange/OUT_wcorr_1pts/']; % Path to results WITH velocity corrections

% outdir_OBSrange = './OUT_OBSrange/OUT_Vshipcorrs/'; % output directory
% if ~exist(outdir_OBSrange)
%     mkdir(outdir_OBSrange);
% end


projname =fliplr(strtok(fliplr(projpath),'/'));
wd = pwd;
% cd(projpath);cd(outdir_OBSrange);
files = dir([outdir_OBSrange_nocorr,'mats']);
Nfiles = length(files);
Nstas = 0;
for ifile = 1:Nfiles
    if any(regexp(files(ifile).name,'mat'))
        Nstas = Nstas+1;
        nocorr = load([outdir_OBSrange_nocorr,'mats/',files(ifile).name]);
        wcorr = load([outdir_OBSrange_wcorr,'mats/',files(ifile).name]);
        if Nstas==1 
            allstas_nocorr = nocorr.datamat;
            allstas_wcorr = wcorr.datamat;
        else
            allstas_nocorr(Nstas) = nocorr.datamat;
            allstas_wcorr(Nstas) = wcorr.datamat;
        end
        Erms_all_nocorr(Nstas) = mean(allstas_nocorr(Nstas).E_rms);
        Erms_all_wcorr(Nstas) = mean(allstas_wcorr(Nstas).E_rms);
        Erms_dtwt_nocorr(Nstas) = rms(allstas_nocorr(Nstas).dtwt_bs);
        Erms_dtwt_wcorr(Nstas) = rms(allstas_wcorr(Nstas).dtwt_bs);
    end
end

dtwt_all_nocorr = vertcat(allstas_nocorr.dtwt_bs);
dtwt_all_wcorr = vertcat(allstas_wcorr.dtwt_bs);
dtwtcorr_all = vertcat(allstas_wcorr.dtwtcorr_bs);
        
%% Plot residuals and Erms
figure(1); clf;
set(gcf,'position',[102   348   845   357]);
dt = 2/1000;

% Divergent colormap
% cmap = cmocean('balance');
cmap = flip(brewermap(100,'RdBu'));
dx = 0.1;
dy = 0.85;

ax1 = subplot(1,2,1);
[R_dtwt,P_dtwt] = corrcoef(dtwt_all_nocorr,dtwt_all_wcorr);
plot([-100 100],[-100 100],'-k','linewidth',2); hold on;
plot([0 0],[-100 100],'--','color',[0.5 0.5 0.5],'linewidth',2);
plot([-100 100],[0 0],'--','color',[0.5 0.5 0.5],'linewidth',2);
% plot(dtwt_all_nocorr*1000,dtwt_all_wcorr*1000,'ok','markerfacecolor',[0 0 0],'markersize',5); hold on;
scatter(dtwt_all_nocorr*1000,dtwt_all_wcorr*1000,60,dtwtcorr_all*1000,'o','filled','markeredgecolor',[0 0 0],'linewidth',1); hold on;
set(gca,'fontsize',16,'linewidth',2,'box','on')
axis equal;
xlabel('No Correction','Interpreter','latex');
ylabel('Corrected','Interpreter','latex');
title('\textbf{TWT Residuals (ms)}','Interpreter','latex');
xlim([-max(abs(dtwt_all_wcorr))-dt max(abs(dtwt_all_wcorr))+dt]*1000);
ylim([-max(abs(dtwt_all_wcorr))-dt max(abs(dtwt_all_wcorr))+dt]*1000);
cb = colorbar('peer',ax1);
ylabel(cb,'TWT corrections (ms)','fontsize',16,'Interpreter','latex');
colormap(cmap);
caxis([-max(abs(dtwtcorr_all)) max(abs(dtwtcorr_all))]*1000);
ax1.XTick = ax1.YTick;
text(ax1,diff(ax1.XLim)*dx+ax1.XLim(1),diff(ax1.YLim)*dy+ax1.YLim(1),...
        sprintf('\\textbf{ %.2f}',R_dtwt(2,1)),'color',[0 0 0],'interpreter','latex','fontsize',14);

ax2 = subplot(1,2,2);
[R_Erms,P_Erms] = corrcoef(Erms_all_nocorr,Erms_all_wcorr);
plot(Erms_all_nocorr*1000,Erms_all_wcorr*1000,'ok','markerfacecolor',[0.5 0.5 0.5],'markersize',8); hold on;
plot([-100 100],[-100 100],'-k','linewidth',2);
set(gca,'fontsize',16,'linewidth',2,'box','on')
axis equal;
xlabel('No Correction','Interpreter','latex');
ylabel('Corrected','Interpreter','latex');
title('\textbf{RMS (ms)}','Interpreter','latex');
xlim([0 max(Erms_all_wcorr)+dt]*1000);
ylim([0 max(Erms_all_wcorr)+dt]*1000);
ax2.XTick = ax2.YTick;
text(ax2,diff(ax2.XLim)*dx+ax2.XLim(1),diff(ax2.YLim)*dy+ax2.YLim(1),...
        sprintf('\\textbf{ %.2f}',R_Erms(2,1)),'color',[0 0 0],'interpreter','latex','fontsize',14);

% resize axis 1
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax1Pos(3:4) = [ax2Pos(3:4)];
ax1Pos(1) = ax1Pos(1)-0.04;
set(ax1,'position',ax1Pos);
set(ax2,'position',[ax2Pos(1)+0.04,ax2Pos(2:4)]);

save2pdf([outdir_OBSrange_wcorr,'analyze_Vshipcorrs.pdf'],1,500)

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


