%% Function to produce Figure 05 
%  Figure 05 compares our algorithm to SIO
function Figure05_test_boot_jr

% ofile = '../Figure05';
ifsave = 0;

%% load 

synth_dirs = {
    '1_OUT_nocorr';
    '2_OUT_wcorr_xrec';
    '3_OUT_wcorr_xrec_TAT';
    '4_OUT_wcorr_xrec_Vp';
    '5_OUT_wcorr_xrec_Z';
    '6_OUT_wcorr_xrec_TAT_Vp_Z';
    };

xlabels = {
    'No Doppler Corr.';
    'Doppler Corr.';
    'No TAT';
    'No V_p';
    'No Z';
    'No TAT, V_p, Z';
    };

synth_path = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_synthtestboot/OUT_OBSrange/';

% Load synthetic
for ifil = 1:length(synth_dirs)
    synth_mat = dir(fullfile(synth_path,synth_dirs{ifil},'mats/*.mat'));
    synth = load(fullfile(synth_path,synth_dirs{ifil},'mats',synth_mat.name));
%     misfit_xsta(ifil) = mean(synth.datamat.x_sta_bs) - trudata.obs_location_xyz(1)*1000;
%     misfit_ysta(ifil) = mean(synth.datamat.y_sta_bs) - trudata.obs_location_xyz(2)*1000;
%     misfit_zsta(ifil) = mean(synth.datamat.z_sta_bs) - (-trudata.obs_location_xyz(3)*1000);
%     if ~isempty(data.datamat.E_rms)
%         misfit_r_xy(ifil) = mean(synth.datamat.E_rms);
%     else
%         misfit_r_xy(ifil) = rms(synth.datamat.dtwt_bs);
%     end

%     misfit_xsta_bs(:,ifil) = synth.datamat.x_sta_bs - trudata.obs_location_xyz(1)*1000;
%     misfit_ysta_bs(:,ifil) = synth.datamat.y_sta_bs - trudata.obs_location_xyz(2)*1000;
%     misfit_zsta_bs(:,ifil) = synth.datamat.z_sta_bs - (-trudata.obs_location_xyz(3)*1000);
%     misfit_TAT_bs(:,ifil) = synth.datamat.TAT_bs - trudata.tat;
%     misfit_Vp_bs(:,ifil) = synth.datamat.V_w_bs - trudata.vp_actual;
%     misfit_r_xy_bs(:,ifil) = sqrt( misfit_xsta_bs(:,ifil).^2 + misfit_ysta_bs(:,ifil).^2 );
%     misfit_r_xy_bs(:,ifil) = sqrt( (misfit_xsta_bs(:,ifil)./trudata.obs_location_xyz(1)*1000).^2 +...
%                                    (misfit_ysta_bs(:,ifil)./trudata.obs_location_xyz(2)*1000).^2 +...
%                                    (misfit_zsta_bs(:,ifil)./(-trudata.obs_location_xyz(3)*1000)).^2 +...
%                                    (misfit_TAT_bs(:,ifil)./trudata.tat).^2 +...
%                                    (misfit_Vp_bs(:,ifil)./trudata.vp_actual).^2);
%     misfit_r_xy(ifil) = mean(misfit_r_xy_bs(:,ifil)); 
    misfit_r_xy_bs(:,ifil) = synth.data(1).misfit_r_xy;
    misfit_r_xy(ifil) = rms(synth.data(1).misfit_r_xy);
    r_xy_95(:,ifil) = abs(mean(synth.data(1).misfit_r_xy) - prctile(synth.data(1).misfit_r_xy,[2.5 97.5]));
end
r_xy_95(:,8) = [nan nan]';
r_xy_95(:,9) = [nan nan]';


%% ---------------------   PLOTTING   ---------------------   

%% Plot comparisons
f905 = figure(906); clf;
set(f905,'position',[159   165   552   540]);
ax1 = subplot(2,1,1); hold(ax1,'on'); box on;
ax2 = subplot(2,1,2); hold(ax2,'on'); box on;
dy = 1.1;
dy_space = 1.6;
ax1.Position = [ax1.Position(1), ax1.Position(2), ax1.Position(3), ax1.Position(4)*dy];
ax2.Position = [ax2.Position(1), ax2.Position(2)*dy_space, ax2.Position(3), ax2.Position(4)*dy];

markersize = 14;
%clr = parula(Nfils);
% clr = [brewermap(7,'Blues'); brewermap(3,'Greys'); brewermap(2,'Greens'); brewermap(1,'Reds'); brewermap(2,'Purples'); brewermap(2,'RdPu')];
% clr = [repmat([0.5 0.5 0.5],7,1); [0.8 0 0]; [0.8 0 0] ];
clr = [repmat([0.7 0.7 0.7],length(synth_dirs),1)];

for ifil = 1:length(synth_dirs)
    h(ifil) = plot(ax2,ifil,misfit_r_xy(ifil),'sk','markerfacecolor',clr(ifil,:),'markersize',markersize); hold on;
    errorbar(ax2,ifil,misfit_r_xy(ifil),r_xy_95(1,ifil),r_xy_95(2,ifil),'.k','markerfacecolor',[0.5 0.5 0.5],'markersize',markersize,'linewidth',1.5); hold on;
    set(ax2,'yscale','log','linewidth',1.5,'fontsize',16);
    ylabel(ax2,'$\mathbf{\delta r_{xy}\, (m)}$','fontsize',18,'Interpreter','latex')
    yticks(ax2,[0.1 1 10 100]);
    xticks(ax2,[1:9]);
    xticklabels(ax2,xlabels);
    xtickangle(ax2,45);
    xlim(ax2,[0 length(synth_dirs)+1]);   
%     ylim(ax2,[0.1 max(misfit_r_xy)+10^(floor(log10(max(misfit_r_xy))))]);
    ylim(ax2,[1 100]);
end

% l = legend(h,lgd,'position',[0.8 0.05 0.1852 0.95],'interpreter','none','fontsize',14,'box','off');

%% Plot histograms
figure(2);

for ii = 1:2
    plot_hist(subplot(1,2,ii),misfit_r_xy_bs(:,ii),10)
    xlabel('\delta r_{xy}');
end

%% SAVE
if ifsave
    save2pdf(ofile,f905)
end

end
