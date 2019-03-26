% Save synthetic bootstrap results in smaller mat files for making figures
% for paper.
%

clear; close all;

%% INPUTS - MAKE SURE THESE ARE 
% path to project
% projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_SynthBoot_surveys_strtmodel/'; 
% projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_SynthBoot_surveys_noTAT/'; 
% path to survey data from the project directory
% datapath = '/Users/russell/Lamont/PROJ_OBSrange/synth_tests_paper/synth_surveys/'; 
% path to output directory from project directory (actually input for this script)
% outdir = [projpath,'/OUT_OBSrange_synthsurveys/OUT_wcorr/']; %

% projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_SynthBoot_surveys_noTAT_REVISION1/'; % Josh PAPER
% % projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_SynthBoot_surveys_noTAT_REVISION1_fixtat/'; % Josh PAPER
% datapath = '/Users/russell/Lamont/PROJ_OBSrange/synth_tests_paper/synth_surveys_paper_REVISION1/'; % Josh
% outdir = [projpath,'/OUT_OBSrange_synthsurveys_REVISION1_noGPS/OUT_wcorr/'];

projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_SynthBoot_REVISION2/'; % Josh PAPER
datapath = '/Users/russell/Lamont/PROJ_OBSrange/synth_tests_paper/synth_surveys_paper_REVISION1/'; % Josh
% outdir = './OUT_OBSrange_synthsurveys_REVISION2_testCovm/';
outdir = [projpath,'/OUT_OBSrange_synthsurveys_REVISION2_testAzimuthal/OUT_wcorr/'];
% Put a string survtion name here to only consider that survtion. 
% Otherwise, to locate all survtions, put ''
surv = 'SynthBoot_PACMAN_rad1.00_z5000m_fr10'; %'SynthBoot_PACMAN_rad1.00_z5000m_fr10'; %'SynthBoot_circle_rad1.00'; %'SynthBoot_line_rad1.00'; %'SynthBoot_PACMAN_rad1.00';
% surv_suffix = 'z5000m_fr10';

% prepend functions directory to MATLAB path
fullMAINpath = mfilename('fullpath');
functionspath = [fullMAINpath(1:regexp(fullMAINpath,mfilename)-1),'functions'];
addpath(functionspath);

%% Load *.mat files
% files = dir([outdir,'/mats/*.mat']);
files = dir([outdir,'/mats/',surv,'_OUT.mat']);

Nfils = length(files);
for ifil = 1:Nfils
    display(files(ifil).name)
    load([outdir,'/mats/',files(ifil).name]);
    
    lat_drop = data(1).drop(1);
    lon_drop = data(1).drop(2);
    z_drop = data(1).drop(3)*-1000;
    lats_ship = data(1).survlats;
    lons_ship = data(1).survlons;
    olon = lon_drop;
    olat = lat_drop;
    [ x_ship, y_ship ] = lonlat2xy_nomap( olon, olat, lons_ship, lats_ship );
    
    
    for ii = 1:length(data)
        x_sta(ii) = data(ii).x_sta; x_sta1(ii) = data(ii).x_sta1;
        y_sta(ii) = data(ii).y_sta; y_sta1(ii) = data(ii).y_sta1;
        z_sta(ii) = data(ii).z_sta; z_sta1(ii) = data(ii).z_sta1;
        V_w(ii) = data(ii).V_w;       V_w1(ii) = data(ii).V_w1;
        drift(ii) = data(ii).drift;
        
        x_sta_tru(ii) = data(ii).obs_loc_xyz(1)*1000;
        y_sta_tru(ii) = data(ii).obs_loc_xyz(2)*1000;
        z_sta_tru(ii) = -data(ii).obs_loc_xyz(3)*1000;
        V_w_tru(ii) = data(ii).Vp_water*1000;
        drift_tru(ii) = data(ii).drift_true;
        
        misfit_xsta(ii) = x_sta(ii) - x_sta_tru(ii);
        misfit_ysta(ii) = y_sta(ii) - y_sta_tru(ii);
        misfit_zsta(ii) = z_sta(ii) - z_sta_tru(ii);
        misfit_r_xy(ii) = sqrt( misfit_xsta(ii).^2 + misfit_ysta(ii).^2 );
        misfit_r_xyz(ii) = sqrt( misfit_xsta(ii).^2 + misfit_ysta(ii).^2 + misfit_zsta(ii).^2 );
        misfit_Vw(ii) = V_w(ii) - V_w_tru(ii);
        misfit_drift(ii) = drift(ii)-drift_tru(ii);
        
        misfit_xsta1(ii) = x_sta1(ii) - x_sta_tru(ii);
        misfit_ysta1(ii) = y_sta1(ii) - y_sta_tru(ii);
        misfit_zsta1(ii) = z_sta1(ii) - z_sta_tru(ii);
        misfit_r_xy1(ii) = sqrt( misfit_xsta1(ii).^2 + misfit_ysta1(ii).^2 );
        misfit_r_xyz1(ii) = sqrt( misfit_xsta1(ii).^2 + misfit_ysta1(ii).^2 + misfit_zsta1(ii).^2 );
        misfit_Vw1(ii) = V_w1(ii) - V_w_tru(ii);
        
        azi_locs{ii} = data(ii).azi_locs;
        azi_sort{ii} = sort(azi_locs{ii});
        azi_gaps{ii} = diff(azi_sort{ii});
        azi_gaps{ii} = [azi_gaps{ii}; azi_sort{ii}(end)-azi_sort{ii}(1)];
        azi_gaps{ii}(azi_gaps{ii}>180) = 360 - azi_gaps{ii}(azi_gaps{ii}>180);
        azi_maxgap(ii) = max(azi_gaps{ii});
        azi_Nge10(ii) = sum(azi_gaps{ii}>=10);
        azi_Nge15(ii) = sum(azi_gaps{ii}>=15);
        azi_Nge20(ii) = sum(azi_gaps{ii}>=20);
        azi_Nge30(ii) = sum(azi_gaps{ii}>=30);
        azi_Nge45(ii) = sum(azi_gaps{ii}>=45);
        azi_Nge90(ii) = sum(azi_gaps{ii}>=90);
        azi_Nge180(ii) = sum(azi_gaps{ii}>=180);
        azi_meangap(ii) = mean(azi_gaps{ii});
        azi_medgap(ii) = median(azi_gaps{ii});
        azi_stdgap(ii) = std(azi_gaps{ii});
        
        % Covariance
        Cm_mat(:,:,ii) = data(ii).Cm_mat;
        covm_95 = sqrt(diag(Cm_mat(:,:,ii)))*2;
        x_std_covm(ii) = covm_95(1);
        y_std_covm(ii) = covm_95(2);
        z_std_covm(ii) = covm_95(3);
        vw_std_covm(ii) = covm_95(4);
        
        x_std_bs(ii) = data(ii).x_sta_std*2; 
        y_std_bs(ii) = data(ii).y_sta_std*2; 
        z_std_bs(ii) = data(ii).z_sta_std*2; 
        vw_std_bs(ii) = data(ii).V_w_std*2;  
    end
    
    
%     for ii = 1:length(data)
%         R_mats(:,:,ii) = data(ii).R_mat;
%         Cm_mats(:,:,ii) = data(ii).Cm_mat;
%         
%         
%     end % loop on synth stations
%     data_summary.R_mat = mean(R_mats,3);
%     data_summary.Cm_mat = mean(Cm_mats,3);
    
    %% Save summary mat files
%     if ~exist([outdir,'/mats_SynthBoot_summary/'])
%         mkdir([outdir,'/mats_SynthBoot_summary/']);
%     end   
%     save([outdir,'/mats_SynthBoot_summary/',files(ifil).name],'data_summary')
%     

    C = corrcoef(x_std_bs,x_std_covm);
    R2x = C(1,2)^2;
    C = corrcoef(y_std_bs,y_std_covm);
    R2y = C(1,2)^2;
    C = corrcoef(z_std_bs,z_std_covm);
    R2z = C(1,2)^2;
    C = corrcoef(vw_std_bs,vw_std_covm);
    R2vw = C(1,2)^2;
end

%% PLOT AZIMUTHAL
f1 = figure(1); clf;
set(gcf,'Position',[235     1   614   704]);

clr = [0.5 0.5 0.5];
subplot(3,2,1);
plot(azi_maxgap,abs(misfit_drift),'.','color',clr,'linewidth',2); hold on;
% plot(azi_maxgap,misfit_r_xy1,'o'); hold on;
y = abs(misfit_drift);
drbin = [15:10:65;25:10:75]; clear('mfdr')
for idr = 1:size(drbin,2), mfdr(idr) = median(y([azi_maxgap]<= drbin(2,idr) & [azi_maxgap]> drbin(1,idr))); end
plot(mean(drbin,1),mfdr,'o-k','linewidth',2,'markerfacecolor',[0 0 0]);
title('Maximum Azimuthal Gap ({\circ})');
ylabel('|\delta r_{xy}| (m)');
set(gca,'fontsize',15,'linewidth',1.5);
axis tight;
ylim([0 max(abs(misfit_drift))]);

subplot(3,2,2);
plot(azi_medgap,abs(misfit_drift),'.','color',clr,'linewidth',2); hold on;
% plot(azi_Nge45,misfit_r_xy1,'o'); hold on;
y = abs(misfit_drift);
drbin = [5:0.25:6 ; 5.25:0.25:6.25]; clear('mfdr')
for idr = 1:size(drbin,2), mfdr(idr) = median(y([azi_medgap]<= drbin(2,idr) & [azi_medgap]> drbin(1,idr))); end
% plot(mean(drbin,1),mfdr,'o-k','linewidth',2,'markerfacecolor',[0 0 0]);
title('Median Azimuthal Gap ({\circ})');
ylabel('|\delta r_{xy}| (m)');
set(gca,'fontsize',15,'linewidth',1.5);
axis tight;
ylim([0 max(abs(misfit_drift))]);
% xlim([4 7]);

subplot(3,2,3);
plot(azi_maxgap,abs(misfit_zsta),'.','color',clr,'linewidth',2); hold on;
% plot(azi_maxgap,misfit_r_xy1,'o'); hold on;
% xlabel('Maximum Azimuthal Gap ({\circ})');
y = abs(misfit_zsta);
drbin = [15:10:65;25:10:75]; clear('mfdr')
for idr = 1:size(drbin,2), mfdr(idr) = median(y([azi_maxgap]<= drbin(2,idr) & [azi_maxgap]> drbin(1,idr))); end
plot(mean(drbin,1),mfdr,'o-k','linewidth',2,'markerfacecolor',[0 0 0]);
ylabel('|\delta Z| (m)');
set(gca,'fontsize',15,'linewidth',1.5);
axis tight;
ylim([0 max(abs(misfit_zsta))]);

subplot(3,2,4);
plot(azi_medgap,abs(misfit_zsta),'.','color',clr,'linewidth',2); hold on;
% plot(azi_Nge45,misfit_r_xy1,'o'); hold on;
% xlabel('N_{gaps} >= 20{\circ}');
y = abs(misfit_zsta);
drbin = [5:0.25:6 ; 5.25:0.25:6.25]; clear('mfdr')
for idr = 1:size(drbin,2), mfdr(idr) = median(y([azi_medgap]<= drbin(2,idr) & [azi_medgap]> drbin(1,idr))); end
% plot(mean(drbin,1),mfdr,'o-k','linewidth',2,'markerfacecolor',[0 0 0]);
ylabel('|\delta Z| (m)');
set(gca,'fontsize',15,'linewidth',1.5);
axis tight;
ylim([0 max(abs(misfit_zsta))]);
% xlim([4 7]);

subplot(3,2,5);
plot(azi_maxgap,abs(misfit_Vw),'.','color',clr,'linewidth',2); hold on;
% plot(azi_maxgap,misfit_r_xy1,'o'); hold on;
% xlabel('Maximum Azimuthal Gap ({\circ})');
y = abs(misfit_Vw);
drbin = [15:10:65;25:10:75]; clear('mfdr')
for idr = 1:size(drbin,2), mfdr(idr) = median(y([azi_maxgap]<= drbin(2,idr) & [azi_maxgap]> drbin(1,idr))); end
plot(mean(drbin,1),mfdr,'o-k','linewidth',2,'markerfacecolor',[0 0 0]);
ylabel('|\delta V_p| (m/s)');
set(gca,'fontsize',15,'linewidth',1.5);
axis tight;
ylim([0 max(abs(misfit_Vw))]);

subplot(3,2,6);
plot(azi_medgap,abs(misfit_Vw),'.','color',clr,'linewidth',2); hold on;
% plot(azi_Nge45,misfit_r_xy1,'o'); hold on;
% xlabel('N_{gaps} >= 20{\circ}');
y = abs(misfit_Vw);
drbin = [5:0.25:6 ; 5.25:0.25:6.25]; clear('mfdr')
for idr = 1:size(drbin,2), mfdr(idr) = median(y([azi_medgap]<= drbin(2,idr) & [azi_medgap]> drbin(1,idr))); end
% plot(mean(drbin,1),mfdr,'o-k','linewidth',2,'markerfacecolor',[0 0 0]);
ylabel('|\delta V_p| (m/s)');
set(gca,'fontsize',15,'linewidth',1.5);
axis tight;
ylim([0 max(abs(misfit_Vw))]);
% xlim([4 7]);

save2pdf('./FigureS12_Synth_azimuthalGaps.pdf',f1,100)

%% PLOT VARIANCE COMPARISON
f300 = figure(300); clf;
set(gcf,'Position',[305   102   680   596]);
clr = [0.5 0.5 0.5];

subplot(2,2,1);
plot(x_std_bs,x_std_covm,'.','color',clr,'linewidth',2); hold on;
plot([0 100],[0 100],'-k','linewidth',1.5);
text(0.1,0.9,['R^2 = ',num2str(R2x,'%.2f')],'Units','normalized','fontsize',14)
xlim([min([x_std_bs, x_std_covm]) max([x_std_bs, x_std_covm])])
ylim([min([x_std_bs, x_std_covm]) max([x_std_bs, x_std_covm])])
title('X (m)')
xlabel('2\sigma Bootstrap');
ylabel('2\sigma Covariance');
set(gca,'fontsize',15,'linewidth',1.5);
set(gca,'XTick',get(gca,'YTick'));
% axis square;

subplot(2,2,2);
plot(y_std_bs,y_std_covm,'.','color',clr,'linewidth',2); hold on;
plot([0 100],[0 100],'-k','linewidth',1.5);
text(0.1,0.9,['R^2 = ',num2str(R2y,'%.2f')],'Units','normalized','fontsize',14)
xlim([min([y_std_bs, y_std_covm]) max([y_std_bs, y_std_covm])])
ylim([min([y_std_bs, y_std_covm]) max([y_std_bs, y_std_covm])])
title('Y (m)')
xlabel('2\sigma Bootstrap');
ylabel('2\sigma Covariance');
set(gca,'fontsize',15,'linewidth',1.5);
set(gca,'XTick',get(gca,'YTick'));
% axis square;

subplot(2,2,3);
plot(z_std_bs,z_std_covm,'.','color',clr,'linewidth',2); hold on;
plot([0 100],[0 100],'-k','linewidth',1.5);
text(0.1,0.9,['R^2 = ',num2str(R2z,'%.2f')],'Units','normalized','fontsize',14)
xlim([min([z_std_bs, z_std_covm]) max([z_std_bs, z_std_covm])])
ylim([min([z_std_bs, z_std_covm]) max([z_std_bs, z_std_covm])])
title('Z (m)')
xlabel('2\sigma Bootstrap');
ylabel('2\sigma Covariance');
set(gca,'fontsize',15,'linewidth',1.5);
set(gca,'XTick',get(gca,'YTick'));
% axis square;

subplot(2,2,4);
plot(vw_std_bs,vw_std_covm,'.','color',clr,'linewidth',2); hold on;
plot([0 100],[0 100],'-k','linewidth',1.5);
text(0.1,0.9,['R^2 = ',num2str(R2vw,'%.2f')],'Units','normalized','fontsize',14)
xlim([min([vw_std_bs, vw_std_covm]) max([vw_std_bs, vw_std_covm])])
ylim([min([vw_std_bs, vw_std_covm]) max([vw_std_bs, vw_std_covm])])
title('V_p (m/s)')
xlabel('2\sigma Bootstrap');
ylabel('2\sigma Covariance');
set(gca,'fontsize',15,'linewidth',1.5);
set(gca,'XTick',get(gca,'YTick'));
% axis square;

save2pdf('./FigureS11_Synth_Bootstrap_Covm_uncertainties.pdf',f300,100)

