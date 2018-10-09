% Save synthetic bootstrap results in smaller mat files for making figures
% for paper.
%

clear; close all;

%% INPUTS - MAKE SURE THESE ARE 
% path to project
% projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_SynthBoot_surveys_strtmodel/'; 
projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_SynthBoot_surveys_noTAT/'; 
% path to survey data from the project directory
datapath = '/Users/russell/Lamont/PROJ_OBSrange/synth_tests_paper/synth_surveys/'; 
% path to output directory from project directory (actually input for this script)
outdir = [projpath,'/OUT_OBSrange_synthsurveys/OUT_wcorr/']; %

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
    data_summary.x_ship = x_ship;
    data_summary.y_ship = y_ship;
    data_summary.v_surv_true = data(1).v_surv_true;
    data_summary.vmag = sqrt(sum(data(1).v_surv_true.^2,2));
    data_summary.survx = data(1).survx;
    data_summary.survy = data(1).survy;
    data_summary.survey = data(1).survey;
    data_summary.radius = data(1).radius;
        
    
    data_summary.misfit_xsta = rms(data(1).misfit_xsta);
    data_summary.misfit_xsta_std = std(data(1).misfit_xsta);
    data_summary.misfit_xsta_mean = mean(data(1).misfit_xsta);
    data_summary.misfit_xsta_absmean = mean(abs(data(1).misfit_xsta));
    data_summary.misfit_xsta_absmed = median(abs(data(1).misfit_xsta));
    data_summary.misfit_xsta_med = median(data(1).misfit_xsta);
    data_summary.misfit_ysta = rms(data(1).misfit_ysta);
    data_summary.misfit_ysta_std = std(data(1).misfit_ysta);
    data_summary.misfit_ysta_mean = mean(data(1).misfit_ysta);
    data_summary.misfit_ysta_med = median(data(1).misfit_ysta);
    data_summary.misfit_ysta_absmean = mean(abs(data(1).misfit_ysta));
    data_summary.misfit_ysta_absmed = median(abs(data(1).misfit_ysta));
    data_summary.misfit_zsta = rms(data(1).misfit_zsta);
    data_summary.misfit_zsta_std = std(data(1).misfit_zsta);
    data_summary.misfit_zsta_mean = mean(data(1).misfit_zsta);
    data_summary.misfit_zsta_absmean = mean(abs(data(1).misfit_zsta));
    data_summary.misfit_zsta_absmed = median(abs(data(1).misfit_zsta));
    data_summary.misfit_zsta_med = median(data(1).misfit_zsta);
    data_summary.misfit_r_xy = rms(data(1).misfit_r_xy);
    data_summary.misfit_r_xy_std = std(data(1).misfit_r_xy);
    data_summary.misfit_r_xy_mean = mean(data(1).misfit_r_xy);
    data_summary.misfit_r_xy_med = median(data(1).misfit_r_xy);
    data_summary.misfit_r_xyz = rms(data(1).misfit_r_xyz);
    data_summary.misfit_r_xyz_std = std(data(1).misfit_r_xyz);
    data_summary.misfit_r_xyz_mean = mean(data(1).misfit_r_xyz);
    data_summary.misfit_r_xyz_med = median(data(1).misfit_r_xyz);
    data_summary.misfit_TAT = rms(data(1).misfit_TAT);
    data_summary.misfit_TAT_std = std(data(1).misfit_TAT);
    data_summary.misfit_TAT_mean = mean(data(1).misfit_TAT);
    data_summary.misfit_TAT_absmean = mean(abs(data(1).misfit_TAT));
    data_summary.misfit_TAT_absmed = median(abs(data(1).misfit_TAT));
    data_summary.misfit_TAT_med = median(data(1).misfit_TAT);
    data_summary.misfit_Vw = rms(data(1).misfit_Vw);
    data_summary.misfit_Vw_std = std(data(1).misfit_Vw);
    data_summary.misfit_Vw_mean = mean(data(1).misfit_Vw);
    data_summary.misfit_Vw_absmean = mean(abs(data(1).misfit_Vw));
    data_summary.misfit_Vw_absmed = median(abs(data(1).misfit_Vw));
    data_summary.misfit_Vw_med = median(data(1).misfit_Vw);
    data_summary.E_rms = mean(data(1).E_rms);
    data_summary.E_rms_std = std(data(1).E_rms);
    data_summary.misfit_v_ship_all = mean(data(1).misfit_v_ship_all,2);
    data_summary.misfit_v_ship_all_std = std(data(1).misfit_v_ship_all,0,2);
    data_summary.misfit_dtwtcorr_all = rms(data(1).misfit_dtwtcorr_all);
    data_summary.misfit_dtwtcorr_all_std = std(data(1).misfit_dtwtcorr_all);
    data_summary.dtwt_all = rms(data(1).dtwt_all);
    data_summary.dtwt_all_std = std(data(1).dtwt_all);
    
    for ii = 1:length(data)
        R_mats(:,:,ii) = data(ii).R_mat;
        Cm_mats(:,:,ii) = data(ii).Cm_mat;
    end % loop on synth stations
    data_summary.R_mat = mean(R_mats,3);
    data_summary.Cm_mat = mean(Cm_mats,3);
    
    %% Save summary mat files
    if ~exist([outdir,'/mats_SynthBoot_summary/'])
        mkdir([outdir,'/mats_SynthBoot_summary/']);
    end   
    save([outdir,'/mats_SynthBoot_summary/',files(ifil).name],'data_summary')
    
end

