%% Script to print out values used in the paper
clear 

%% Section 3.1

% Success of one-station theoretical test
one_sta_syn = 'syn10';
truedat = load(['../figs/figdata/trudata_',one_sta_syn,'.mat']);
load(['../figs/figdata/',one_sta_syn,'_data.mat']);
% misfit in location
dr = norm(truedat.obs_location_xyz(1:2) - datamat.loc_xyz(1:2)/1e3);
fprintf('Station %s is mis-located by %.2f m\n',one_sta_syn,1e3*dr)
fprintf('Total true drift %.2f m\n',1e3*norm(truedat.obs_location_xyz(1:2)))
fprintf('True water depth %.0f m\n',1e3*truedat.obs_location_xyz(3))
fprintf('Survey diameter %.0f m\n',2*truedat.radius*truedat.nm2km*1e3)

% % Statistics for survey SynthBoot_PACMAN_rad1.00
% % 1e4 iterations of station locations, each with 1e3 bootstrap iterations
% % 
% % Mean x-misfit   = 0.04 m
% % Mean y-misfit   = 0.15 m
% % Mean z-misfit   = -0.60 m
% % Mean Vw-misfit  = 0.17 m
% % Mean tat-misfit = 0.00 m
% % 
% % Abs 2D error = 2.30 ± 1.22 m
% % Abs 3D error = 8.29 ± 5.43 m
% % 95% percentile 2D error = 4.58 m
% % Fraction of x-misfit outside bootstrap 2sig = 5.18 %
% % Fraction of y-misfit outside bootstrap 2sig = 5.42 %
% % Fraction of z-misfit outside bootstrap 2sig = 6.54 %
% % Fraction of Vw-misfit outside bootstrap 2sig = 6.18 %
% % Fraction of tat-misfit outside bootstrap 2sig = 85.62 %
% % Fraction of 2d-misfit outside bootstrap 2sig = 5.67 %

%% Section 3.2
fprintf('\n ============ yORCA locations ============\n\n')

% results of all station inversions for young Pac ORCA
load('../data/yORCA_locations.mat');
Nstas = length(allstas);
drop_lolaz = reshape([allstas.drop_lonlatz],3,Nstas)';
loc_lolaz = reshape([allstas.loc_lolaz],3,Nstas)';
mean_drift_az = reshape([allstas.mean_drift_az],2,Nstas)';
% some stats
Erms = zeros(Nstas,1); 
for is = 1:Nstas
    Vp(is) = mean(allstas(is).V_w_bs);
    Erms(is) = mean(allstas(is).E_rms);
	r_sta_std(is) = (abs(mean(allstas(is).x_sta_bs))*std(allstas(is).x_sta_bs) + abs(mean(allstas(is).y_sta_bs))*std(allstas(is).y_sta_bs))./mean_drift_az(is,1);
    F_sta_std(is) = allstas(is).Ftest_res.uncertainties.xy68;
end
% 
fprintf('Mean drift = %.1f\n',mean(mean_drift_az(:,1)));
fprintf('Mean RMS error = %.2f ms\n',1e3*mean(Erms));
fprintf('Mean 2d location std (boots) = %.2f m\n',mean(r_sta_std));
fprintf('Mean 2d location std (ftest) = %.2f m\n',mean(F_sta_std));
