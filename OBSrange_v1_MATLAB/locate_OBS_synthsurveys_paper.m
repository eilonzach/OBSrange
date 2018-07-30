%% Testing various survey geometries
% Performs inversion for 1000 realizations of each survey geometry to
% evaluate which one is best able to recover the true model parameters.
% 
% Requires: 
%   Acoustic survey files in some directory 
%   Inputs below
% 
% Uses two-way travel time information and ship coordinates from an OBS 
% survey to invert for survtion location on the seafloor (Lat, Lon, Depth),
% turn-around time (TAT), survtic correction to the sound velocity through the
% water column (dvp), and the velocity of the ship in the radial
% direction of the survey circle (vr0).
%
%
%

clear; close all;

%% INPUTS - MAKE SURE THESE ARE 
% path to project
% projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA/'; % Josh
projpath = '~/Work/OBSrange/synthetics/'; 

% path to survey data from the project directory
% datapath = '/Users/russell/Lamont/PROJ_OBSrange/synth_tests_paper/synth_surveys/'; % Josh
datapath = 'synth_surveys_paper/';
% path to output directory from project directory(will be created if it does not yet exist)
outdir = './OUT_OBSrange_synthsurveys/'; 

% % path to survey data from the project directory
% %datapath = '/Users/russell/Lamont/PROJ_OBSrange/synth_tests_paper/synth_surveys/';
% datapath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/synthetics/synth_surveys_paper/noise4ms/';
% % path to output directory from project directory(will be created if it does not yet exist)
% outdir = './OUT_OBSrange_synthsurveys_noise4ms/'; 

% Put a string survtion name here to only consider that survtion. 
% Otherwise, to locate all survtions, put ''
onesurvey = '';

%% Parameters
ifsave = 1; % Save results to *.mat?
ifplot = 1; % Plot results?

par = struct([]);
par(1).vp_w = 1500; % Assumed water velocity (m/s)
par.N_bs = 500; %500; % Number of bootstrap iterations (= 1 for these tests)
par.E_thresh = 1e-5; % RMS reduction threshold for inversion

% Traveltime correction parameters
par.if_twtcorr = 1; % Apply a traveltime correction to account for ship velocity?
par.if_perfectcorr = 0;
par.npts_movingav = 1; %5; % number of points to include in moving average smoothing of ship velocity (1 = no smoothing);

% Ping QC -- Remove pings > ping_thresh ms away from neighbor
ifQC_ping = 1; % Do quality control on pings?
res_thresh = 500; % (ms) Will filter out pings with residuals > specified magnitude

% TAT - Define turnaround time to damp towards in the inversion
par.TAT_start = 0.013; % (s)
par.TAT_bounds = [0.005 0.025]; % (s) Bounds allowed for TAT

% Norm damping for each model parameter (damping towards survrting model)
% Larger values imply more damping towards the survrting model.
par.dampx = 0;
par.dampy = 0;
par.dampz = 0; %0
par.dampTAT = 2e-1; %2e-1; %2e-1
par.dampdvp = 5e-8; %5e-8

% Global norm damping for stabilization
par.epsilon = 1e-10;

%% ===================================================================== %%
%% ================ NOT ADVISED TO EDIT BELOW THIS LINE ================ %%
%% ===================================================================== %%

% prepend functions directory to MATLAB path
fullMAINpath = mfilename('fullpath');
functionspath = [fullMAINpath(1:regexp(fullMAINpath,mfilename)-1),'functions'];
addpath(functionspath);

% output directory
if par.if_twtcorr
    if par.if_perfectcorr
        modified_outdir = [outdir,'OUT_wperfectcorr/'];
        display('*APPLYING PERFECT CORRECTIONS*');
    elseif ~par.if_perfectcorr
        modified_outdir = [outdir,'OUT_wcorr/'];
    end
elseif ~par.if_twtcorr
    modified_outdir = [outdir,'OUT_nocorr/'];
end

%% Load 2-way Travel Time Data
wd = pwd;
cd(projpath);
files = dir([datapath,'/*.mat']);
surveys = extractBefore({files.name},'.mat');
Nsurveys = length(surveys);
for is = 1:Nsurveys
surv = surveys{is};
if ~isempty(onesurvey), if ~strcmp(surv,onesurvey),continue; end; end
fprintf('===========================\nWorking on %s\n\n',surv);
%% load all the data    
%     rawdatafile = sprintf('%s%s.txt',datapath,surv);
rawdatfile = dir([datapath,surv,'*']);data = [];
if exist([modified_outdir,'/mats/',surv,'_OUT.mat'],'file') == 2
    fprintf('\n%s already processed. Skipping...\n',files(is).name);
    continue
else
    load([rawdatfile.folder,'/',rawdatfile.name]);    
%     fprintf('\nWorking on: %s\n',surv);
end  
if isempty(data)
	continue;
end

% Store "global" variables for survey
survey = data(1).survey;
radius = data(1).radius;
rmsnoise = data(1).rmsnoise;
vship_kn = data(1).vship_kn;
dt_survey = data(1).dt_survey;
survx = data(1).survx*1000;
survy = data(1).survy*1000;
survt = data(1).survt;
lat_drop = data(1).drop(1);
lon_drop = data(1).drop(2);  
z_drop = data(1).drop(3)*-1000;
fprintf('%s\n',survey);

% prepare bootstrap result vectors

x_sta = nan(par.N_bs,1); y_sta = nan(par.N_bs,1); z_sta = nan(par.N_bs,1);
lon_sta = nan(par.N_bs,1); lat_sta = nan(par.N_bs,1);
TAT = nan(par.N_bs,1);
dvp = nan(par.N_bs,1);
V_w = nan(par.N_bs,1);
E_rms = nan(par.N_bs,1);
v_eff = nan(par.N_bs,1);
dx_drift = nan(par.N_bs,1);
dy_drift = nan(par.N_bs,1);
dz_sta = nan(par.N_bs,1);
drift = nan(par.N_bs,1);
azi  = nan(par.N_bs,1);
dist = cell(par.N_bs,1);
azi_locs = cell(par.N_bs,1);
models = cell(par.N_bs,1);

R_mats = zeros(5,5,length(data));
Cm_mats = zeros(5,5,length(data));

for ii = 1:length(data)
    if mod(ii,100)==0, fprintf('Iteration %u\n',ii); end
    twt = data(ii).tot_dt;
    [Is0nan, ~] = find(~isnan(twt));
    if length(Is0nan) < 2
        fprintf('\nOnly %d pings for iteration %d. Skipping...\n',length(Is0nan),ii);
        continue
    end
    twt = twt(Is0nan);
    lats_ship = data(ii).survlats(Is0nan);
    lons_ship = data(ii).survlons(Is0nan);
    t_ship = data(ii).survts(Is0nan)*24*60*60;
    corr_dt = data(ii).corr_dt(Is0nan);
    v_surv_true = data(ii).v_surv_true(Is0nan,:);
    
    % Set origin of coordinate system to be lat/lon of drop point
    olon = lon_drop;
    olat = lat_drop;

    Nobs = length(twt);
    z_ship = zeros(Nobs,1); % ship is always at surface

    % Convert Lon/Lat to x/y
    [ x_ship, y_ship ] = lonlat2xy_nomap( olon, olat, lons_ship, lats_ship );
    [ x_drop, y_drop ] = lonlat2xy_nomap( olon, olat, lon_drop, lat_drop );

    % Calculate velocity of ship
    v_ship = pt_veloc( x_ship, y_ship, z_ship, t_ship );
    v_ship = [moving_average(v_ship(1,:),par.npts_movingav)'; moving_average(v_ship(2,:),par.npts_movingav)'; moving_average(v_ship(3,:),par.npts_movingav)'];

    %% Set up initial model
    m0_strt(1,1) = x_drop; %x0;
    m0_strt(2,1) = y_drop; %y0;
    m0_strt(3,1) = z_drop; %z0;
    m0_strt(4,1) = par.TAT_start; %TAT;
    m0_strt(5,1) = 0; %dvp;

    %% Do Bootstrap Inversion
    Nobs = length(twt);
    M = length(m0_strt);
    % damping matrix
    H = eye(M, M) .* diag([par.dampx, par.dampy, par.dampz, par.dampTAT, par.dampdvp]);

    [xmat_ship_bs, ymat_ship_bs, zmat_ship_bs, vmat_ship_bs, twtmat_bs, indxs] = bootstrap(x_ship, y_ship, z_ship, v_ship, twt, par.N_bs-1);
    dtwt_mat = zeros(size(indxs));
    twtcorr_mat = zeros(size(indxs));
    dtwtcorr_mat = zeros(size(indxs));
    vr_mat = zeros(size(indxs));
    for ibs = 1:par.N_bs
        x_ship_bs = xmat_ship_bs(:,ibs);
        y_ship_bs = ymat_ship_bs(:,ibs);
        z_ship_bs = zmat_ship_bs(:,ibs);
        v_ship_bs = [vmat_ship_bs(1).amatbs(:,ibs)'; vmat_ship_bs(2).amatbs(:,ibs)'; vmat_ship_bs(3).amatbs(:,ibs)'];
        twt_bs = twtmat_bs(:,ibs);
        
        if ~par.if_perfectcorr
            [ m_final,models,v,N,R,Cm ] = ...
                inv_newtons( par,m0_strt,twt_bs,...
                            x_ship_bs,y_ship_bs,z_ship_bs,...
                            v_ship_bs,H);
        elseif par.if_perfectcorr
            [ m_final,models,v ] = ...
                inv_newtons_perfectcorr( par,m0_strt,twt_bs,...
                            x_ship_bs,y_ship_bs,z_ship_bs,...
                            v_ship_bs,H,corr_dt);
        end

        x_sta(ibs) = m_final(1);
        y_sta(ibs) = m_final(2);
        z_sta(ibs) = m_final(3);
        TAT(ibs) = m_final(4);
        dvp(ibs) = m_final(5);
        V_w(ibs) = par.vp_w + dvp(ibs);
        E_rms(ibs) = models(end).E;
        v_eff(ibs) = v;
        dtwt_mat(:,ibs) = models(end).dtwt;
        twtcorr_mat(:,ibs) = models(end).twt_corr;
        dtwtcorr_mat(:,ibs) = models(end).dtwtcorr;
        vr_mat(:,ibs) = models(end).vr;
        [lon_sta(ibs), lat_sta(ibs)] = xy2lonlat_nomap(olon, olat, x_sta(ibs), y_sta(ibs));

        % Calculate OBS drift distance and azimuth
        dx_drift(ibs) = m_final(1) - x_drop;
        dy_drift(ibs) = m_final(2) - y_drop;
        dz_sta(ibs) = m_final(3) - z_drop;
        drift(ibs) = sqrt(dx_drift(ibs).^2 + dy_drift(ibs).^2);
        azii = -atan2d(dy_drift(ibs),dx_drift(ibs))+90;
        azii(azii<0) = azii(azii<0) + 360;
        azi(ibs) = azii;

        % Calculate ship distance azimuth from OBS
        dx_ship = x_ship_bs-m_final(1);
        dy_ship = y_ship_bs-m_final(2);
        dist{ibs} = sqrt(dx_ship.^2 + dy_ship.^2);
        azi_locss = -atan2d(dy_ship , dx_ship) + 90;
        azi_locss(azi_locss<0) = azi_locss(azi_locss<0) + 360;
        azi_locs{ibs} = azi_locss;
    end % loop on bootstrap (if any)
    dtwt_bs = mean(unscramble_randmat(dtwt_mat,indxs),2);
    twtcorr_bs = mean(unscramble_randmat(twtcorr_mat,indxs),2);
    dtwtcorr_bs = mean(unscramble_randmat(dtwtcorr_mat,indxs),2);
    vr_bs = mean(unscramble_randmat(vr_mat,indxs),2);

    range = sqrt( (mean(x_sta)-x_ship).^2 + (mean(y_sta)-y_ship).^2 + (mean(z_sta)-z_ship).^2 );
    
    % Store output
    data(ii).lat_sta = mean(lat_sta);
    data(ii).lon_sta = mean(lon_sta);
    data(ii).x_sta = mean(x_sta);           if par.N_bs>1, data(ii).x_sta_std = std(x_sta); end
    data(ii).y_sta = mean(y_sta);           if par.N_bs>1, data(ii).y_sta_std = std(y_sta); end
    data(ii).z_sta = mean(z_sta);           if par.N_bs>1, data(ii).z_sta_std = std(z_sta); end
    data(ii).TATsta = mean(TAT);            if par.N_bs>1, data(ii).TAT_std = std(TAT); end
    data(ii).V_w = mean(V_w);               if par.N_bs>1, data(ii).V_w_std = std(V_w); end
    data(ii).drift = mean(drift);           if par.N_bs>1, data(ii).drift_std = std(drift); end
    data(ii).drift_true = sqrt(data(ii).obs_loc_xyz(1).^2 + data(ii).obs_loc_xyz(2).^2)*1000;
    misfit_latsta(ii) = data(ii).lat_sta - data(ii).obs_loc_laloz(1);
    misfit_lonsta(ii) = data(ii).lon_sta - data(ii).obs_loc_laloz(2);
    misfit_xsta(ii) = data(ii).x_sta - data(ii).obs_loc_xyz(1)*1000;
    misfit_ysta(ii) = data(ii).y_sta - data(ii).obs_loc_xyz(2)*1000;
    misfit_zsta(ii) = data(ii).z_sta - (-data(ii).obs_loc_xyz(3)*1000);
    misfit_r_xy(ii) = sqrt( misfit_xsta(ii).^2 + misfit_ysta(ii).^2 );
    misfit_r_xyz(ii) = sqrt( misfit_xsta(ii).^2 + misfit_ysta(ii).^2 + misfit_zsta(ii).^2 );
    misfit_TAT(ii) = data(ii).TATsta - data(ii).TAT;
    misfit_Vw(ii) = data(ii).V_w - data(ii).Vp_water*1000;
    
    data(ii).Is0nan = Is0nan(:);
    data(ii).twtcorr = twtcorr_bs(:);
    data(ii).dtwt = dtwt_bs(:)';
    data(ii).dtwtcorr = dtwtcorr_bs(:);
    data(ii).misfit_dtwtcorr = data(ii).dtwtcorr(:)' - corr_dt(:)'; % correction time misfit
    data(ii).v_ship = v_ship;
    data(ii).misfit_v_ship = v_ship(1:2,:) - v_surv_true';
    

    R_mats(:,:,ii) = R;
    Cm_mats(:,:,ii) = Cm;
    data(ii).R_mat = R;
    data(ii).Cm_mat = Cm;
end % loop on synth stations
R_mat = mean(R_mats,3);
Cm_mat = mean(Cm_mats,3);

% Store output
data(1).misfit_latsta = misfit_latsta(:);
data(1).misfit_lonsta = misfit_lonsta(:);
data(1).misfit_xsta = misfit_xsta(:);
data(1).misfit_ysta = misfit_ysta(:);
data(1).misfit_zsta = misfit_zsta(:);
data(1).misfit_r_xy = misfit_r_xy(:);
data(1).misfit_r_xyz = misfit_r_xyz(:);
data(1).misfit_TAT = misfit_TAT(:);
data(1).misfit_Vw = misfit_Vw(:);
data(1).E_rms = E_rms(:);
data(1).misfit_v_ship_all = [data(:).misfit_v_ship];
data(1).misfit_dtwtcorr_all = [data(:).misfit_dtwtcorr];
data(1).dtwt_all = [data(:).dtwt];

%% Print out some important parameters
fprintf('\nStatistics for survey %s\n\n',surv);
fprintf('Mean x-misfit   = %.3f m\n',mean(data(1).misfit_xsta));
fprintf('Mean y-misfit   = %.3f m\n',mean(data(1).misfit_ysta));
fprintf('Mean z-misfit   = %.3f m\n',mean(data(1).misfit_zsta));
fprintf('Mean Vw-misfit  = %.3f m\n',mean(data(1).misfit_Vw));
fprintf('Mean tat-misfit = %.3f m\n\n',mean(data(1).misfit_TAT));
fprintf('Abs 2D error = %.3f ± %.3f m\n',mean(data(1).misfit_r_xy),std(data(1).misfit_r_xy));
fprintf('Abs 3D error = %.3f ± %.3f m\n',mean(data(1).misfit_r_xyz),std(data(1).misfit_r_xyz));
sor2De = sort(data(1).misfit_r_xy);
fprintf('95%% percentile 2D error = %.2f m\n',sor2De(round(0.95*length(sor2De))));
fprintf('Abs depth error = %.3f m\n',mean(abs(data(1).misfit_zsta)));
if par.N_bs>1
fprintf('Fraction of x-misfit outside bootstrap 2sig = %.2f %%\n',  100*sum(abs(data(1).misfit_xsta)>2*[data.x_sta_std]')/length(data(1).misfit_xsta));
fprintf('Fraction of y-misfit outside bootstrap 2sig = %.2f %%\n',  100*sum(abs(data(1).misfit_ysta)>2*[data.y_sta_std]')/length(data(1).misfit_ysta));
fprintf('Fraction of z-misfit outside bootstrap 2sig = %.2f %%\n',  100*sum(abs(data(1).misfit_zsta)>2*[data.z_sta_std]')/length(data(1).misfit_zsta));
fprintf('Fraction of Vw-misfit outside bootstrap 2sig = %.2f %%\n', 100*sum(abs(data(1).misfit_Vw)  >2*[data.V_w_std]')  /length(data(1).misfit_Vw));
fprintf('Fraction of tat-misfit outside bootstrap 2sig = %.2f %%\n',100*sum(abs(data(1).misfit_TAT) >2*[data.TAT_std]')  /length(data(1).misfit_TAT));
r_sta_std = (abs([data.x_sta]).*[data.x_sta_std] + abs([data.y_sta]).*[data.y_sta_std])./[data.drift];
fprintf('Fraction of 2d-misfit outside bootstrap 2sig = %.2f %%\n', 100*sum(abs(data(1).misfit_r_xy)>2*r_sta_std')/length(data(1).misfit_r_xy));
end


%% HISTOGRAMS OF MODEL PARAMETERS
if ifplot
PLOT_histograms_all_synthsurveys
end

%% PLOTTING
if ifplot
	%% Geographic Plotting
    
    % Plot in map view and depth
    f1 = figure(102); clf;
    set(gcf,'position',[53 292 1157 413]);
    clr = lines(5);
    
    [minE_rms, IminE_rms] = min(E_rms);
    [maxE_rms, ImaxE_rms] = max(E_rms);
    subplot(1,2,1);
    plot(data(IminE_rms).survlons(data(IminE_rms).Is0nan),data(IminE_rms).survlats(data(IminE_rms).Is0nan),'ok','markerfacecolor',clr(2,:),'markersize',13,'linewidth',1); hold on;
    plot(data(1).drop(2),data(1).drop(1),'sk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1); hold on;
    plot(data(IminE_rms).lon_sta,data(IminE_rms).lat_sta,'pk','markerfacecolor',[1 1 0],'markersize',25,'linewidth',1)
    plot(data(IminE_rms).obs_loc_laloz(2),data(IminE_rms).obs_loc_laloz(1),'p','color',clr(5,:),'markersize',25,'linewidth',2)
    grid on; box on; set(gca,'fontsize',16,'linewidth',2);
    title(['Min RMS: ',num2str(E_rms(IminE_rms)*1000,'%.2f'),' ms'],'fontsize',18);
    xlabel('Longitude','fontsize',18);
    ylabel('Latitude','fontsize',18);
    axis equal;
    
    subplot(1,2,2);
    plot(data(ImaxE_rms).survlons(data(ImaxE_rms).Is0nan),data(ImaxE_rms).survlats(data(ImaxE_rms).Is0nan),'ok','markerfacecolor',clr(2,:),'markersize',13,'linewidth',1); hold on;
    plot(data(1).drop(2),data(1).drop(1),'sk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1); hold on;
    plot(data(ImaxE_rms).lon_sta,data(ImaxE_rms).lat_sta,'pk','markerfacecolor',[1 1 0],'markersize',25,'linewidth',1)
    plot(data(ImaxE_rms).obs_loc_laloz(2),data(ImaxE_rms).obs_loc_laloz(1),'p','color',clr(5,:),'markersize',25,'linewidth',2)
    grid on; box on; set(gca,'fontsize',16,'linewidth',2);
    title(['Max RMS: ',num2str(E_rms(ImaxE_rms)*1000,'%.2f'),' ms'],'fontsize',18);
    xlabel('Longitude','fontsize',18);
    ylabel('Latitude','fontsize',18);
    axis equal;
    
    %% Model resolution and Covariance
    PLOT_resolution_covariance
    %%

    % drift vs. misfit
    figure(107), clf, hold on
    plot([data.drift_true],data(1).misfit_r_xy,'.')
    drbin = [0:50:350;50:50:400]; clear('mfdr')
    for idr = 1:size(drbin,2), mfdr(idr) = median(data(1).misfit_r_xy([data.drift_true]<= drbin(2,idr) & [data.drift_true]> drbin(1,idr))); end
    plot(mean(drbin,1),mfdr,'o-r')
    xlabel('Drift distance (m)'); ylabel('Horizontal misfit (m)')
    % depth vs. tat
    figure(108)
    plot(data(1).misfit_zsta,data(1).misfit_Vw,'.')
    xlabel('Z misfit (m)'); ylabel('Vw misfit (m)')
end


%% Save output

% % Save textfile
% fid = fopen([modified_outdir,'/',data.surv,'_location.txt'],'w');
% fprintf(fid,'Bootstrap inversion results (2sigma uncertainty)');
% fprintf(fid,'\nsurvtion: %s',data.surv);
% fprintf(fid,'\nLat:   %.5f deg (%f) \nLon:   %.5f deg (%f) \nX:     %f m (%f) \nY:    %f m (%f) \nDepth: %f m (%f) \nTAT:   %f ms (%f) \nWater Vel.: %f m/s (%f)',mean(lat_sta),std(lat_sta)*2,mean(lon_sta),std(lon_sta)*2,mean(x_sta),std(x_sta)*2,mean(y_sta),std(y_sta)*2,mean(z_sta),std(z_sta)*2,mean(TAT)*1000,std(TAT)*1000*2,mean(V_w),std(V_w)*2);
% fprintf(fid,'\nDrift Lon: %f m (%f) \nDrift Lat: %f m (%f) \nDrift:    %f m (%f) \nDrift Azi: %f deg (%f)\ndz: %f m (%f)\n',mean(dx_drift),std(dx_drift)*2,mean(dy_drift),std(dy_drift)*2,mean(drift),std(drift)*2,mean(azi),std(azi)*2,mean(dz_sta),std(dz_sta)*2);
% fprintf(fid,'\nRMS:  %f ms (%f)\n',mean(E_rms)*1000,std(E_rms)*2*1000);
% fprintf(fid,'\nBad pings Removed: %d',N_badpings);
% fprintf(fid,'\n===================================================\n');
% fprintf(fid,'%10s %10s %15s %15s %15s %15s \n','Lat','Lon','Range (m)','Residual (s)','Doppler Vel. (m/s)','TWT corr. (ms)');
% for ii = 1:Nobs
% 	fprintf(fid,'%3d: %10f %10f %10f %10f %10f %10f\n',ii,lats_ship(ii),lons_ship(ii),range(ii),dtwt_bs(ii)*1000,vr_bs(ii),dtwtcorr_bs(ii)*1000);
% end
% fclose(fid);

if ifsave && ifplot
% Save plots
if ~exist([modified_outdir,'/plots/'])
	mkdir([modified_outdir,'/plots/']);
end
save2pdf([modified_outdir,'/plots/',surv,'_1_ModelHists.pdf'],100,500)
save2pdf([modified_outdir,'/plots/',surv,'_2_MisfitHists.pdf'],101,500)
save2pdf([modified_outdir,'/plots/',surv,'_3_Survey.pdf'],102,500)
save2pdf([modified_outdir,'/plots/',surv,'_4_Resolution.pdf'],103,500)
end

if ifsave
	if ~exist([modified_outdir,'/mats'])
		mkdir([modified_outdir,'/mats']);
	end
	save([modified_outdir,'/mats/',surv,'_OUT.mat'],'data');
end


end

% message if no success
if is==Nsurveys && ~exist('rawdatfile','var')
    fprintf('*********************\nsurvtion %s does not seem to exist in that folder\n',onesurvey);
end
cd (wd)
