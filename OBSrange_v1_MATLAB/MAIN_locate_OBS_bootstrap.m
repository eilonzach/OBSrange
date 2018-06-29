% MAIN script to run location inversion for an OBS deployment
% 
% Requires: 
%   Acoustic survey files in some directory 
%   Inputs below
% 
% Uses two-way travel time information and ship coordinates from an OBS 
% survey to invert for station location on the seafloor (Lat, Lon, Depth),
% turn-around time (TAT), static correction to the sound velocity through the
% water column (dvp), and the velocity of the ship in the radial
% direction of the survey circle (vr0).
%
% Josh Russell & Zach Eilon 4/16/18

clear; close all;

%% INPUTS - MAKE SURE THESE ARE 
% path to project
projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA/'; 
% path to survey data from the project directory
datapath = './'; 
% path to output directory from project directory(will be created if it does not yet exist)
outdir = './OUT_OBSrange/'; 
% Put a string station name here to only consider that station. 
% Otherwise, to locate all stations, put ''
onesta = 'CC01';

%% Parameters
ifsave = 1; % Save results to *.mat?
ifplot = 1; % Plot results?

par = struct([]);
par(1).vp_w = 1500; % Assumed water velocity (m/s)
par.N_bs = 500; % Number of bootstrap iterations
par.E_thresh = 1e-5; % RMS reduction threshold for inversion

% Traveltime correction parameters
par.if_twtcorr = 0; % Apply a traveltime correction to account for ship velocity?
par.npts_movingav = 1; %5; % number of points to include in moving average smoothing of ship velocity (1 = no smoothing);

% Ping QC -- Remove pings > ping_thresh ms away from neighbor
ifQC_ping = 1; % Do quality control on pings?
res_thresh = 500; % (ms) Will filter out pings with residuals > specified magnitude

% TAT - Define turnaround time to damp towards in the inversion
par.TAT_start = 0.013; % (s)

% Norm damping for each model parameter (damping towards starting model)
% Larger values imply more damping towards the starting model.
par.dampx = 0;
par.dampy = 0;
par.dampz = 0; %0
par.dampTAT = 2e-1; %2e-1
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

%% Load 2-way Travel Time Data
wd = pwd;
cd(projpath);
files = dir([datapath,'/*.txt']);
stas = unique(strtok({files.name},{'_','.txt'}));
Nstas = length(stas);
for is = 1:Nstas
sta = stas{is};
if ~isempty(onesta), if ~strcmp(sta,onesta),continue; end; end
fprintf('===========================\nLocating OBS %s\n\n',sta);
%% load all the data    
%     rawdatafile = sprintf('%s%s.txt',datapath,sta);
rawdatfile = dir([datapath,sta,'*']);data = [];
for irf = 1:length(rawdatfile)
	if any(regexp(rawdatfile(irf).name,'orrect')) 
		continue;
	end
data = load_pings([datapath,rawdatfile(irf).name]);        
end
if isempty(data) || isempty(data.lon_drop) || isempty(data.lat_drop)
	continue;
end
if ifQC_ping
	[dataQC, data_bad] = pingQC(data, par.vp_w, res_thresh);
	data = dataQC;
    N_badpings = length(data_bad.twt);
    fprintf('\nNumber of bad pings removed: %d\n',N_badpings);
end
lat_drop = data.lat_drop;
lon_drop = data.lon_drop;
z_drop = data.z_drop;
lats_ship = data.lats;
lons_ship = data.lons;
t_ship = data.t_ship;
twt = data.twt;
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
% prepare bootstrap result vectors
dtwt_mat = zeros(size(indxs));
twtcorr_mat = zeros(size(indxs));
dtwtcorr_mat = zeros(size(indxs));
vr_mat = zeros(size(indxs));
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
for ibs = 1:par.N_bs
    x_ship_bs = xmat_ship_bs(:,ibs);
    y_ship_bs = ymat_ship_bs(:,ibs);
    z_ship_bs = zmat_ship_bs(:,ibs);
    v_ship_bs = [vmat_ship_bs(1).amatbs(:,ibs)'; vmat_ship_bs(2).amatbs(:,ibs)'; vmat_ship_bs(3).amatbs(:,ibs)'];
    twt_bs = twtmat_bs(:,ibs);
    
    [ m_final,models,v ] = ...
        inv_newtons( par,m0_strt,twt_bs,...
                    x_ship_bs,y_ship_bs,z_ship_bs,...
                    v_ship_bs,H);

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
end
dtwt_bs = mean(unscramble_randmat(dtwt_mat,indxs),2);
twtcorr_bs = mean(unscramble_randmat(twtcorr_mat,indxs),2);
dtwtcorr_bs = mean(unscramble_randmat(dtwtcorr_mat,indxs),2);
vr_bs = mean(unscramble_randmat(vr_mat,indxs),2);

range = sqrt( (mean(x_sta)-x_ship).^2 + (mean(y_sta)-y_ship).^2 + (mean(z_sta)-z_ship).^2 );

%% HISTOGRAMS OF MODEL PARAMETERS
if ifplot
PLOT_histograms_all
end

%% F-test for uncertainty using grid search
% Set grid size
ngridpts = 40;
D = max([ std(x_sta) std(y_sta) std(z_sta) ]*4);
Dx = D;
Dy = D; 
Dz = D;
dx = 2*Dx/ngridpts;
dy = 2*Dy/ngridpts;
dz = 2*Dz/ngridpts;
x_grid = [mean(x_sta)-Dx:dx:mean(x_sta)+Dx];
y_grid = [mean(y_sta)-Dy:dy:mean(y_sta)+Dy];
z_grid = [mean(z_sta)-Dz:dz:mean(z_sta)+Dz];
[lon_grid, lat_grid] = xy2lonlat_nomap(olon, olat, x_grid, y_grid);
Nx = length(x_grid);
Ny = length(y_grid);
Nz = length(z_grid);

% Bootstrap residual
twt_pre_bs = calcTWT(mean(x_sta), mean(y_sta), mean(z_sta), mean(dvp), mean(TAT), x_ship, y_ship, z_ship, par.vp_w);
resid_bs = twtcorr_bs-twt_pre_bs;

% Determine the eigenvectors for z_sta, V_w, and TAT
X = [z_sta, V_w, TAT];
[V, ~] = eig(X'*X);
eigvec1 = V(:,1); % Closest to TAT axis
eigvec2 = V(:,2); % Closest to V_w axis
eigvec3 = V(:,3); % Closest to z_sta axis
eig3_z = eigvec3(1);
eig3_vw = eigvec3(2);
eig3_TAT = eigvec3(3);

% Grid search
P = zeros(Nx,Ny,Nz);
E_gs = zeros(Nx,Ny,Nz);
[Xgrd,Ygrd,Zgrd] = meshgrid(x_grid,y_grid,z_grid);
[LONgrd,LATgrd,~] = meshgrid(lon_grid,lat_grid,z_grid);
Npara = length(m_final);
for ix = 1:Nx
	for iy = 1:Ny
		for iz = 1:Nz
			% Apply scaling to vp_w and TAT to account for tradeoffs with Z
			dz = Zgrd(ix,iy,iz) - mean(z_sta);
			dvw = (eig3_vw/eig3_z)*dz; % perturbation to water velocity to account for dz
			dTAT = (eig3_TAT/eig3_z)*dz; % perturbation to TAT to account for dz

			% Grid search residual;
			twt_pre_gs = calcTWT(Xgrd(ix,iy,iz), Ygrd(ix,iy,iz), Zgrd(ix,iy,iz), mean(dvp)+dvw, mean(TAT)+dTAT, x_ship, y_ship, z_ship, par.vp_w);
			resid_gs = twtcorr_bs-twt_pre_gs;

			% Calculate P statistic
			P(ix,iy,iz) = ftest_dof( resid_gs,mean(v_eff),resid_bs,mean(v_eff) );

			E_gs(ix,iy,iz) = sqrt(resid_gs'*resid_gs/length(resid_gs));
		end
	end
end

[Pz_max, Iz_max] = max(max(max(P)));
[Py_max, Iy_max] = max(max(P(:,:,Iz_max)));
[Px_max, Ix_max] = max(P(:,Iy_max,Iz_max));

Ftest_res = struct('x_grid',x_grid,'y_grid',y_grid,'z_grid',z_grid,...
                  'Pstat',P,'Erms',E_gs);

fprintf('\nStation: %s',data.sta);
fprintf('\nlat:   %.5f deg (%f) \nlon:   %.5f deg (%f) \nx:     %f m (%f) \ny:    %f m (%f) \ndepth: %f m (%f) \nTAT:   %f ms (%f) \nv_H20: %f m/s (%f)',mean(lat_sta),std(lat_sta)*2,mean(lon_sta),std(lon_sta)*2,mean(x_sta),std(x_sta)*2,mean(y_sta),std(y_sta)*2,mean(z_sta),std(z_sta)*2,mean(TAT)*1000,std(TAT)*1000*2,mean(V_w),std(V_w)*2);
fprintf('\nDrift Lon: %f m (%f) \nDrift Lat: %f m (%f) \nDrift:    %f m (%f) \nDrift Azi: %f deg (%f)\ndz: %f m (%f)\n',mean(dx_drift),std(dx_drift)*2,mean(dy_drift),std(dx_drift)*2,mean(drift),std(drift)*2,mean(azi),std(azi)*2,mean(dz_sta),std(dz_sta)*2);
fprintf('\nRMS:  %f ms (%f)\n',mean(E_rms)*1000,std(E_rms)*2*1000);

%% PLOTTING
if ifplot
	%% F-plot plots
	PLOT_Ftest_all
	%% Plot Misfit
	PLOT_misfit
	%% Geographic PLOTTING
	PLOT_survey
	PLOT_twt_corr
end


%% Save output
% output directory
if par.if_twtcorr
    modified_outdir = [outdir,'OUT_wcorr/'];
elseif ~par.if_twtcorr
    modified_outdir = [outdir,'OUT_nocorr/'];
end  
if ~exist(modified_outdir)
	mkdir(modified_outdir);
end

% Save textfile
fid = fopen([modified_outdir,'/',data.sta,'_location.txt'],'w');
fprintf(fid,'Bootstrap inversion results (2sigma uncertainty)');
fprintf(fid,'\nStation: %s',data.sta);
fprintf(fid,'\nLat:   %.5f deg (%f) \nLon:   %.5f deg (%f) \nX:     %f m (%f) \nY:    %f m (%f) \nDepth: %f m (%f) \nTAT:   %f ms (%f) \nWater Vel.: %f m/s (%f)',mean(lat_sta),std(lat_sta)*2,mean(lon_sta),std(lon_sta)*2,mean(x_sta),std(x_sta)*2,mean(y_sta),std(y_sta)*2,mean(z_sta),std(z_sta)*2,mean(TAT)*1000,std(TAT)*1000*2,mean(V_w),std(V_w)*2);
fprintf(fid,'\nDrift Lon: %f m (%f) \nDrift Lat: %f m (%f) \nDrift:    %f m (%f) \nDrift Azi: %f deg (%f)\ndz: %f m (%f)\n',mean(dx_drift),std(dx_drift)*2,mean(dy_drift),std(dy_drift)*2,mean(drift),std(drift)*2,mean(azi),std(azi)*2,mean(dz_sta),std(dz_sta)*2);
fprintf(fid,'\nRMS:  %f ms (%f)\n',mean(E_rms)*1000,std(E_rms)*2*1000);
fprintf(fid,'\nBad pings Removed: %d',N_badpings);
fprintf(fid,'\n===================================================\n');
fprintf(fid,'%10s %10s %15s %15s %15s %15s \n','Lat','Lon','Range (m)','Residual (s)','Doppler Vel. (m/s)','TWT corr. (ms)');
for ii = 1:Nobs
	fprintf(fid,'%3d: %10f %10f %10f %10f %10f %10f\n',ii,lats_ship(ii),lons_ship(ii),range(ii),dtwt_bs(ii)*1000,vr_bs(ii),dtwtcorr_bs(ii)*1000);
end
fclose(fid);

if ifsave && ifplot
% Save plots
if ~exist([modified_outdir,'/plots/'])
	mkdir([modified_outdir,'/plots/']);
end
save2pdf([modified_outdir,'/plots/',data.sta,'_1_OBSlocation.pdf'],f1,500)
save2pdf([modified_outdir,'/plots/',data.sta,'_2_misfit.pdf'],f2,500)
save2pdf([modified_outdir,'/plots/',data.sta,'_3_VelCorrs.pdf'],f3,500)
save2pdf([modified_outdir,'/plots/',data.sta,'_4_bootstrap.pdf'],f100,500)
save2pdf([modified_outdir,'/plots/',data.sta,'_5_Ftest.pdf'],f101,500)
end

if ifsave
    datamat.sta = data.sta;
	datamat.drop_lonlatz = [data.lon_drop,data.lat_drop,data.z_drop];
	datamat.lons_ship = lons_ship;
	datamat.lats_ship = lats_ship;
	datamat.x_ship = x_ship;
	datamat.y_ship = y_ship;
	datamat.z_ship = z_ship;
	datamat.v_ship = v_ship;
	datamat.dtwt_bs = dtwt_bs;
	datamat.twtcorr_bs = twtcorr_bs;
	datamat.dtwtcorr_bs = dtwtcorr_bs;
	datamat.lon_sta_bs = lon_sta;
	datamat.lat_sta_bs = lat_sta;
	datamat.x_sta_bs = x_sta;
	datamat.y_sta_bs = y_sta;
	datamat.z_sta_bs = z_sta;
	datamat.drift_bs = drift;
	datamat.azi_bs = azi;
	datamat.TAT_bs = TAT;
	datamat.V_w_bs = V_w;
	datamat.E_rms = E_rms;
    datamat.Ftest_res = Ftest_res;
	datamat.loc_xyz = [mean(x_sta),mean(y_sta),mean(z_sta)];
	datamat.loc_lolaz = [mean(lon_sta),mean(lat_sta),mean(z_sta)];
	datamat.mean_drift_az = [mean(drift) r2d(mean_ang(d2r(azi)))];
	if ~exist([modified_outdir,'/mats'])
		mkdir([modified_outdir,'/mats']);
	end
	save([modified_outdir,'/mats/',data.sta,'_data.mat'],'datamat');
end


end

% message if no success
if is==Nstas && ~exist('rawdatfile','var')
    fprintf('*********************\nStation %s does not seem to exist in that folder\n',onesta);
end
cd (wd)
