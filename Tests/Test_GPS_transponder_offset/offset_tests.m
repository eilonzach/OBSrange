addpath('../../OBSrange_v1_MATLAB/functions/')
addpath('../../synthetics/')
%% Some simple tests of location
% apparent location (GPS location)
Eg = 100;
Ng = 100;

% playing with these shows that this works perfectly
dforward = -10;
dstarboard = -10;
cog = 340;

[dE,dN] = GPS_transp_correction(dforward,dstarboard,cog);

Et = Eg + dE;
Nt = Ng + dN;

figure(1); clf; hold on
plot(Eg,Ng,'*b','markersize',15);
plot(Et,Nt,'^r','markersize',15);
set(gca,'xlim',[80,120],'ylim',[80,120]);

%% synthetic survey (copied from synth_survey.m) with/without offset
water_depth = 5; % in km
drop_location = [-7.54 -133.62 water_depth]; % [lat,lon,z]
noise = 0.004; %0.004; % std of timing error

datN = 10; 
obs_location_xyz = [+.1 +.3 5.050]; % x,y,z (km)
vp_actual = 1.520; % km/s
tat = 0.014; %s

%% survey parameters
survey = 'PACMAN'; % 'PACMAN' or 'circle' or 'diamond' or 'tri_edge' or 'tri_center' or 'cross(2)' or 'line(2)'
radius = 1.0; % radius of survey, in Nm
fprintf('Survey %s rad=%.00f\n',survey,radius);
survstart = now;
survey_dt = max([10,60*radius/1.3]); % time lapse of ranging

% system/default parameters
nm2km = 1.852;
obs_default_xyz = [ 0 0 water_depth ]; % x,y,z (km)
vp_default = 1.5; % km/s
tat_default = 0.013; %s

%% set up geometry and project to map coords
% centpt = drop_location(1:2); % [lat lon]
az = 0; % azimuth of +y direction
% proj.origin  = centpt;   % model origin
% mstruct = defaultm('mercator');
% mstruct.origin = [proj.origin -az];
% mstruct = defaultm( mstruct ); 
% proj.map_proj = mstruct;

if niter==1
% obs location
% [obs_location_laloz(1),obs_location_laloz(2)] = project_xy(proj,obs_location_xyz(1),obs_location_xyz(2),'inverse');
[obs_location_laloz(2),obs_location_laloz(1)] = xy2lonlat_nomap(drop_location(2),drop_location(1),1e3*obs_location_xyz(1),1e3*obs_location_xyz(2));

obs_location_laloz(3) = obs_location_xyz(3);

fprintf('True OBS location: \nlat = %.5f \nlon = %.5f \ndep = %.4f\n',obs_location_laloz)

end
%% make survey
dl = 0.005;
[ fsurvx,fsurvy ] = synth_survey_pattern( survey, radius,dl );

% add noise + smooth to round the corners (as the ships do..
fsurvx =  moving_average(fsurvx + normrnd(0,1e-2,size(fsurvx)),50); 
fsurvy =  moving_average(fsurvy + normrnd(0,1e-2,size(fsurvy)),50);

fsurvl = [0;cumsum(sqrt(diff(fsurvx).^2 + diff(fsurvy).^2))];

Nf = length(fsurvx);
% ship speed at each point - keep constant
vship_kn = 8*[1,1];
            

vship_ms = vship_kn*1.852*1000/60/60; % in m/s
% noisify ship speed
% vship_ms = vship_ms + normrnd(0,0.05,size(vship_ms));
% times at each point on survey
fsurvt = [0;cumsum(1000*diff(fsurvl)./midpts(vship_ms))];
clf, plot(fsurvt/3600); axis equal

% interogate every survey_dt seconds
send_survt = [zeros(4,1);[0:survey_dt:fsurvt(end)-10]'];
% + add in some random extra times
send_survt = sort([send_survt;fsurvt(end)*rand(3,1)]);

%% loop and do forward model  different OBS locations and different data shadows
if niter > 1
    data = struct('survey',survey,'radius',radius,'rmsnoise',noise,...
                  'drop',drop_location,'vship_kn',unique(vship_kn),'dt_survey',survey_dt,...
                  'survx',fsurvx,'survy',fsurvy,'survt',fsurvt,...
                  'obs_loc_xyz',[],'obs_loc_laloz',[],...
                  'TAT',[],'Vp_water',[],...
                  'Nobs',[],'survlats',[],'survlons',[],'survts',[],'tot_dt',[],...
                  'corr_dt',[],'v_surv_true',[]);
end
                             

for jj = 1:niter
    if jj == 1 || rem(jj,500)==0,fprintf('Iteration %.0f\n',jj);end

% randomise model parameters if bootstrapping
if niter > 1
    obs_location_xyz = obs_default_xyz + [normrnd(0,x_std),normrnd(0,y_std),normrnd(0,z_std)];
    tat = tat_default + normrnd(0,tat_std);
    vp_actual = vp_default + normrnd(0,vp_std);
    [obs_location_laloz(2),obs_location_laloz(1)] = xy2lonlat_nomap(drop_location(2),drop_location(1),1e3*obs_location_xyz(1),1e3*obs_location_xyz(2));
    obs_location_laloz(3) = obs_location_xyz(3);
    % save these true values
    data(jj).obs_loc_xyz = obs_location_xyz;
    data(jj).obs_loc_laloz = obs_location_laloz;
    data(jj).TAT = tat;
    data(jj).Vp_water = vp_actual;
end

send_survx = linterp(fsurvt,fsurvx,send_survt);
send_survy = linterp(fsurvt,fsurvy,send_survt);
send_dr = sqrt( (send_survx-obs_location_xyz(1)).^2 ...
               +(send_survy-obs_location_xyz(2)).^2 ...
               + obs_location_xyz(3).^2);

prelim_dt = 2*send_dr./vp_actual;
shift_dt = prelim_dt;
delt = 10;
% cycle to make time correction as accurate as possible.
while delt>1e-6
rec_survt = send_survt + shift_dt;
rec_survx = linterp(fsurvt,fsurvx,rec_survt);
rec_survy = linterp(fsurvt,fsurvy,rec_survt);
rec_dr = sqrt(  (rec_survx-obs_location_xyz(1)).^2 ...
               +(rec_survy-obs_location_xyz(2)).^2 ...
               + obs_location_xyz(3).^2);
send_dt = send_dr./vp_actual;
rec_dt  = rec_dr./vp_actual;
tot_dt = send_dt+rec_dt + tat;
delt = rms(tot_dt-shift_dt);
shift_dt = tot_dt;
end

corr_dt = rec_dt-send_dt;

% calc instantaneous velocities
v_surv = [1i*(rec_survx-send_survx)+(rec_survy-send_survy)]./(send_dt+rec_dt); % in m/s
[abs(v_surv)*1000,r2d(angle(v_surv)) + az];
v_surv_true = abs(v_surv).*1000.*[cosd(r2d(angle(v_surv)) + az),sind(r2d(angle(v_surv)) + az)]; % in m/s

% add a little noise
tot_dt = tot_dt + normrnd(0,noise,size(tot_dt));

% assign rec_survx/y/t to survx/y/t
survx = rec_survx;
survy = rec_survy;
survt = rec_survt/24/60/60 + survstart;

% project back to lon,lat
% [survlat,survlon] = project_xy(proj,survx,survy,'inverse');
[ survlon, survlat ] = xy2lonlat_nomap( drop_location(2),drop_location(1), survx*1.e3, survy*1.e3);
[survgc,survaz] = distance(drop_location(1), drop_location(2),survlat,survlon);


%% randomly drop out some points
okpt = rand(size(survlat))>0.2;

% drop out all within some distance of 3 bazs (if not line...
% if ~any(regexp(survey,'line'))
% badaz = 360*rand(3,1);
% badazw = abs(normrnd(0,20,3,1)); % width of drop-out bin
% for ibz = 1:3
%     okpt(ang_diff(survaz,badaz(ibz))<badazw(ibz)) = false;
% end
% end
% reinstate points immediately on top of (<100m from) station
okpt(survgc*111.1<0.1) = true;
% kill not-ok points
tot_dt(~okpt) = nan;


%% save synth data values
data(jj).Nobs = length(survlon);
data(jj).survlats = survlat;
data(jj).survlons = survlon;
data(jj).survts = survt;
data(jj).tot_dt = tot_dt;
data(jj).corr_dt = corr_dt;
data(jj).v_surv_true = v_surv_true;


end % loop on iterations
