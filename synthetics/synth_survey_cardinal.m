clear
close all
% special version of synthetic survey creation script that just has the
% ship hold station at 4 cardinal points and over the drop location. 
addpath('../OBSrange_v1_MATLAB_clean/functions');

%% TRUE VALUES
for water_depth = [5,2,0.5] % in km
drop_location = [-7.54 -133.62 water_depth]; % [lat,lon,z]
noise = 0.004; %0.004; % std of timing error

dforward = 10; % in m
dstarboard = 10; % in m
gps_offset_str = 'fr10';

niter = 1e4; %1;%1e4; % if niter>1, will not make plots or save output file in SIO format

ifsave = false;

%% system/default parameters
nm2km = 1.852;
obs_default_xyz = [ 0 0 water_depth ]; % x,y,z (km)
vp_default = 1.5; % km/s
tat_default = 0.013; %s

%% survey parameters
survey = 'cardinal';
for radius = [1]; % radius of survey, in Nm
fprintf('Survey %s rad=%4.2f depth=%4.0f\n',survey,radius,water_depth*1e3);
survstart = now;
survey_dt = max([10,60*radius/1.3]); % time lapse of ranging

%% model parameters - if one iteration
% % syn10
% datN = 10; 
% obs_location_xyz = [+.1 +.3 5.050]; % x,y,z (km)
% vp_actual = 1.520; % km/s
% tat = 0.014; %s

% % syn11
% datN = 11; 
% obs_location_xyz = [+.2 -.3 5.100]; %[+.1 +.3 5.050]; % x,y,z (km)
% vp_actual = 1.530; %1.520; % km/s
% tat = 0.0145; %0.014; %s

% syn12
datN = 12; 
obs_location_xyz = [+.2 -.4 5.050]; %[+.1 +.3 5.050]; % x,y,z (km)
vp_actual = 1.520; %1.520; % km/s
tat = 0.014; %0.014; %s

%% model parameters - if doing many iteration
x_std = 0.100; % in km
y_std = 0.100; % in km
z_std = 0.050; % in km
vp_std = 0.010; % in km/s
tat_std = 3e-3; %1e-3; %3e-3;


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

%% loop and do forward model  different OBS locations and different data shadows
if niter > 1
    data = struct('survey',survey,'radius',radius,'rmsnoise',noise,...
                  'drop',drop_location,'vship_kn',0,'dt_survey',survey_dt,...
                  'TG_dforward',dforward,'TG_dstarboard',dstarboard,...
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

%% make survey
% send
send_survx_gps = radius*nm2km*[zeros(1,4),zeros(1,4),ones(1,4),zeros(1,4),-ones(1,4)]';
send_survy_gps = radius*nm2km*[zeros(1,4),ones(1,4),zeros(1,4),-ones(1,4),zeros(1,4)]';
Nping = length(send_survx_gps);
send_surv_cog = 120*ones(Nping,1) + normrnd(0,1,[Nping 1]);
send_survx_gps = send_survx_gps + normrnd(0,1e-2,[Nping 1]);
send_survy_gps = send_survy_gps + normrnd(0,1e-2,[Nping 1]);
    % account for gps-transp offset
    [dx,dy] = GPS_transp_correction(dforward,dstarboard,send_surv_cog);
    send_survx = send_survx_gps + dx/1000;
    send_survy = send_survy_gps + dy/1000;
send_dr = sqrt( (send_survx-obs_location_xyz(1)).^2 ...
               +(send_survy-obs_location_xyz(2)).^2 ...
               + obs_location_xyz(3).^2);
% receive
rec_survx_gps = send_survx_gps + normrnd(0,2e-3,[Nping 1]);
rec_survy_gps = send_survy_gps + normrnd(0,2e-3,[Nping 1]);
rec_surv_cog  = send_surv_cog + normrnd(0,1,[Nping 1]);
    % account for gps-transp offset
    [dx,dy] = GPS_transp_correction(dforward,dstarboard,rec_surv_cog);
    rec_survx = rec_survx_gps + dx/1000;
    rec_survy = rec_survy_gps + dy/1000;
rec_dr = sqrt(  (rec_survx-obs_location_xyz(1)).^2 ...
               +(rec_survy-obs_location_xyz(2)).^2 ...
               + obs_location_xyz(3).^2);
send_dt = send_dr./vp_actual;
rec_dt  = rec_dr./vp_actual;
tot_dt = send_dt+rec_dt + tat;

corr_dt = rec_dt-send_dt;

% calc instantaneous velocities
v_surv = [1i*(rec_survx-send_survx)+(rec_survy-send_survy)]./(send_dt+rec_dt); % in m/s
[abs(v_surv)*1000,r2d(angle(v_surv)) + az];
sog = abs(v_surv).*1000; % m/s
cog = send_surv_cog;
surv_vel_true = sog.*[cosd(cog),sind(cog)]; % in m/s [N,E]

% add a little noise
tot_dt = tot_dt + normrnd(0,noise,size(tot_dt));

% assign rec_survx/y/t to survx/y/t
survx = rec_survx_gps;
survy = rec_survy_gps;
survt = now + [1:Nping];

% project back to lon,lat
% [survlat,survlon] = project_xy(proj,survx,survy,'inverse');
[ survlon, survlat ] = xy2lonlat_nomap( drop_location(2),drop_location(1), survx*1.e3, survy*1.e3);
[survgc,survaz] = distance(drop_location(1), drop_location(2),survlat,survlon);


%% save synth data values
data(jj).Nobs = length(survlon);
data(jj).survlats = survlat;
data(jj).survlons = survlon;
data(jj).tot_dt = tot_dt;
data(jj).corr_dt = corr_dt;
data(jj).v_surv_true = surv_vel_true;
data(jj).cog_true = cog;


end % loop on iterations


%% some plots to check
if niter==1
figure(1);
plot(fsurvx,fsurvy)

figure(2);
plot(send_survx,send_survy,'xr',rec_survx,rec_survy,'xb'); axis equal


figure(3), clf; hold on
scatter(survlon,survlat,60,tot_dt,'filled');
plot(survlon(~okpt),survlat(~okpt),'xk','markersize',8,'linewidth',2);
plot(obs_location_laloz(2),obs_location_laloz(1),'*r','markersize',14)

end

% profile viewer

%% output file
if ifsave
    sta = sprintf('syn%u',datN);

% many iterations - output data structure
if niter>1
    save(sprintf('synth_surveys_paper/SynthBoot_%s_rad%.2f_z%.0fm_%s.mat',survey,radius,water_depth*1e3,gps_offset_str),'data');
end

if niter==1
% single station - make output file to mimic SIO survey files
survlah = ['N','S'];
survloh = ['E','W'];

fid = fopen([sta,'.txt'],'w');
% header like real files
fprintf(fid,'Ranging data taken on:  %s.%s\n',datestr(survstart,31),datestr(survstart,'fff'));
fprintf(fid,'Cruise:                 synthetic_survey_%s\n',survey);
fprintf(fid,'Site:                   %s\n',sta);
fprintf(fid,'Instrument:\n');
fprintf(fid,'Drop Point (Latitude):  %.5f\n',drop_location(1));
fprintf(fid,'Drop Point (Longitude): %.5f\n',drop_location(2));
fprintf(fid,'Depth (meters):         %.5f\n',drop_location(3)*1e3);
fprintf(fid,'Comment:\n');
fprintf(fid,'==================================================\n\n');

for ii = 1:length(survx)
    if okpt(ii)
        fprintf(fid,'%4.0f msec. Lat: %1.0f %07.4f %1s  Lon: %1.0f %07.4f %1s  Alt: %5.2f Time(UTC): %s:%3.0f:%s\n',...
                tot_dt(ii)*1000,...
                abs(fix(survlat(ii))),60*abs(rem(survlat(ii),1)),survlah(1.5-0.5*sign(survlat(ii))),...
                abs(fix(survlon(ii))),60*abs(rem(survlon(ii),1)),survloh(1.5-0.5*sign(survlon(ii))),...
                25 + normrnd(0,0.5),...
                datestr(survt(ii),'yyyy'),doy(datestr(survt(ii),'yyyy'),datestr(survt(ii),'mm'),datestr(survt(ii),'dd')),...
                datestr(survt(ii),'HH:MM:SS'));
    else
        fprintf(fid,'Event skipped - bad luck\n');
    end
end
fclose(fid);

save(sprintf('trudata_%s',sta));
% copyfile(sprintf('SynthSurvey%.0f.mat',datN),'synth_surveys_paper/');
end             

end





end
end

