addpath('../../OBSrange_v1_MATLAB_clean/functions/')
addpath('../../synthetics/')
%% Some simple tests of location
% apparent location (GPS location)
Eg = 100;
Ng = 100;

% playing with these shows that this works perfectly
dforward = -10;
dstarboard = -10;
cog = 340;

[dE,dy] = GPS_transp_correction(dforward,dstarboard,cog);

Et = Eg + dE;
Nt = Ng + dy;

figure(1); clf; hold on
plot(Eg,Ng,'*b','markersize',15);
plot(Et,Nt,'^r','markersize',15);
set(gca,'xlim',[80,120],'ylim',[80,120]);

%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
%% synthetic survey (copied from synth_survey.m) with/without offset
water_depth = 5; % in km
drop_location = [-7.54 -133.62 water_depth]; % [lat,lon,z]
noise = 0.004; %0.004; % std of timing error

obs_location_xyz = [.1 .055 water_depth+0.00]; % x,y,z (km)
vp_actual = 1.51; % km/s
tat = 0.014; %s

vship_kn = 8; % constant ship velocity

dforward = 10; % m
dstarboard = -10; % m

nm2km = 1.852;

niter = 1;

fprintf('\n\nTRUE\n X=%5.1fm Y=%5.1fm Z = %6.1fm, dVp=%4.1fm/s\n',...
         obs_location_xyz*1e3,(vp_actual-1.5)*1e3)

%% start loop on survey, random noise and semi-random patterns
for iter = 1:niter
    iter

%% survey parameters
survey = 'PACMAN'; % 'PACMAN' or 'circle' or 'diamond' or 'tri_edge' or 'tri_center' or 'cross(2)' or 'line(2)'
radius = 1.0; % radius of survey, in Nm
fprintf('Survey %s rad=%.00f\n',survey,radius);
survstart = now;
survey_dt = max([10,60*radius/1.3]); % time lapse of ranging

%% set up geometry and project to map coords
% centpt = drop_location(1:2); % [lat lon]
az = 0; % azimuth of +y direction
% proj.origin  = centpt;   % model origin
% mstruct = defaultm('mercator');
% mstruct.origin = [proj.origin -az];
% mstruct = defaultm( mstruct ); 
% proj.map_proj = mstruct;

%% make survey
dl = 0.005;
[ fsurvx,fsurvy ] = synth_survey_pattern( survey, radius,dl );

% add noise + smooth to round the corners (as the ships do..
fsurvx =  moving_average(fsurvx + normrnd(0,1e-2,size(fsurvx)),50); 
fsurvy =  moving_average(fsurvy + normrnd(0,1e-2,size(fsurvy)),50);

fsurvl = [0;cumsum(sqrt(diff(fsurvx).^2 + diff(fsurvy).^2))];

Nf = length(fsurvx);
% ship speed at each point - keep constant
vship_kn = vship_kn.*[1,1];
            
vship_ms = vship_kn*nm2km*1000/60/60; % in m/s
% noisify ship speed
% vship_ms = vship_ms + normrnd(0,0.05,size(vship_ms));
% times at each point on survey
fsurvt = [0;cumsum(1000*diff(fsurvl)./midpts(vship_ms))];
fsurvsog = sqrt(diff(fsurvx).^2 + diff(fsurvy).^2)./diff(fsurvt)*1000; % in m/s
fsurvcog = atan2d(diff(fsurvx),diff(fsurvy));
% make these just as long as position vectors by extending
fsurvsog = midpts(fsurvsog([1,1:end,end]));
fsurvcog = midpts(fsurvcog([1,1:end,end]));

clf, plot(fsurvt/3600); axis equal

% interogate every survey_dt seconds
send_survt = [zeros(4,1);[0:survey_dt:fsurvt(end)-10]'];
% + add in some random extra times
send_survt = sort([send_survt;fsurvt(end)*rand(3,1)]);


%% loop and do forward model  different OBS locations and different data shadows                        
[ obs_location_laloz(2), obs_location_laloz(1) ] = xy2lonlat_nomap( drop_location(2),drop_location(1), obs_location_xyz(1), obs_location_xyz(2));

send_survx_gps = linterp(fsurvt,fsurvx,send_survt);
send_survy_gps = linterp(fsurvt,fsurvy,send_survt);
send_surv_cog = linterp(fsurvt,fsurvcog,send_survt);

% account for gps-transp offset
[dE,dy] = GPS_transp_correction(dforward,dstarboard,send_surv_cog);
send_survx = send_survx_gps + dE/1000;
send_survy = send_survy_gps + dy/1000;

send_dr = sqrt( (send_survx-obs_location_xyz(1)).^2 ...
               +(send_survy-obs_location_xyz(2)).^2 ...
               + obs_location_xyz(3).^2);

prelim_dt = 2*send_dr./vp_actual;
shift_dt = prelim_dt;
delt = 10;
% cycle to make time correction as accurate as possible.
while delt>1e-6
rec_survt = send_survt + shift_dt;
rec_survx_gps = linterp(fsurvt,fsurvx,rec_survt);
rec_survy_gps = linterp(fsurvt,fsurvy,rec_survt);
rec_surv_cog = linterp(fsurvt,fsurvcog,rec_survt);
% account for gps-transp offset
[dE,dy] = GPS_transp_correction(dforward,dstarboard,rec_surv_cog);
rec_survx = rec_survx_gps + dE/1000;
rec_survy = rec_survy_gps + dy/1000;

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
sog = abs(v_surv).*1000;
cog = r2d(angle(v_surv)) + az;
survsog = sog.*[cosd(cog),sind(cog)]; % in m/s [N,E]

% add a little noise
tot_dt = tot_dt + normrnd(0,noise,size(tot_dt));

% assign rec_survx/y/t to survx/y/t
survx = rec_survx_gps;
survy = rec_survy_gps;
survt = rec_survt/24/60/60 + survstart;

% project back to lon,lat
% [survlat,survlon] = project_xy(proj,survx,survy,'inverse');
[ survlon, survlat ] = xy2lonlat_nomap( drop_location(2),drop_location(1), survx*1.e3, survy*1.e3);
% [survgc,survaz] = distance(drop_location(1), drop_location(2),survlat,survlon);


%% randomly drop out some points
okpt = rand(size(survlat))>0.2;

% drop out all within some distance of 3 bazs (if not line...
if ~any(regexp(survey,'line'))
badaz = 360*rand(3,1);
badazw = abs(normrnd(0,20,3,1)); % width of drop-out bin
for ibz = 1:3
    okpt(ang_diff(survaz,badaz(ibz))<badazw(ibz)) = false;
end
end
% reinstate points immediately on top of (<100m from) station
okpt(survgc*111.1<0.1) = true;
% kill not-ok points
tot_dt = tot_dt(okpt);
survlat = survlat(okpt);
survlon = survlon(okpt);
survt = survt(okpt);
survsog = survsog(okpt);

%% data values in synthetic dataset (i.e. available to algorithm)
survlat;
survlon;
survt;
tot_dt;
% and
corr_dt;
survsog;

%% plots of send, receive points
if iter==1
figure(4), clf; hold on
plot(send_survx_gps,send_survy_gps,'og',...
     send_survx,send_survy,'*g',...     
     rec_survx_gps,rec_survy_gps,'or',...
     rec_survx,rec_survy,'*r'...     
     )
 
%  figure(1);
% plot(fsurvx,fsurvy)
% 
% figure(2);
% plot(send_survx,send_survy,'xr',rec_survx,rec_survy,'xb'); axis equal
% 
% 
figure(3), clf; hold on
scatter(survlon,survlat,60,tot_dt,'filled');
% plot(survlon(~okpt),survlat(~okpt),'xk','markersize',8,'linewidth',2);
plot(obs_location_laloz(2),obs_location_laloz(1),'*r','markersize',14)
end
 
%% Estimates of sog, cog for corrections a posteriori (based only on data)
% Parameters

par = struct([]);
par(1).vp_w = 1500; % Assumed water velocity (m/s)
par.E_thresh = 1e-8; % RMS reduction threshold for inversion
par.if_twtcorr = 0;
par.npts_movingav = 1; %5; % number of points to include in moving average smoothing of ship velocity (1 = no smoothing);
par.TAT = tat; % (s)
% Norm damping for each model parameter (damping towards starting model)
% Larger values imply more damping towards the starting model.
par.dampx = 0;
par.dampy = 0;
par.dampz = 0; %1e3; %0
par.dampdvp = 5e-8; %1e3; %5e-8
H = diag([par.dampx, par.dampy, par.dampz, par.dampdvp]);
% Global norm damping for stabilization
par.epsilon = 1e-10;
m0_strt = [0,0,water_depth*1000,0]';
% shots
[ x_ship, y_ship ] = lonlat2xy_nomap( drop_location(2), drop_location(1), survlon, survlat );
z_ship = zeros(size(x_ship));
t_ship = (survt - survt(1))*24*60*60;
% Calculate velocity of ship
v_ship = pt_veloc( x_ship, y_ship, z_ship, t_ship );
v_ship = [moving_average(v_ship(1,:),par.npts_movingav)'; moving_average(v_ship(2,:),par.npts_movingav)'; moving_average(v_ship(3,:),par.npts_movingav)'];

% dumb estimate
par.if_twtcorr = 0;
[ m_final_1(:,iter),models,v,N,R,Cm ] = inv_newtons( par,m0_strt,tot_dt,...
                    x_ship,y_ship,z_ship,...
                    v_ship,H);
if niter==1
fprintf('Dumb estimate\n X=%5.1fm Y=%5.1fm Z = %6.1fm, dVp=%4.1fm/s\n',...
         m_final_1(:,iter))
end
     
% estimate accounting for doppler
par.if_twtcorr = 1;
[ m_final_2(:,iter),models,v,N,R,Cm ] = inv_newtons( par,m0_strt,tot_dt,...
                    x_ship,y_ship,z_ship,...
                    v_ship,H);
if niter==1
fprintf('w/ Doppler\n X=%5.1fm Y=%5.1fm Z = %6.1fm, dVp=%4.1fm/s\n',...
         m_final_2(:,iter))
end
     
% estimate accounting for GPS offset
par.if_twtcorr = 1;
survcog = atan2d(v_ship(1,:),v_ship(2,:));
[dx,dy] = GPS_transp_correction(dforward,dstarboard,survcog');
x_ship = x_ship + dx;
y_ship = y_ship + dy;

[ m_final_3(:,iter),models,v,N,R,Cm ] = inv_newtons( par,m0_strt,tot_dt,...
                    x_ship,y_ship,z_ship,...
                    v_ship,H);
if niter==1
fprintf('w/ Doppler+GPSoffset\n X=%5.1fm Y=%5.1fm Z = %6.1fm, dVp=%4.1fm/s\n',...
         m_final_3(:,iter))
end


end; clear iter % loop on iter

