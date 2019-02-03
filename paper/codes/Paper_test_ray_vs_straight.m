clear all
close all
%% parameters 
dr = 0.0001; % for ray tracing - ray discretisation. Units of km
vdz = 0.0001; % for ray tracing - depth discretisation. Units of km
nm2km = 1.852;
xmax = 1.5; % units of nm

% profile location
glat = 46;
glon = -133;
% glat = -6.29;
% glon = -131.91;
gmonth = 0;

% Depths
zmaxs = flip([0.5:0.5:5.5]); %[5,2,0.5]; % units of km

% Ray parameters
ps = 0.01:0.01:1.5;

ifsave = true;
ofile = '../figs/FigureS11.jpg';

%% prep plot
figure(12); clf, set(gcf,'pos',[1151 40 861 960])
ax1 = axes('pos',[0.09    0.4100    0.61    0.5500]); hold on
ax2 = axes('pos',[0.09    0.2300    0.61    0.1577]); hold on
ax3 = axes('pos',[0.73    0.4100    0.23    0.5500]); hold on
ax4 = axes('pos',[0.09    0.0500    0.61    0.1577]); hold on

% cols = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];
cols = jet(length(zmaxs));


%% Load V profile
addpath ../../OBSrange_v2_MATLAB_clean/functions/ocean_profiles_obsrange
wd = pwd;
cd ../../OBSrange_v2_MATLAB_clean/functions/ocean_profiles_obsrange
[v_profile.ssp,v_profile.z]=getlev_obsrange(glat,glon,gmonth);
cd(wd)

v_profile.ssp = v_profile.ssp/1000;
v_profile.z = v_profile.z/1000;

plot(ax3,v_profile.ssp,v_profile.z,'b','linewidth',2)
% plot(ax3,V_xbt(:,3)/1000,V_xbt(:,1)/1000,'linewidth',2)


%% loop over depths
Nz = length(zmaxs);

%% Loop over ray parameters
Np = length(ps);

plot(ax3,v_profile.ssp,v_profile.z,'b','linewidth',2)
for ip=1:Np
    p = ps(ip);    
    % compute rays
    try
        [rx,rz,Dr,rt,rv] = shootrays(p, v_profile, max(zmaxs), dr,vdz);
    end
    %plot true rays
    plot(ax1,rx/nm2km,rz*1000,'-k','linewidth',1.5); hold on;
    
    for iz = 1:Nz
        zmax = zmaxs(iz);
        %% Ray vs. straight-line approx
        % length of respective rays
        Dr_ray_z(ip,iz) = interp1(rz,Dr,zmax); % bent ray
        Dx_ray_z(ip,iz) = interp1(rz,rx,zmax); % bent ray
        t_ray(ip,iz) = interp1(rz,rt,zmax); % bent ray
        Sr_apx(ip,iz) = sqrt(Dx_ray_z(ip,iz).^2 + zmax.^2); % straight line start to end
        t_apx(ip,iz) = Sr_apx(ip,iz) / harmmean(rv(rz<=zmax));
        if Dx_ray_z(ip,iz)/nm2km>xmax, continue; end
        % plot approx rays
        plot(ax1,[0;Dx_ray_z(ip,iz)]/nm2km,[0;zmaxs(iz)]*1000,'color',cols(iz,:),'linewidth',1.5);
        
        % traveltimes
        plot(ax2,Dx_ray_z(:,iz)/nm2km,t_ray(:,iz)*2,'-','color',cols(iz,:),'linewidth',1.5); hold on;
        plot(ax2,Dx_ray_z(:,iz)/nm2km,t_apx(:,iz)*2,'^','color',cols(iz,:));     
    end
end
% difference in ray travel time (ms) (one-way)
dtime_ms = 1000*(t_ray-t_apx);
for iz = 1:Nz
        plot(ax4,Dx_ray_z(:,iz)/nm2km,2*dtime_ms(:,iz),'-o','color',cols(iz,:),'linewidth',1.5);  
end

set(ax1,'ydir','rev','xaxislocation','top','xticklabel','','xlim',[0 xmax],'ylim',[0 1000*max(zmaxs)],...
    'box','on','layer','top','linewidth',2,'fontsize',14)
set(ax2,'xlim',[0 xmax],'ylim',[0 max(t_ray(Dx_ray_z<xmax*nm2km))*2],...
    'box','on','layer','top','linewidth',2,'XTickLabel',[],'fontsize',14)
set(ax3,'ydir','rev','yticklabels','',...
    'ylim',[0 max(zmaxs)],'xlim',[1.47 1.55],...
    'box','on','layer','top','linewidth',2,'fontsize',14)    
set(ax4,'xlim',[0 xmax],'ylim',[min(dtime_ms(Dx_ray_z<xmax*nm2km))*2 max(dtime_ms(Dx_ray_z<xmax*nm2km))*2],...
    'box','on','layer','top','linewidth',2,'fontsize',14)

xlabel(ax4,'Horizontal distance from ship (Nm)','fontweight','bold')
xlabel(ax3,'Soundspeed (km/s)','fontweight','bold')

ylabel(ax1,'Depth (m)','fontweight','bold')
ylabel(ax2,'2-way times (s)','fontweight','bold')
ylabel(ax4,{'\delta T_{bend}','Ray - Line (ms)'},'fontweight','bold')


if ifsave
    save2jpg(12,ofile)
end