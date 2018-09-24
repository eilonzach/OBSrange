% Plot geostropic flow data
clear; close all;

datadir = './data/dataset-duacs-nrt-global-merged-allsat-phy-l4-v3/unzipped/';
outdir = './matfigs/';
if ~exist(outdir)
    mkdir(outdir);
end

fils = dir([datadir,'*.nc']);

ui_sum = 0;
vi_sum = 0;
hi_sum = 0;
for ifil = 1:length(fils)-6 %[1 8]
    
    filename = fils(ifil).name;
    datapath = [datadir,filename];
    
    strs = strsplit(strtok(filename,'.'),'_');
    startdate{ifil} = strs{end-1};
    enddate{ifil} = strs{end};

    ugos = ncread(datapath,'ugos')'; % N-S velocity
    vgos = ncread(datapath,'vgos')'; % E-W velocity
    h = ncread(datapath,'adt')'; % Dynamic Sea Level
    lat = ncread(datapath,'latitude');
    lon = ncread(datapath,'longitude');

    [LON, LAT] = meshgrid(lon,lat);
    % LON(LON>180) = LON(LON>180)-360;
    
    %% Plot
    lats = [-9 -3];
    lons = [360-136 360-130];
    ilon = lon>=lons(1)-1 & lon<=lons(2)+1;
    ilat = lat>=lats(1)-1 & lat<=lats(2)+1;

    LONi = LON(ilat,ilon);
    LATi = LAT(ilat,ilon);
    hi = h(ilat,ilon);
    ui = vgos(ilat,ilon);
    vi = ugos(ilat,ilon);

    fig1 = figure(1); clf;
    contourf(LONi,LATi,hi,20,'linestyle','none'); hold on;
    quiver(LONi,LATi,vi,ui,'k','filled','linewidth',1)
    xlim(lons);
    ylim(lats);
    ax = gca;
    xlabel('Longitude');
    ylabel('Latitude');
    title([startdate{ifil},'-',enddate{ifil}]);
    set(gca,'linewidth',1.5,'fontsize',16, ...
        'XTickLabel',num2str(str2num(cell2mat(ax.XTickLabel))-360));
    % ax = gca
    cb = colorbar('linewidth',1.5,'fontsize',16);
    ylabel(cb,'Dynamic Sea Level (m)');
    % colormap('inferno');
    colormap('parula');

    save2pdf([outdir,filename,'.pdf'],fig1,500);
    
    ui_sum = ui_sum + ui;
    vi_sum = vi_sum + vi;
    hi_sum = hi_sum + hi;
end
ui_mean = ui_sum / ifil;
vi_mean = vi_sum / ifil;
hi_mean = hi_sum / ifil;
%% Plot deployment average
fig2 = figure(2); clf;
contourf(LONi,LATi,hi_mean,20,'linestyle','none'); hold on;
quiver(LONi,LATi,vi_mean,ui_mean,'k','filled','linewidth',1)
xlim(lons);
ylim(lats);
ax = gca;
xlabel('Longitude');
ylabel('Latitude');
title([startdate{1},'-',enddate{end}]);
set(gca,'linewidth',1.5,'fontsize',16, ...
    'XTickLabel',num2str(str2num(cell2mat(ax.XTickLabel))-360));
% ax = gca
cb = colorbar('linewidth',1.5,'fontsize',16);
ylabel(cb,'Dynamic Sea Level (m)');
% colormap('inferno');
colormap('parula');

save2pdf([outdir,'mean_',startdate{1},'_',enddate{end},'.pdf'],fig2,500);