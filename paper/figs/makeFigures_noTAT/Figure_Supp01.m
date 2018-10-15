%% Function to produce Supplemental figure 01
%  Figure S1 geostrophic flow data
function Figure_Supp01
addpath('/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/OBSrange_v1_MATLAB/functions/');
ofile = '../FigureS1';
datadir = '/Users/russell/Lamont/PROJ_OBSrange/geostrophicflow_plots/data/dataset-duacs-nrt-global-merged-allsat-phy-l4-v3/unzipped/';
ifsave = 1;

filenames = {
            'nrt_global_allsat_phy_l4_20180421_20180427.nc';
            'nrt_global_allsat_phy_l4_20180429_20180505.nc';
            };

f20 = figure(20); clf;
set(gcf,'position',[25   242   913   436],'color','w')
ax(1) = subplot(1,2,1); hold(ax(1),'on'); box on;
ax(2) = subplot(1,2,2); hold(ax(2),'on'); box on;

for ifil = 1:length(filenames) %1:length(fils)-6 %[1 8]
    
    filename = filenames{ifil};
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
    
    % cut out area for legend
    if ifil == 1
        minlon = -131+360;
        maxlon = -120+360;
        minlat = -15;
        maxlat = -8.4;
        Ind = LONi>minlon & LONi<maxlon & LATi>minlat & LATi<maxlat;
        LONi(Ind) = nan;
        LATi(Ind) = nan;
        hi(Ind) = nan;
        ui(Ind) = nan;
        vi(Ind) = nan;
    end

    contourf(ax(ifil),LONi,LATi,hi,20,'linestyle','none'); hold on;
    quiver(ax(ifil),LONi,LATi,vi,ui,'k','filled','linewidth',1,'AutoScale','off')
    xlim(ax(ifil),lons);
    ylim(ax(ifil),lats);
    xlabel(ax(ifil),'\bf{Longitude}','interpreter','latex');
    if ifil == 1
        ylabel(ax(ifil),'\bf{Latitude}','interpreter','latex');
    end
    title(ax(ifil),['\bf{',startdate{ifil},'-',enddate{ifil},'}'],'interpreter','latex');
    set(ax(ifil),'linewidth',1.5,'fontsize',14,'TickDir','out','layer','top', ...
        'XTickLabel',num2str(str2num(cell2mat(ax(ifil).XTickLabel))-360));
    % ax = gca
    if ifil == 2
        pos = ax(ifil).Position;
        cb = colorbar(ax(ifil),'linewidth',1.5,'fontsize',14);
        ylabel(cb,'\bf{Dynamic Sea Level (m)}','interpreter','latex');
        ax(ifil).Position = pos;
        set(ax(ifil),'yticklabel',[]);
    end
    % colormap('inferno');
    colormap('parula');
%     axis(ax(ifil),'equal');
    if ifil == 1
        scale_size = 0.5; % m/s
        quiver(ax(ifil),-130.75+360,-8.85,scale_size,0,'k','filled','linewidth',1,'AutoScale','off','MaxHeadSize',1)
        txt = sprintf('%.f cm/s',scale_size*100);
        text(ax(ifil),-130.9+360,-8.6,txt,'fontsize',12);
    end
    
    min_hi(ifil) = min(min(hi));
    max_hi(ifil) = max(max(hi));
    
end

caxis(ax(1),[min(min_hi) max(max_hi)]);
caxis(ax(2),[min(min_hi) max(max_hi)]);

dx_shift = -0.055;
dx = 1.15;
ax(1).Position = [ax(1).Position(1)+dx_shift, ax(1).Position(2), ax(1).Position(3)*dx, ax(1).Position(4)];
dx_shift = -0.07;
ax(2).Position = [ax(2).Position(1)+dx_shift, ax(2).Position(2), ax(2).Position(3)*dx, ax(2).Position(4)];


save2pdf(ofile,f20,500);


end

