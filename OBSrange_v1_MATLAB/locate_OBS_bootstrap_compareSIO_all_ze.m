% Uses two-way travel time information and ship coordinates from an OBS 
% survey to invert for station location on the seafloor (Lat, Lon, Depth),
% turn-around time (TAT), static correction to the sound velocity through the
% water column (dvp), and the velocity of the ship in the radial
% direction of the survey circle (vr0).
%
% Josh R. & Zach E. 4/16/18
%
%
clear; close all;
addpath('./functions');

%% Parameters
% projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA/'; % path to project
% projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_EC03/'; % EC03
% projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_synthtest/'; % SYNTH
% projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_synthtest3/'; % SYNTH 3
projpath = '/Users/russell/Lamont/PROJ_OBSrange/working/OBSrange/projects/PacificORCA_synthtest4_REVISION1_GPScorr/'; % REVISION1
datapath = './'; % path to survey data
datapathSIO = './'; % path to survey data
outdir_OBSrange = './OUT_OBSrange/2_OUT_wcorr_xrec/'; % output directory
% outdir_SIOcomp = './OUT_OBSrange/SIO_compare_wbads/'; % output directory
is_savemat = 1; % Save *.mat file of results?
Nbins = 15; % Bins for histogram plots
is_nobads = 1; % use SIO results with all bad pings removed

if is_nobads
    outdir_SIOcomp = './OUT_OBSrange/8_SIO_compare_nobads/';
else
    outdir_SIOcomp = './OUT_OBSrange/9_SIO_compare_wbads/';
end

wd = pwd;
cd(projpath);
files = dir([datapath,'/*.txt']);
% stas = unique(strtok({files.name},{'_','.txt'}));
stas = {'syn12_z5000m_fr10'};
Nstas = length(stas);
for is = 1:Nstas
    sta = stas{is};
    fprintf('===========================\nComparing OBS %s\n\n',sta);
    %% load all the data
    rawdatafile = sprintf('%s%s.txt',datapath,sta);
    invdatafile = sprintf('%smats/%s_data.mat',outdir_OBSrange,sta);
    data = load_pings(rawdatafile);
    invdata = load(invdatafile); invdata = invdata.datamat;
    if is_nobads
%         dataSIO = loadSIO([datapathSIO,data.sta,'_SIOcorrected_nobads.txt']);
        dataSIO = loadSIO([datapathSIO,data.sta,'_nobads_SIOcorrected.txt']);
    else
        dataSIO = loadSIO([datapathSIO,data.sta,'_SIOcorrected.txt']);
    end
    
    %% basic facts for this station
    % drop point
    olon = data.lon_drop;
    olat = data.lat_drop;
    

%     fprintf('\nlat:   %.5f deg (%f) \nlon:   %.5f deg (%f) \nx:     %f m (%f) \ny:    %f m (%f) \ndepth: %f m (%f) \nTAT:   %f ms (%f) \nv_H20: %f m/s (%f)',mean(lat_sta),std(lat_sta)*2,mean(lon_sta),std(lon_sta)*2,mean(x_sta),std(x_sta)*2,mean(y_sta),std(y_sta)*2,mean(z_sta),std(z_sta)*2,mean(TAT)*1000,std(TAT)*1000*2,mean(v_w),std(v_w)*2);
%     fprintf('\nDrift Lon: %f m (%f) \nDrift Lat: %f m (%f) \nDrift:    %f m (%f) \nDrift Azi: %f deg (%f)\n',mean(dx_drift),std(dx_drift)*2,mean(dy_drift),std(dx_drift)*2,mean(drift),std(drift)*2,mean(azi),std(azi)*2);
%     fprintf('\nRMS:  %f ms (%f)\n',mean(E_rms)*1000,std(E_rms)*2*1000);
    
    %% Compare with SIO
  
    [ dataSIO.x_sta, dataSIO.y_sta ] = lonlat2xy_nomap( olon, olat, dataSIO.lon_sta, dataSIO.lat_sta);
    dloc = sqrt( (dataSIO.x_sta-mean(invdata.x_sta_bs)).^2 + (dataSIO.y_sta-mean(invdata.y_sta_bs)).^2 + (dataSIO.z_sta-mean(invdata.z_sta_bs)).^2 );
    
    % Plot in map view and depth
    f6 = figure(6); clf;
    set(gcf,'position',[53 292 1157 413]);
    % clr = parula(length(models));
    col = lines(length(data.lats));

    subplot(1,2,1);
    plot(data.lons,data.lats,'ok','markerfacecolor',col(2,:),'markersize',13,'linewidth',1); hold on;
%     plot(data_bad.lons,data_bad.lats,'xk','markersize',13,'linewidth',2);
    plot(olon,olat,'sk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1); hold on;
    plot(mean(invdata.lon_sta_bs),mean(invdata.lat_sta_bs),'pk','markerfacecolor',[1 1 0],'markersize',25,'linewidth',1)
    plot(dataSIO.lon_sta,dataSIO.lat_sta,'p','color',col(1,:),'markersize',25,'linewidth',2);
    grid on; box on; set(gca,'fontsize',16,'linewidth',2);
    xlabel('Longitude','fontsize',18);
    ylabel('Latitude','fontsize',18);
    title(['Drift: ',num2str(mean(invdata.drift_bs),'%.2f'),' m (SIO: ',num2str(dataSIO.drift,'%.2f'),' m) dloc=',num2str(dloc,'%.3f'),' m'],'fontsize',18);

    subplot(1,2,2);
    h1(1) = plot(invdata.lats_ship,invdata.z_ship,'ok','markerfacecolor',col(2,:),'markersize',13,'linewidth',1); hold on;
%     plot(data_bad.lats,zeros(size(data_bad.lats)),'xk','markersize',13,'linewidth',2);
    h1(2) = plot(data.lat_drop,data.z_drop,'sk','markerfacecolor',[0.5 0.5 0.5],'markersize',15,'linewidth',1); hold on;
    h1(3) = plot(mean(invdata.lat_sta_bs),mean(invdata.z_sta_bs),'pk','markerfacecolor',[1 1 0],'markersize',25,'linewidth',1);
    h1(4) = plot(dataSIO.lat_sta,dataSIO.z_sta,'p','color',col(1,:),'markersize',25,'linewidth',2);
    grid on; box on; set(gca,'fontsize',16,'linewidth',2);
    xlabel('Latitude','fontsize',18);
    ylabel('Depth','fontsize',18);
    z_diff = mean(invdata.z_sta_bs) - dataSIO.z_sta;
    title(['Z - Z_{SIO} = ',num2str(z_diff,'%.3f'),' m'],'fontsize',18);
    legend(h1,{'Ship','Drop Pt.','Solution','SIO'},'location','southeast','fontsize',13);

    % Plot Misfit
    f7 = figure(7); clf;
    set(gcf,'position',[3         251        1187         447]);
    subplot(1,2,1);
    obs_num = 1:length(invdata.dtwt_bs);
    plot(obs_num,invdata.dtwt_bs*1000,'ok','markerfacecolor',col(2,:),'markersize',13,'linewidth',1); hold on;
    plot([obs_num(1),obs_num(end)],[mean(invdata.dtwt_bs)+std(invdata.dtwt_bs)*2 mean(invdata.dtwt_bs)+std(invdata.dtwt_bs)*2]*1000,'--','color',col(2,:),'linewidth',2); hold on;
    plot([obs_num(1),obs_num(end)],[mean(invdata.dtwt_bs)-std(invdata.dtwt_bs)*2 mean(invdata.dtwt_bs)-std(invdata.dtwt_bs)*2]*1000,'--','color',col(2,:),'linewidth',2); hold on;
    plot(dataSIO.iobs,dataSIO.res*1000,'ok','markerfacecolor',col(1,:),'markersize',13,'linewidth',1); hold on;
    plot([dataSIO.iobs(1),dataSIO.iobs(end)],[mean(dataSIO.res)+std(dataSIO.res)*2 mean(dataSIO.res)+std(dataSIO.res)*2]*1000,'--','color',col(1,:),'linewidth',2); hold on;
    plot([dataSIO.iobs(1),dataSIO.iobs(end)],[mean(dataSIO.res)-std(dataSIO.res)*2 mean(dataSIO.res)-std(dataSIO.res)*2]*1000,'--','color',col(1,:),'linewidth',2); hold on;
    grid on; box on; set(gca,'fontsize',16,'linewidth',2);
    title(['RMS: ',num2str(mean(invdata.E_rms)*1000,'%.2f'),' ms (SIO: ',num2str(dataSIO.E_rms*1000,'%.2f'),' ms)'],'fontsize',18);
    xlabel('Observation #','fontsize',18);
    ylabel('Residuals (ms)','fontsize',18);


    [ Ncount_Erms, cent_Erms ] = plot_hist(subplot(122),invdata.E_rms*1000,Nbins); hold on;
    plot([dataSIO.E_rms*1000 dataSIO.E_rms*1000],[0 max(Ncount_Erms)],'-','color',col(1,:),'linewidth',3)
    set(gca,'fontsize',16,'linewidth',2); box on; grid off
    title('Misfit','fontsize',18);
    xlabel('RMS (ms)','fontsize',18);
    
    %% Ftest
    x_grid = invdata.Ftest_res.x_grid;
    y_grid = invdata.Ftest_res.y_grid;
    z_grid = invdata.Ftest_res.z_grid;
    [Xgrd,Ygrd,Zgrd] = meshgrid(x_grid,y_grid,z_grid);
    P = invdata.Ftest_res.Pstat;
    [Pz_max, Iz_max] = max(max(max(P)));
    [Py_max, Iy_max] = max(max(P(:,:,Iz_max)));
    [Px_max, Ix_max] = max(P(:,Iy_max,Iz_max));

    PLOT_Ftest_all
    f3 = gcf;
    plot(ax1,dataSIO.x_sta,dataSIO.y_sta,'p','color','yellow','markersize',25,'linewidth',2);
    plot(ax2,dataSIO.x_sta,dataSIO.z_sta,'p','color','yellow','markersize',25,'linewidth',2);
    plot(ax3,dataSIO.y_sta,dataSIO.z_sta,'p','color','yellow','markersize',25,'linewidth',2);
    col = parula;
%     % set background to match P=0 colour
%     set([ax1,ax2,ax3],'Color',col(1,:))
%     fill(ax1,[-2000,2000,2000,-2000,-2000],[-2000,-2000,2000,2000,-2000],col(1,:),'layer','bottom')
%     fill(ax2,[-2000,2000,2000,-2000,-2000],[-8000,-8000,8000,8000,-8000],col(1,:),'layer','bottom')
%     fill(ax3,[-2000,2000,2000,-2000,-2000],[-8000,-8000,8000,8000,-8000],col(1,:),'layer','bottom')
    
    % expand bounds and ticks if needed to include SIO station
    set(ax1,'xlim',[min([x_grid,dataSIO.x_sta-5]),max([x_grid,dataSIO.x_sta+5])],...
            'ylim',[min([y_grid,dataSIO.y_sta-5]),max([y_grid,dataSIO.y_sta+5])])
    set(ax2,'xlim',[min([x_grid,dataSIO.x_sta-5]),max([x_grid,dataSIO.x_sta+5])],...
            'ylim',[min([z_grid,dataSIO.z_sta-5]),max([z_grid,dataSIO.z_sta+5])])
    set(ax3,'xlim',[min([y_grid,dataSIO.y_sta-5]),max([y_grid,dataSIO.y_sta+5])],...
            'ylim',[min([z_grid,dataSIO.z_sta-5]),max([z_grid,dataSIO.z_sta+5])])
    dtick = [1,2,4,5,10,20,40,50,100];
    dtickx = dtick(mindex(abs((axlim(ax1,2)-axlim(ax1,1))/7 - dtick)));
    dticky = dtick(mindex(abs((axlim(ax1,4)-axlim(ax1,3))/7 - dtick)));
    dtickz = dtick(mindex(abs((axlim(ax2,4)-axlim(ax2,3))/7 - dtick)));
    set(ax1,'xtick',[-1000:dtickx:1000],'xticklabel',[-1000:dtickx:1000],...
            'ytick',[-1000:dticky:1000],'yticklabel',[-1000:dticky:1000]);
    set(ax2,'xtick',[-1000:dtickx:1000],'xticklabel',[-1000:dtickx:1000],...
            'ytick',[-6000:dtickz:6000],'yticklabel',[-6000:dtickz:6000]);        
    set(ax3,'xtick',[-1000:dticky:1000],'xticklabel',[-1000:dticky:1000],...
            'ytick',[-6000:dtickz:6000],'yticklabel',[-6000:dtickz:6000]);
    % kill text
    delete([htx1,hty1,htx2,htz2,hty3,htz3])
    % plot range
    x12 = [dataSIO.x_sta,invdata.loc_xyz(1)];
    y12 = [dataSIO.y_sta,invdata.loc_xyz(2)];
    z12 = [dataSIO.z_sta,invdata.loc_xyz(3)];
    %dxy
    plot(ax1,x12,y12,'-y','linewidth',2)
    text(ax1,mean(x12),mean(y12),num2str(sqrt(diff(x12).^2 + diff(y12).^2),'%.1f m'),...
        'fontsize',14,'color','yellow','fontweight','bold',...
        'rotation',atand(vertexag(ax1)*diff(y12)/diff(x12)),...
        'horizontalalignment','center','verticalalignment','top')
    %dxz
    plot(ax2,x12,z12,'-y','linewidth',2)
    text(ax2,mean(x12),mean(z12),num2str(sqrt(diff(x12).^2 + diff(z12).^2),'%.1f m'),...
        'fontsize',14,'color','yellow','fontweight','bold',...
        'rotation',atand(vertexag(ax2)*diff(z12)/diff(x12)),...
        'horizontalalignment','center','verticalalignment','top')
    %dyz
    plot(ax3,y12,z12,'-y','linewidth',2)
    text(ax3,mean(y12),mean(z12),num2str(sqrt(diff(z12).^2 + diff(y12).^2),'%.1f m'),...
        'fontsize',14,'color','yellow','fontweight','bold',...
        'rotation',-atand(vertexag(ax3)*diff(z12)/diff(y12)),...
        'horizontalalignment','center','verticalalignment','top')

    
        
    %% Save output
    if ~exist(outdir_SIOcomp)
        mkdir(outdir_SIOcomp);
    end
    
%     % Save textfile
%     fid = fopen([outdir,'/',data.sta,'_location.txt'],'w');
%     fprintf(fid,'Bootstrap inversion results (2sigma uncertainty)');
%     fprintf(fid,'\nStation: %s',data.sta);
%     fprintf(fid,'\nLat:   %.5f deg (%f) \nLon:   %.5f deg (%f) \nX:     %f m (%f) \nY:    %f m (%f) \nDepth: %f m (%f) \nTAT:   %f ms (%f) \nWater Vel.: %f m/s (%f)',mean(lat_sta),std(lat_sta)*2,mean(lon_sta),std(lon_sta)*2,mean(x_sta),std(x_sta)*2,mean(y_sta),std(y_sta)*2,mean(z_sta),std(z_sta)*2,mean(TAT)*1000,std(TAT)*1000*2,mean(v_w),std(v_w)*2);
%     fprintf(fid,'\nDrift Lon: %f m (%f) \nDrift Lat: %f m (%f) \nDrift:    %f m (%f) \nDrift Azi: %f deg (%f)\n',mean(dx_drift),std(dx_drift)*2,mean(dy_drift),std(dx_drift)*2,mean(drift),std(drift)*2,mean(azi),std(azi)*2);
%     fprintf(fid,'\nRMS:  %f ms (%f)\n',mean(E_rms)*1000,std(E_rms)*2*1000);
%     fprintf(fid,'\nBad pings Removed: %d',N_badpings);
%     fprintf(fid,'\n===================================================\n');
%     fprintf(fid,'%10s %10s %15s %15s %15s %15s \n','Lat','Lon','Range (m)','Residual (s)','Ship Vel. (m/s)','TWT corr. (ms)');
%     for i = 1:Nobs
%         fprintf(fid,'%3d: %10f %10f %10f %10f %10f %10f\n',i,lats_ship(i),lons_ship(i),range(i),dtwt_bs(i)*1000,sqrt(sum(v_ship(:,i).^2)),dtwtcorr_bs(i)*1000);
%     end
%     fclose(fid);
    
    % Save plots
    if ~exist([outdir_SIOcomp,'/plots/'])
        mkdir([outdir_SIOcomp,'/plots/']);
    end
    save2pdf([outdir_SIOcomp,'/plots/',data.sta,'_1_OBSlocations.pdf'],f6,500)
    save2pdf([outdir_SIOcomp,'/plots/',data.sta,'_2_misfits.pdf'],f7,500)
    save2pdf([outdir_SIOcomp,'/plots/',data.sta,'_3_Ftest.pdf'],f3,500)
%     save2eps(f3,[outdir_SIOcomp,'/plots/',data.sta,'_3_Ftest.eps'])
    

    if is_savemat
        datamat.par = [];
        datamat.sta = sta;
        datamat.drop_lonlatz = [dataSIO.lon_drop,dataSIO.lat_drop,dataSIO.z_drop];
        datamat.lons_ship = dataSIO.lat;
        datamat.lats_ship = dataSIO.lon;
        datamat.x_ship = [];
        datamat.y_ship = [];
        datamat.z_ship = [];
        datamat.v_ship = [];
        datamat.databad = [];
        datamat.dtwt_bs = dataSIO.res;
        datamat.twtcorr_bs = [];
        datamat.dtwtcorr_bs = [];
        datamat.lon_sta_bs = dataSIO.lon_sta;
        datamat.lat_sta_bs = dataSIO.lat_sta;
        datamat.x_sta_bs = dataSIO.x_sta;
        datamat.y_sta_bs = dataSIO.y_sta;
        datamat.z_sta_bs = dataSIO.z_sta;
        datamat.drift_bs = dataSIO.drift;
        datamat.azi_bs = dataSIO.azi;
        datamat.TAT_bs = [];
        datamat.V_w_bs = [];
        datamat.E_rms = dataSIO.E_rms;
        datamat.Ftest_res = [];
        datamat.loc_xyz = [mean(dataSIO.x_sta),mean(dataSIO.y_sta),mean(dataSIO.z_sta)];
        datamat.loc_lolaz = [mean(dataSIO.lon_sta),mean(dataSIO.lat_sta),mean(dataSIO.z_sta)];
        datamat.mean_drift_az = [mean(dataSIO.drift) r2d(mean_ang(d2r(dataSIO.azi)))];
        datamat.R_mat = [];
        datamat.Cm_mat = [];
        if ~exist([outdir_SIOcomp,'/mats'])
            mkdir([outdir_SIOcomp,'/mats']);
        end
        save([outdir_SIOcomp,'/mats/',data.sta,'_data.mat'],'datamat');
    end
end
