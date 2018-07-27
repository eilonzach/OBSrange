clear; close all;

%% Parameters
% projpath = '/Users/zeilon/Documents/MATLAB/OBS_range_dist/projects/Worcester_OBS-Relocations'; % path to project
projpath = '/Users/zeilon/Work/OBSrange/projects/PacificORCA/'; % path to project
outdir_OBSrange = './OUT_OBSrange/OUT_nocorr/'; % output directory

projname =fliplr(strtok(fliplr(projpath),'/'));
wd = pwd;
cd(projpath);cd(outdir_OBSrange);
files = dir('mats');
Nfiles = length(files);
Nstas = 0;
for ifile = 1:Nfiles
    if any(regexp(files(ifile).name,'mat'))
        Nstas = Nstas+1;
        load(['mats/',files(ifile).name]);
        if Nstas==1 
            allstas = datamat;
        else
            allstas(Nstas,1) = datamat;
        end
    end
end

        
%% water speed map
Erms = zeros(Nstas,1);
r_sta_stds = zeros(Nstas,1);
for is = 1:Nstas
    Vp(is) = mean(allstas(is).V_w_bs);
    Erms(is) = mean(allstas(is).E_rms);
    F(is) = mean(allstas(is).E_rms);
end
drop_lolaz = reshape([allstas.drop_lonlatz],3,Nstas)';
loc_lolaz = reshape([allstas.loc_lolaz],3,Nstas)';
mean_drift_az = reshape([allstas.mean_drift_az],2,Nstas)';
fprintf('Multibeam vs. relocated (+ means relocated shallower\n')
abs(drop_lolaz(:,3))-abs(loc_lolaz(:,3));
normErms = Erms./max(Erms);


%% plot water velocity
figure(1);clf
scatter(loc_lolaz(:,1),loc_lolaz(:,2),2e5*Erms,Vp,'filled')
hc = colorbar;

set(gca,'fontsize',15,'linewidth',2,'box','on')

%% plot drift
figure(2), set(gcf,'pos',[440 172 911 626]);clf, hold on
scale = 0.24e-2;
for is = 1:Nstas
%     plot(allstas(is).drop_lonlatz(1),allstas(is).drop_lonlatz(2),'.r','markersize',30)
%     plot(allstas(is).loc_lolaz(1),allstas(is).loc_lolaz(2),'ro','markersize',10)
%     plot(allstas(is).drop_lonlatz(1) + scale*allstas(is).mean_drift_az(1)*[0 1 0.9 1 0.9].*sind(allstas(is).mean_drift_az(2)+ [0 0 2 0 -2]),...
%          allstas(is).drop_lonlatz(2) + scale*allstas(is).mean_drift_az(1)*[0 1 0.9 1 0.9].*cosd(allstas(is).mean_drift_az(2)+ [0 0 2 0 -2]),...
%          'k','linewidth',2)     
    plot(allstas(is).drop_lonlatz(1) + scale*allstas(is).mean_drift_az(1)*[0 1].*sind(allstas(is).mean_drift_az(2)+ [0 0]),...
         allstas(is).drop_lonlatz(2) + scale*allstas(is).mean_drift_az(1)*[0 1].*cosd(allstas(is).mean_drift_az(2)+ [0 0]),...
         'k','linewidth',2)     
%     quiver(allstas(is).drop_lonlatz(1),allstas(is).drop_lonlatz(2),...
%         scale*allstas(is).mean_drift_az(1)*sind(allstas(is).mean_drift_az(2)),...
%         scale*allstas(is).mean_drift_az(1)*cosd(allstas(is).mean_drift_az(2)))
    scatter(loc_lolaz(:,1),loc_lolaz(:,2),80./normErms,-abs(loc_lolaz(:,3)),'filled','marker','p','markeredgecolor','k','linewidth',1)
end
plot(-131.5+scale*200*[0 1],-3.5* [1 1],'color',[0.3 0.3 0.3],'linewidth',2)  ;
text(-131.5+scale*100,-3.35,'200m drift','horizontalalignment','center','fontweight','bold','fontsize',14,'color',[0.3 0.3 0.3]);
set(gca,'fontsize',15,'linewidth',2,'box','on')
hc = colorbar;
set(get(hc,'title'),'string','Water depth (m)','fontweight','bold')
xlabel('Longitude'); ylabel('Latitude')


[abs(loc_lolaz(:,3)),sqrt(mean_drift_az(:,1).^2 + loc_lolaz(:,3).^2)]

save2pdf(2,'ORCAdrift.pdf');

%% find mean predicted error in phase velocity
Npair = handshake(Nstas);
truL = nan(Npair,1); dropL = nan(Npair,1);
ipair = 0;
for is1 = 1:Nstas
for is2 = (is1+1):Nstas
    ipair = ipair+1;
    truL(ipair) = distance_km_nomap(loc_lolaz(is1,2),loc_lolaz(is1,1),loc_lolaz(is2,2),loc_lolaz(is2,1));
    dropL(ipair) = distance_km_nomap(drop_lolaz(is1,2),drop_lolaz(is1,1),drop_lolaz(is2,2),drop_lolaz(is2,1));
end
end
gdpair = truL<200;
dL_L = 100*abs(truL(gdpair)-dropL(gdpair))./truL(gdpair);
dc = dL_L*4.0;
fprintf('Mean error in phase velocity of %.2f km/s\n',mean(dc))

return
%% print to csv file:
fid = fopen([projname,'_allstas.csv'],'w');
fprintf(fid,'#, STA,  LAT(deg.),  LON(deg.),  depth(m), driftX, driftY, TAT, V_water, Erms  \n');
for is = 1:length(deployarray.slats)
    if deployarray.selevs(is)==0, depstr = ' '; else, depstr = sprintf('%4.0f',-1000*deployarray.selevs(is)); end
    fprintf(fid,'%2.0f, %s, %2.0f%s%06.3f%sS, %4.0f%s%06.3f%sW, %9.5f, %9.5f,   %s,    %6.1f,      %6.1f \n',...
        deployarray.order(is),deployarray.stas{is},...
        abs(fix(deployarray.slats(is))),'  ',60*abs(rem(deployarray.slats(is),1)),39,...
        abs(fix(deployarray.slons(is))),'  ',60*abs(rem(deployarray.slons(is),1)),39,...
        [deployarray.slats(is),deployarray.slons(is)],...
        depstr,deployarray.nm_from_last(is),deployarray.nm_to_next(is));
end
fclose(fid);


