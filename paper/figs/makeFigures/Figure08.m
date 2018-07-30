%% Function to produce Figure 08
%  Figure 08 shows the meso scale eddy
function Figure08

ofile = '../Figure08';
ifsave = 0;

col = colormap(parula); % colormap for f test plots

%% load 
load('../figdata/junk_allstas_wcorr.mat');

%% process into nicer structures
Nstas = length(allstas);
Erms = zeros(Nstas,1);
for is = 1:Nstas
    Vp(is) = mean(allstas(is).V_w_bs);
    Erms(is) = mean(allstas(is).E_rms);
end
normErms = Erms./max(Erms);
drop_lolaz = reshape([allstas.drop_lonlatz],3,Nstas)';
loc_lolaz = reshape([allstas.loc_lolaz],3,Nstas)';
drift_az = reshape([allstas.mean_drift_az],2,Nstas)';


%% ---------------------   PLOTTING   ---------------------   

%% Map with drift on
f908 = figure(908); clf; hold on
set(f908,'position',[440 6 900 750]);

% plot drifts
scale = 0.26e-2;
p1s = drop_lolaz(:,1:2);
p2s = p1s + scale*(drift_az(:,1)*[1 1]).*[sind(drift_az(:,2)),cosd(drift_az(:,2))];

px = reshape([p1s(:,1),p2s(:,1),nan(Nstas,1)]',3*Nstas,1);
py = reshape([p1s(:,2),p2s(:,2),nan(Nstas,1)]',3*Nstas,1);
plot(px,py,'k-','linewidth',2)

% plot drop points, size by rms, colour by depth
scatter(loc_lolaz(:,1),loc_lolaz(:,2),80./normErms,-abs(loc_lolaz(:,3)),'filled','p','linewidth',1,'markeredgecolor','k')

% plot arrows to show movement
% p3s = p1s + 0.92*scale*(drift_az(:,1)*[1 1]).*[sind(drift_az(:,2)+3),cosd(drift_az(:,2)+3)];
% p4s = p1s + 0.92*scale*(drift_az(:,1)*[1 1]).*[sind(drift_az(:,2)-3),cosd(drift_az(:,2)-3)];
% px = reshape([p3s(:,1),p2s(:,1),p4s(:,1),nan(Nstas,1)]',4*Nstas,1);
% py = reshape([p3s(:,2),p2s(:,2),p4s(:,2),nan(Nstas,1)]',4*Nstas,1);
% plot(px,py,'k-','linewidth',2)
% plot(px,py,'k-','linewidth',2)

% pretty axes
set(gca,'fontsize',15,'linewidth',2,'box','on','ylim',[-9 -3])
xlabel('Longitude'); ylabel('Latitude')
axis equal

% colourbar
hc = colorbar; set(hc,'linewidth',2)
set(get(hc,'title'),'string','Water depth (m)','fontweight','bold')

% key/scale
plot(-131.5+scale*200*[0 1],-3.5* [1 1],'color',[0.3 0.3 0.3],'linewidth',2)  ;
text(-131.5+scale*100,-3.35,'200m drift','horizontalalignment','center','fontweight','bold','fontsize',14,'color',[0.3 0.3 0.3]);


%% SAVE
if ifsave
    save2pdf(ofile,f908)
end

end
