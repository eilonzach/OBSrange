%% Function to produce Figure 08
%  Figure 08 shows the meso scale eddy
function Figure08

ofile = '../Figure08';
ifsave = 1;

col = colormap(parula); close(gcf) % colormap for f test plots
lonlim = [-136.3 -129.9];
latlim = [-9 -3];

%% load 
load('../../data/yORCA_locations.mat');

%% process into nicer structures
Nstas = length(allstas);
Erms = zeros(Nstas,1);
for is = 1:Nstas
    Vp(is,1) = mean(allstas(is).V_w_bs);
    Erms(is,1) = mean(allstas(is).E_rms);
end
normErms = Erms./max(Erms);
drop_lolaz = reshape([allstas.drop_lonlatz],3,Nstas)';
loc_lolaz = reshape([allstas.loc_lolaz],3,Nstas)';
drift_az = reshape([allstas.mean_drift_az],2,Nstas)';


%% ---------------------   PLOTTING   ---------------------   

%% Map with drift on
f908 = figure(908); clf; hold on
set(f908,'position',[440 6 900 750]);

% % water vel surface
% [Xgrd,Ygrd] = meshgrid(linspace(lonlim(1),lonlim(2),43),linspace(latlim(1),latlim(2),45));
% F = scatteredInterpolant(loc_lolaz(:,1),loc_lolaz(:,2),Vp(:),'linear');
% Vgrd = F(Xgrd,Ygrd);
% % nan out edges
% Vgrd(sqrt((Xgrd+133.0).^2/10 + (Ygrd+6).^2/8)>1) = nan;
% % Vgrd = griddata(loc_lolaz(:,1),loc_lolaz(:,2),Vp(:),Xgrd,Ygrd);
% % plot surface
% contourf(Xgrd,Ygrd,Vgrd,[1495:1510],'linestyle','none');
% % overlay white box to damp colours
% patch(axlim(gca,[1 2 2 1 1]),axlim(gca,[3 3 4 4 3]),'w','FaceAlpha',0.7,'linestyle','none')


% plot drifts
scale = 0.25e-2;
p1s = drop_lolaz(:,1:2);
p2s = p1s + scale*(drift_az(:,1)*[1 1]).*[sind(drift_az(:,2)),cosd(drift_az(:,2))];

px = reshape([p1s(:,1),p2s(:,1),nan(Nstas,1)]',3*Nstas,1);
py = reshape([p1s(:,2),p2s(:,2),nan(Nstas,1)]',3*Nstas,1);
plot(px,py,'k-','linewidth',2)

% plot arrows to show movement
p3s = p2s + 30*scale*[sind(drift_az(:,2)+165),cosd(drift_az(:,2)+165)];
p4s = p2s + 30*scale*[sind(drift_az(:,2)+195),cosd(drift_az(:,2)+195)];
px = reshape([p3s(:,1),p2s(:,1),p4s(:,1),nan(Nstas,1)]',4*Nstas,1);
py = reshape([p3s(:,2),p2s(:,2),p4s(:,2),nan(Nstas,1)]',4*Nstas,1);
plot(px,py,'k-','linewidth',2.)
plot(px,py,'k-','linewidth',2.)

% plot drop points, size by rms, colour by depth
% scatter(loc_lolaz(:,1),loc_lolaz(:,2),80./normErms,-abs(loc_lolaz(:,3)),'filled','p','linewidth',1,'markeredgecolor','k')
scatter(loc_lolaz(:,1),loc_lolaz(:,2),80./normErms,Vp,'filled','p','linewidth',1,'markeredgecolor','k')


% colourbar
hc = colorbar; set(hc,'linewidth',2)
set(get(hc,'title'),'string','Water velocity (m/s)','fontweight','bold')

% key/scale for drift lines are, and how RMS translates to symbol size
set(gca,'xlim',lonlim,'ylim',latlim);
% key box
fill(axlim(gca,1)+1.1*[0 1 1 0 0],axlim(gca,4)-2.1*[0 0 1 1 0],'k','linewidth',2,'facecolor','w')
% RMS scales
rmss = [6,3,2,1];
for ir = 1:length(rmss)
    scatter(axlim(gca,1)+0.35,4.5+axlim(gca,3)+ir/4,80*6./rmss(ir),'filled','p','linewidth',1,'markeredgecolor','k','markerfacecolor',[0.3 0.3 0.3])
    text(axlim(gca,1)+0.55,   4.5+axlim(gca,3)+ir/4,sprintf('%u ms',rmss(ir)),'verticalalignment','middle','fontweight','bold','fontsize',14,'color',[0.3 0.3 0.3]);
end
text(axlim(gca,1)+0.55,axlim(gca,4)-0.2,'Data RMS','horizontalalignment','center','fontweight','bold','fontsize',14,'color',[0.3 0.3 0.3]);

plot(axlim(gca,1)+0.55+scale*200*[-0.5 0.5],axlim(gca,4)-1.89* [1 1],'color',[0.3 0.3 0.3],'linewidth',2)  ;
text(axlim(gca,1)+0.55,axlim(gca,4)-1.7,'200m drift','horizontalalignment','center','fontweight','bold','fontsize',14,'color',[0.3 0.3 0.3]);

% pretty axes
xlabel('\textbf{Longitude}','interpreter','latex','fontsize',18); 
ylabel('\textbf{Latitude}','interpreter','latex','fontsize',18); 
axis equal
set(gca,'fontsize',15,'linewidth',2,'box','on','xlim',lonlim,'ylim',latlim);



%% SAVE
if ifsave
    save2pdf(ofile,f908)
end

end
