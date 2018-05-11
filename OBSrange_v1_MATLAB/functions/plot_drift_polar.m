function [ mean_drift_az ] = plot_drift_polar( ax,driftaz )
%[ mean_drift_az ] = plot_drift_polar( ax,driftaz )
%   Plot polar histogram of drift direction

% polar histogram
polarhistogram(ax,driftaz,36,'Normalization','pdf','FaceColor',[0 0.1 0.8],'linewidth',0.1);

% mean drift direction
mean_drift_az = mean_ang(driftaz);
polarplot(ax,mean_ang(driftaz)*[1 1],get(ax,'rlim'),'linewidth',2,'color',[0.8500 0.3250 0.0980]);

% bounding circle
polarplot(ax,d2r([0:1:360]'),ones(361,1)*max(get(ax,'rlim')),'k','linewidth',3);

% station
polarplot(ax,0,0,'^','linewidth',2,'markersize',12,'markerfacecolor',[ 0.0980 0.9500 0.4250],'markeredgecolor','k');

set(ax,'ThetaZeroLocation','top','ThetaDir','clockwise',...
    'thetatick',[0:90:270],'thetaticklabel',{'N','E','S','W'},...
    'rticklabel',{},'fontweight','bold','fontsize',14);




end

