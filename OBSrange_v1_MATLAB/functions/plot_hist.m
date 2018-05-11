function [ Ncount, cent ] = plot_hist( ax,x,Nbins)
% [ Ncount, cent ] = plot_hist( ax,x,Nbins )
% Plot histograms

x_mean = mean(x);
x_std = std(x);

hold(ax,'on')

[Ncount, cent] = hist(x,Nbins);
Ncount_norm = Ncount./sum(Ncount);
hb = bar(ax,cent,Ncount_norm,1,'facecolor',[0 0.1 0.9]);

% plot([x_mean, x_mean],[0 max(Ncount)],'-','color',[0.8500 0.3250 0.0980],'linewidth',2);
% plot([x_mean+1*x_std, x_mean+1*x_std],[0 max(Ncount)],'--','color',[0.8500 0.3250 0.0980],'linewidth',3);
% plot([x_mean-1*x_std, x_mean-1*x_std],[0 max(Ncount)],'--','color',[0.8500 0.3250 0.0980],'linewidth',3);

x_95_50 = prctile(x,[2.5 50 97.5]);

plot(ax,[x_95_50(1) x_95_50(1)],[0 1.2*max(Ncount_norm)],'--','color',[0.8500 0.3250 0.0980],'linewidth',3);
plot(ax,[x_95_50(3) x_95_50(3)],[0 1.2*max(Ncount_norm)],'--','color',[0.8500 0.3250 0.0980],'linewidth',3);
plot(ax,[x_95_50(2) x_95_50(2)],[0 1.2*max(Ncount_norm)],'-','color',[0.8500 0.3250 0.0980],'linewidth',3);

set(ax,'fontsize',16,'linewidth',2,'Xgrid','on','box','on',...
    'ylim',[0 1.05*max(Ncount_norm)]);


end

