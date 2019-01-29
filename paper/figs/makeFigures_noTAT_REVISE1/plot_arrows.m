function [ ] = plot_arrows( ax, X, Y )

% PLOT ARROWS
Xdata = [X X];
Ydata = [Y-15 Y];
%get axes drawing area position in normalized units
axpos = get(ax, 'Position');
%get axes drawing area in data units
ax_xlim = xlim(ax);
ax_ylim = ylim(ax);
ax_norm_per_xdata = axpos(3) ./ diff(ax_xlim);
ax_norm_per_ydata = axpos(4) ./ diff(ax_ylim);
%these are figure-relative
Xnorm = (Xdata - ax_xlim(1)) .* ax_norm_per_xdata + axpos(1);
Ynorm = (Ydata - ax_ylim(1)) .* ax_norm_per_ydata + axpos(2);
annotation('arrow',Xnorm,Ynorm,'LineStyle','-','color',[0 0 0],'linewidth',2);

end

