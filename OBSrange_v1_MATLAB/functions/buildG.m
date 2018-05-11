function [ G ] = buildG( x, y, z, dvp, x_ship, y_ship, z_ship, vp_w, Nobs, M)
% Build the G matrix for 2-way travel time inversion

G = zeros(Nobs,M);

% Shortcut for distance...
D = sqrt((x_ship-x).^2 + (y_ship-y).^2 + (z_ship-z).^2);

% Setup G Matrix  
% dti/dx
G(:,1) = -(x_ship-x) .* 2 ./ (vp_w+dvp) ./ D;
% dti/dy
G(:,2) = -(y_ship-y) .* 2 ./ (vp_w+dvp) ./ D;
% dti/dz
G(:,3) = -(z_ship-z) .* 2 ./ (vp_w+dvp) ./ D;
% dti/dTAT
G(:,4) = ones(Nobs,1);
% dti/ddvp
G(:,5) = -2 .* D ./ (vp_w+dvp).^2;

end

