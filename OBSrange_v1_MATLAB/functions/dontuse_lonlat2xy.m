function [ x, y ] = lonlat2xy( olon, olat, lon, lat )
% [ x, y ] = lonlat2xy( olon, olat, lon, lat )
% Convert lon/lat (deg) to x/y (m) using the reference ellipsoid

% e = referenceEllipsoid('GRS80');
% [x, y, ~] = geodetic2enu(lat, lon, zeros(size(lon)), olat, olon, 0, e);

% % If no mapping toolbox, use this version
% x = -(olon - lon)*cos(olat*pi/180)*111.1 * 1000; % m
% y = -(olat - lat)*111.1 * 1000; % ms

az_rot = 0;
proj.origin  = [olat olon];   % model origin
mstruct = defaultm('mercator');
mstruct.origin = [proj.origin -az_rot];
mstruct = defaultm( mstruct ); 
proj.map_proj = mstruct;
[x, y] = dontuse_project_xy(proj, lat, lon, 'forward');
x = x*1000;
y = y*1000;

end

