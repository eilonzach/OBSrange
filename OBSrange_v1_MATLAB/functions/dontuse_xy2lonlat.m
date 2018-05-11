function [ lon, lat ] = xy2lonlat( olon, olat, x, y )
% [ lon, lat ] = xy2lonlat( olon, olat, x, y )
% Convert x/y (m) to lat/lon (degrees) using the reference ellipsoid

e = referenceEllipsoid('GRS80');
[lat, lon, ~] = enu2geodetic(x, y, zeros(size(x)), olat, olon, 0, e);

end

