function [result1,result2] = project_xy( par, input1, input2, direction )
% project_xy( par, input1, input2, direction )
%
% PROJECT_XY
%
% projects station lat/lon to cartesian
% and back again
%
% if direction=='forward'
%         input1 = lat
%         input2 = lon
%         result1 = x
%         result2 = y
% if direction=='inverse'
%         input1 = x
%         input2 = y
%         result1 = lat
%         result2 = lon
%

% error( nargchk(3,4,nargin) );
narginchk(3,4)

if ( nargin == 3 )
	direction = 'forward';
elseif ( ~ischar(direction) || (~strcmpi(direction,'forward') && ...
		~strcmpi(direction,'inverse')) )
	error( 'incorrect 4th argument' );
end


switch ( lower(direction) )

	% lat,lon -> x,y
	case 'forward'

		lat = input1;
		lon = input2;
		
		if ( isnaive(par) )
			y = deg2km(lat - par.origin(1));
			x = deg2km(lon-par.origin(2)) * cos( deg2rad(par.origin(1)) );
		elseif ( isflat(par) )
			y = deg2km(lat - par.origin(1));
			x = deg2km(lon - par.origin(2));
		elseif ( ismaptoolbox(par) )
			[x,y] = mfwdtran( par.map_proj, lat, lon );
			x = rad2km( x, 6371 );
			y = rad2km( y, 6371 );
		else
			error( 'invalid choice for map_proj' );
		end

		result1 = x;
		result2 = y;
		
	% x,y -> lat,lon
	case 'inverse'

		x = input1;
		y = input2;
		
		if ( isnaive(par) )
			
			error( 'haven''t implemented naive inverse projection' );
			
		elseif ( ismaptoolbox(par) )
			
			x = km2rad( x, 6371 );
			y = km2rad( y, 6371 );
			[lat,lon] = minvtran( par.map_proj, x, y );
			
		end
		
		result1 = lat;
		result2 = lon;
		
	otherwise

		error( 'direction must be forwards or backwards' );

end





%
% ISFLAT
%
function truefalse = isflat( par )

truefalse = ~isstruct(par.map_proj) && strcmpi('flat',par.map_proj);



%
% ISNAIVE
%
function truefalse = isnaive( par )

truefalse = ~isstruct(par.map_proj) && strcmpi('naive',par.map_proj);


%
% ISMAPTOOLBOX
%
function truefalse = ismaptoolbox( par )

truefalse = isstruct(par.map_proj) && isfield(par.map_proj,'mapprojection');


