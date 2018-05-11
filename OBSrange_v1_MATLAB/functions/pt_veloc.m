function [ v ] = pt_veloc( x, y, z, t )
% Calculate velocity between points in 3-D
%
% Josh Russell : 4/18/18
%

% Check for stationary positions
dt = diff(t);
Idt0 = (dt==0);

% Calculate velocity
v_half = [diff(x)./dt, diff(y)./dt, diff(z)./dt]';
v_half(:,Idt0') = 0; % replace stationary positions with zero velocity;
v = (v_half(:,1:end-1) + v_half(:,2:end))/2;
v = [v_half(:,1), v, v_half(:,end)]; % Replace endpoints

end

