function [rx,rz,Dr,rt,rv] = shootrays_jbr(p, v_profile, zmax, dr,vdz)
% [rx,rz,Dr] = ray_trace_1D(p, v_profile, zmax, dr,vdz)
% 
% 1D ray tracing function
% all units in km or km/s
if nargin < 4
    dr = 0.001;
end
if nargin<5
    vdz = 0.005;
end

zz1 = [0:vdz:(zmax-vdz)]';
zz2 = [vdz:vdz:zmax]';

v1  = linterp(v_profile.z,v_profile.vwater,zz1);
v2  = linterp(v_profile.z,v_profile.vwater,zz2);
% assume flat Earth (no radius terms)
u1 = 1./v1;
u2 = 1./v2; 

dv = v2-v1;

dz = zz2-zz1;
b = dv ./ dz;
const_indx = (b == 0);

X = zeros( length(v1), 1 );
X(const_indx) = constv_dist( u1(const_indx), dz(const_indx), p );
X(~const_indx) = gradv_dist(b(~const_indx),u1(~const_indx),u2(~const_indx),p);
X = [0; X];
X(imag(X)>0) = NaN;
X = X(~isnan(X));
Xd = cumsum(X);
rayz = [zz1(1); zz2]; 
rayz = rayz(1:length(X));
dz = dz(1:length(X)-1);

% calc Dr
Dr1 = (X(2:end).^2 + dz.^2).^(1/2);
if any(~isreal(Dr1)), error; end
Dr1 = [0; Dr1];
Dr2 = cumsum(Dr1);
Dr = [0:dr:max(Dr2)]';
Dr = [Dr; max(Dr2)];
% cDr = Dr;
rx = interp1(Dr2,Xd,Dr);
rz = interp1(Dr2,rayz,Dr);
% nn = 1:round(length(Dr)/length(rz)):length(Dr)+10;
% nn = nn(1:length(rz));
% Dr = interp1(1:length(Dr),Dr,nn)'; % this will be exactly the same Dr as above
% cDr = interp1(1:length(cDr),cDr,nn)'; % this is the same as Dr

% calculate travel time (JBR)
dDr = Dr(2:length(Dr))-Dr(1:length(Dr)-1);
dDr = [0; dDr];
v = [v1(1); v2];
rv = interp1(Dr2,v,Dr);
dt = dDr./rv;
rt = cumsum(dt);

%-------------------- sub-functions --------------------------------------
%
function x = constv_dist( u, dz, p )

x = (p * dz) ./ eta(u,p);


function x = gradv_dist( b, u1, u2, p )

x1 = eta(u1,p) ./ (b.*u1*p);
x2 = eta(u2,p) ./ (b.*u2*p);

x = x1 - x2;

%
function val = eta( u, p )

val = sqrt( u.^2 - p.^2 );
