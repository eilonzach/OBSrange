%
% interpolate SSPS derived from Levitus DataBase of typeSSP (as described below)
%  to get Sound Speed Profiles at (lons, lats)
%
% typeSSP = 0 indicates annual Levitus,
% typeSSP = 1:12 indicates Levitus for month 1:12
% typeSSP = 13:16 indicates Levitus for Winter, Spring, Summer, or Fall  

% typeVAR = 1 indicates sound speed
%           2 indicates temperature
%           3 indicates salinity
%           4 indicates buoyancy

% originally set up to read only sound speed variables, now modified to work
% on these other variables.

function levSSPs = getLevSSPs (typeSSP, lons, lats, dataBaseDir,typeVAR);

if typeVAR==1,
   qnq='';
elseif typeVAR==2,
   qnq='T';
elseif typeVAR==3,
   qnq='S';
elseif typeVAR==4,
   qnq='N';
end

if isunix, slsh = '/'; else, slsh = '\'; end

if typeSSP == 0,
levFileName = sprintf ('%s%slev%s_ann.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 1,
levFileName = sprintf ('%s%slev%s_jan.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 2,
levFileName = sprintf ('%s%slev%s_feb.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 3,
levFileName = sprintf ('%s%slev%s_mar.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 4,
levFileName = sprintf ('%s%slev%s_apr.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 5,
levFileName = sprintf ('%s%slev%s_may.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 6,
levFileName = sprintf ('%s%slev%s_jun.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 7,
levFileName = sprintf ('%s%slev%s_jul.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 8,
levFileName = sprintf ('%s%slev%s_aug.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 9,
levFileName = sprintf ('%s%slev%s_sep.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 10,
levFileName = sprintf ('%s%slev%s_oct.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 11,
levFileName = sprintf ('%s%slev%s_nov.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 12,
levFileName = sprintf ('%s%slev%s_dec.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 13,
levFileName = sprintf ('%s%slev%s_win.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 14,
levFileName = sprintf ('%s%slev%s_spr.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 15,
levFileName = sprintf ('%s%slev%s_sum.mat',dataBaseDir, slsh, qnq);
elseif typeSSP == 16,
levFileName = sprintf ('%s%slev%s_fal.mat',dataBaseDir, slsh, qnq);
end

eval (['load ' levFileName ])
% this loads in sound speed: (1000 + c*0.01) is sound speed
%               temperature: t/1000 is temperature
%               salinity:    s/1000 is salinity
% OR            buoyancy frequency: vf/1000 is buoyancy in cy/hour

if typeVAR==2,
   c=t;
elseif typeVAR==3,
   c=s;
elseif typeVAR==4,
   c=vf;
end

levDataName = sprintf ('%s%slev_latlonZ.mat',dataBaseDir, slsh);
eval (['load ' levDataName])
lon=reshape(lon,360,180);
lat=reshape(lat,360,180);
lat=lat*0.1;
lon=lon*0.1;
% also loads in z, the levitus depths.
% this is a function so z stays in this subroutine.

% interpolate one depth at a time.
bb=[];

[n m]=size(lons);
if m>n,
lons=lons';
lats=lats';
end

for i=1:33,
c1=c(:,i);
c1=reshape(c1,360,180);
c1=reshape(c1,360,180);
bb=[bb interp2(lon',lat',c1',lons,lats,'*cubic')];
end

% convert bb to real sound speed values
% rows of bb are the sound speed profiles
if typeVAR==1,
   bb=1000 + bb*0.01;
else
   bb=bb*0.001;
end

levSSPs = cell (size (lons));
for k=1:length(lons),
   levSSPs{k} = num2cell(bb(k,:));
end
