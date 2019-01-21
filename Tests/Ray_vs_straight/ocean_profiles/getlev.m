% function [P,z]=getlev(glat,glon,'param',typessd,ifofile)
%
% A routine to extract SSPs, temperature, salinity, or buoyancy profiles 
% along a given point or points along a path (glat,glon)
%
% Uses "dist.m", which uses WGS84 ellipsoid by default, and
% "getLevSSPs.m", a program to load the variables and interpolate
% onto the desired location(s).
%
% Gets the World Ocean Atlas files from directory '/home/dushaw/mapprogram/dataBase'
% Modify this program to set the directory where files lev_ann.mat, 
% lev_latlonZ.mat, levN_ann.mat, levS_ann.mat, and/or levT_ann.mat are located.
% These files are preprocessed to save disk space, and have all shortened
% profiles of the original data filled in by nearest neighbor data.
% All profiles extend to 5500 m depth, and there is world wide coverage
% (yes, including continents)
%
%   See:  http://staff.washington.edu/dushaw/WOA/
%   Data originally from:  http://www.nodc.noaa.gov/OC5/indprod.html
%
% Data at points where the original WOA has data should be unchanged from the WOA.
% The annual mean World Ocean Atlas is used.
%
%
% glat,glon are a single point, or points along a path, perhaps found by dist.m
% 'param' is one of 'c', 'T', 'S', or 'N' for
% sound speed, temperature, salinity, or buoyancy, respectively.
%
% Returns P, a matrix of profiles: 33X(No. of points)
% and z, the standard 33 World Ocean Atlas depths
%
% Interpolation from World Ocean Atlas grid points to 
% desired points by cubic spline interpolation horizontally.
%
% Writes out file templev.ssp; if sound speed is called, 
% this can be used in the RAY or EIGENRAY programs.
%
% The format of file templev.ssp is two columns: 
% first,            "-1  range(m)', 
% followed by   "depth(n) variable(n,m)" for n=1:33
% and repeated for points glat(m), glon(m)
%
% UNITS:  sound speed in m/s, temperature in degrees, salinity in ppt,
% buoyancy in cy/hr, depth in positive meters.

function [P,z]=getlev(glat,glon,param,typeSSP,ifofile)

if nargin < 5
    ifofile = false;
end

dataBaseDir='./allmats09';
%
%A sample line for Windows: 
%dataBaseDir='c:\Documents and Settings\portermic.SAIC-US-WEST\Desktop\SSP database';

if glon(1) < 0,
   glon=glon+360;
end

stdDpts = [0 10 20 30 50 75 100 125 150 200 250 300 400 500 600 700 800 900 ...
           1000 1100 1200 1300 1400 1500 1750 2000 2500 3000 3500 4000 4500 ...
           5000 5500]; 

typeStr = ['Annual   '; 'January  '; 'February '; 'March    '; 'April    '; ...
                        'May      '; 'June     '; 'July     '; 'August   '; ...
                        'September'; 'October  '; 'November '; 'December '; ...
                        'Winter   '; 'Spring   '; 'Summer   '; 'Fall     '];

   % typeSSP = 0 indicates annual Levitus,
   % typeSSP = 1:12 indicates Levitus for month 1:12
   % typeSSP = 13:16 indicates Levitus for Winter, Spring, Summer, or Fall
%typeSSP=0;
%types 1-16 not implemented yet, if ever.

   % typeVAR = 1 indicates sound speed
   %           2 indicates temperature
   %           3 indicates salinity
   %           4 indicates buoyancy
if param=='c'
    typeVAR=1;
elseif param=='T'
    typeVAR=2;
elseif param=='S'
    typeVAR=3;
elseif param=='N'
    typeVAR=4;
else
disp('Variable type must be one of ''c'', ''T'', ''S'', or ''N''')
return
end

% obtains ssps (or T, S, N) from a levitus database
predSSPs = getLevSSPs (typeSSP, glon, glat, dataBaseDir,typeVAR);

lats=glat(1);
lons=glon(1);
rtmp=[];

for k=1:length(glon)
       rtmp = [rtmp dist( [lats glat(k)], [lons glon(k)] ) ];
end

% r is returned in meters.
rtmp = rtmp*0.001;


%  first load variables into var
var = zeros (33, length(rtmp));

for k = 1:length(rtmp)
   var1 = predSSPs{k}; 
   var2 = cat (1, var1{:});
   var(1:length(var2),k) = var2;
end

% Finally print SSPs
if ifofile
    disp ('Sound speed, or other variable, along path saved in file templev.ssp.')
    fid = fopen ('templev.ssp', 'w');
    for k = 1:length(rtmp)
       fprintf (fid, '-1 %11.3f\n', rtmp(k));
       for l=1:33
          if typeVAR==1
            fprintf (fid, '%5d %10.2f\n', stdDpts(l), var(l,k));
          else
            fprintf (fid, '%5d %10.3f\n', stdDpts(l), var(l,k));
          end
       end
    end
    fclose (fid);
end


P=var(:);
z=stdDpts(:);

