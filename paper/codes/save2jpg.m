function [filepath] = save2jpg(fignumber,filename,odir,renderer)
% function [filepath] = save2jpg(fignumber,filename,odir,renderer)
% 
% function to save a figure to a good-quality jpg
%
% INPUTS:
% fignumber - (required) number of figure to save
% filename - (optional, default is figN) name of saved file
% odir - (optional, default is current dir) dir in which to save file
% renderer (optional, default is -painters)
% 
% ZCE  7/2016


if nargin < 3 || isempty(odir)
    odir = pwd;
end
if nargin < 2 || isempty(filename)
    filename = sprintf('fig%u',fignumber);
end
if nargin < 4 || isempty(renderer)
    renderer = '-painters';
end

% thisdir = pwd;
% odir = cd(odir);
% cd(thisdir)

if strcmp(odir(end),'/')~=1
    odir = strcat(odir,'/');
end

if ~strcmp(filename(end-3:end),'.jpg')
    filename = [filename '.jpg'];
end

h = figure(fignumber);

%% Figure size sorting
% scrn_wd = 11.28*2.54; % inches to cm
% scrn_ht =  7.05*2.54; % inches to cm
% scrn = get(0,'ScreenSize');
pos = get(h,'position');
switch get(h,'Units')
    case 'pixels'
        ppcm = get(0,'ScreenPixelsPerInch')/2.54;
        wd = pos(3)/ppcm; % in cm
        ht = pos(4)/ppcm; % in cm
    case 'centimeters'
        wd = pos(3); % in cm
        ht = pos(4); % in cm
end
 
% if wd > ht
%     orientation = 'landscape';
% else
%     orientation = 'portrait';
% end
% set(h,'paperorientation',orientation);


set(h,'PaperUnits','centimeters');
set(h,'papersize',[wd ht]);
set(h,'PaperPosition', [0 0 wd ht]);

print(h,renderer,'-djpeg80','-r300',strcat(odir,filename))

filepath = strcat(odir,filename);
end