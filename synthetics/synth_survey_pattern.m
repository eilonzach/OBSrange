function [ fsurvx,fsurvy ] = synth_survey_pattern( survey, radius,dl )
% [ fsurvx,fsurvy ] = synth_survey_patter( survey, rad,dl )

% all internal units in km
% radii in Nm

if nargin < 2 || isempty(radius)
    radius = 1.3;
end
if nargin < 3 || isempty(dl)
    dl = 0.005;
end

nm2km = 1.852;

switch survey

    case 'PACMAN'
    %% PACMAN
    % move-out to circle
    moy = [0:dl:radius*nm2km];
    mox = zeros(size(moy));
    %way round circle 
    npercirc = round(0.8333*2*pi*radius*nm2km./dl);
    thetas = -linspace(0,2*pi*5/6,npercirc);
    cix = radius*nm2km*sin(thetas);
    ciy = radius*nm2km*cos(thetas);
    %way back in to circle 
    mix = -[-radius*nm2km:dl:-0.5*radius*nm2km]*sin(thetas(end));
    miy = -[-radius*nm2km:dl:-0.5*radius*nm2km]*cos(thetas(end));
    % full survey path
    fsurvx = [mox,cix,mix]'; % in km
    fsurvy = [moy,ciy,miy]'; % in km
    case 'circle'
    %% CIRCLE
    % round circle 
    npercirc = round(2*pi*radius*nm2km./dl);
    thetas = -linspace(0,2*pi,npercirc);
    fsurvx = radius*nm2km*sin(thetas); % in km
    fsurvy = radius*nm2km*cos(thetas); % in km

    case 'diamond'
    %% diamond
    % move-out to diamond
    moy = [0:dl:radius*nm2km];
    mox = zeros(size(moy));
    % 3 sides
    nperside = radius*nm2km*sqrt(2)./dl; 
    dix = radius*nm2km*[linspace(0,1,nperside),linspace(1,0,nperside),linspace(0,-1,nperside)];
    diy = radius*nm2km*[linspace(1,0,nperside),linspace(0,-1,nperside),linspace(-1,0,nperside)];
    %way back in to cenbter 
    mix = [-radius*nm2km:dl:-0.5*radius*nm2km];
    miy = zeros(size(mix));
    % full survey path
    fsurvx = [mox,dix,mix]'; % in km
    fsurvy = [moy,diy,miy]'; % in km   

    case 'tri_edge'
    %% triangle with the instrument at the long edge
    % move-out to apex
    moy = [0:dl:radius*nm2km];
    mox = zeros(size(moy));
    % 2 sides
    nperside = radius*nm2km*sqrt(2)./dl; 
    dix = radius*nm2km*[linspace(0,1,nperside),linspace(1,0,nperside)];
    diy = radius*nm2km*[linspace(1,0,nperside),linspace(0,-1,nperside)];
    %way back in to station 
    miy = [-radius*nm2km:dl:-0.25*radius*nm2km];
    mix = zeros(size(miy));
    % full survey path
    fsurvx = [mox,dix,mix]'; % in km
    fsurvy = [moy,diy,miy]'; % in km   

    case 'tri_center'
    %% triangle with the instrument at the center
    vertices = radius*nm2km*[sind([330 90 210]'),cosd([330 90 210]')]; %[x,y]
    nperside = radius*nm2km/dl;
    % move-out to apex at heading 330
    mox = linspace(0,vertices(1,1),nperside);
    moy = linspace(0,vertices(1,2),nperside);
    % 2 sides
    dix = [linspace(vertices(1,1),vertices(2,1),nperside),linspace(vertices(2,1),vertices(3,1),nperside)];
    diy = [linspace(vertices(1,2),vertices(2,2),nperside),linspace(vertices(2,2),vertices(3,2),nperside)];
    %way back in to station 
    mix = linspace(vertices(3,1),0.5*radius*sind(210),0.5*nperside);
    miy = linspace(vertices(3,2),0.5*radius*cosd(210),0.5*nperside);
    % full survey path
    fsurvx = [mox,dix,mix]'; % in km
    fsurvy = [moy,diy,miy]'; % in km  
 
    case 'cross'
    %% cross, with the instrument at the center
    % move-out to point at heading 0
    moy = [0:dl:radius*nm2km];
    mox = zeros(size(moy));
    % long axis, across, long axis sides
    nperside = radius*nm2km*sqrt(2)./dl; 
    dix = radius*nm2km*[linspace(0,0,2*nperside),linspace(0,1,nperside),linspace(1,-1,2*nperside)];
    diy = radius*nm2km*[linspace(1,-1,2*nperside),linspace(-1,0,nperside),linspace(0,0,2*nperside)];
    % full survey path
    fsurvx = [mox,dix]'; % in km
    fsurvy = [moy,diy]'; % in km  
    case 'cross2'
    %% cross, with the instrument at the center, but doing two diagonals
    % move-out to point at heading 0
    moy = [0:dl:radius*nm2km];
    mox = zeros(size(moy));
    % long axis, across, long axis sides
    nperside = radius*nm2km*sqrt(2)./dl; 
    dix = radius*nm2km*[linspace(0,1,nperside),linspace(1,-1,2*nperside),linspace(-1,0,nperside)];
    diy = radius*nm2km*[linspace(1,0,nperside),linspace(0,0,2*nperside),linspace(0,-1,nperside)];
    % full survey path
    fsurvx = [mox,dix]'; % in km
    fsurvy = [moy,diy]'; % in km  
    case 'line'
    %% line, with the instrument at the center
    fsurvx = [-radius*nm2km:dl:radius*nm2km];
    fsurvy = zeros(size(fsurvx));
    case 'line2'
    %% line offset from the instrument 
    % move-out to point at heading 0
    fsurvx = [-radius*nm2km:dl:radius*nm2km];
    fsurvy = 0.3*radius*nm2km*ones(size(fsurvx));
end

end

