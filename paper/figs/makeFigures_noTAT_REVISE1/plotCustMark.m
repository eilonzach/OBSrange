%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS SIMPLE CODE ENABLES YOU TO PLOT USING CUSTOM MARKERS.             %
% AS THE INPUT, YOU NEED TO PROVIDE X AND Y FOR THE DATA POINTS, LIKE    %
% YOU DO IN REGULAR PLOTS. YOU ALSO NEED TO PROVIDE X AND Y DATA FOR     %
% YOUR CUSTOM MARKER, ASSUMING THE CENTER OF THE MARKER IS AT (0,0).     %
% YOU NEED TO ALSO PROVDE THE SCALE FOR YOUR MARKER, WHICH IS ANALOGOUS  %
% TO MATLAB'S MARKER SIZE.                                               %
% SEE THE EXAMPLE FOR CLARIFICATION! :D                                  % 
% PREPARED BY: SALMAN MASHAYEKH                                          %
% DATE: DECEMBER 2012                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [patchHndl,lineHndl] = plotCustMark(ax,xData,yData,markerDataX,markerDataY,markerSize)


%% make sure the inputs are OK
if ~isvector(xData) || ~isvector(yData) || length(xData)~=length(yData)
    fprintf('Error! xData and yData mar must be vectors of the same length!\n');
    return
end
if ~isvector(markerDataX) || ~isvector(markerDataY) || length(markerDataX)~=length(markerDataY)
    fprintf('Error! markerDataX and markerDataY mar must be vectors of the same length!\n');
    return
end
% ------
xData = reshape(xData,length(xData),1) ;
yData = reshape(yData,length(yData),1) ;
markerDataX = markerSize * reshape(markerDataX,1,length(markerDataX)) ;
markerDataY = markerSize * reshape(markerDataY,1,length(markerDataY)) ;
% -------------------------------------------------------------


%% prepare and plot the patches
markerEdgeColor = [0 0 0] ;
markerFaceColor = [0 1 1] ;
lineStyle = '--' ;
lineColor = [1 0 0] ;
% ------
vertX = repmat(markerDataX,length(xData),1) ; vertX = vertX(:) ;
vertY = repmat(markerDataY,length(yData),1) ; vertY = vertY(:) ;
% ------
vertX = repmat(xData,length(markerDataX),1) + vertX ;
vertY = repmat(yData,length(markerDataY),1) + vertY ;
% ------
faces = 0:length(xData):length(xData)*(length(markerDataY)-1) ;
faces = repmat(faces,length(xData),1) ;
faces = repmat((1:length(xData))',1,length(markerDataY)) + faces ;
% ------
% figHndl = figure ; box on ; 
lineHndl = plot(ax,xData,yData,'lineStyle',lineStyle,'Color',lineColor) ;
patchHndl = patch(ax,'Faces',faces,'Vertices',[vertX vertY]);
set(patchHndl,'FaceColor',markerFaceColor,'EdgeColor',markerEdgeColor) ;
% -------------------------------------------------------------

