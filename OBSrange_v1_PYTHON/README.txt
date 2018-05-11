To Do List from V1:

1) I tried to implement the eigenvector scaling for the Ftest, but something
   isn't working quite right. Can't figure it out at this point. However, the 
   original version of the ftest seems to work fine.

2) Still haven't added output functionality, but that'll be simple enough.

3) I changed the structure of the code from the previous python verions I made,
   but I introduced an error somewhere. Many of the residuals are now too high,
   notably the misfit RMS. I can't find the error at this point.

4) I made my own smoothing function for velocity smoothing since I'm not sure
   what the MATLAB curve fitting toolbox does.

5) The plots have been tidied up, but there's a number of annoyances I can't
   work out for the time being:
     - the tick label formats for the lat and lon plots are terrible.
     - plotting specific contours on the Ftest error plots.

6) Check drift azimuth.

Notes from V2:

1) I fixed the errors I had introduced into V1 be re-arranging my code. :) Versions give comparable results now. RMSs of 1.5 ms for EE04!!

2) I updated the smoothing function in python by translating Zach's moving average code in matlab.

3) The only other "big thing" to do now is re-work the ftest and eigenvector scaling. But I anticipate this will be smooth now that the other aforementioned errors have been worked out.

4) Still need to add output functionality. Will be straightforward but tedious.

5) Still need to update plots, but they're workable for now.

6) I added a row (or column, whatever it was) to the resampled bootstrap matrices which corresponds to unscrambeld raw data. Makes comparisons much easier! 

