#!/usr/bin/env python
##############################################################################
# Copyright (c) 2004 U.C. San Diego
#      
# Filename:    obslocate
# Version:     0.2
# Date:        Tue Jan 18 11:50:02 2005
# Author:      Paul Georgief (Based on Tammy Katherine Baldwin's Matlab code)
# Email:       paul@ucsd.edu, tkb@indiana.edu
# Purpose:     Using a grid-based least squares fit, find the location of an 
#              OBS on the bottom of the ocean floor
#
##############################################################################
import math, sys, re

##############################################################################
# Name          : printLogo
# Purpose       : Display the 'obslocate' logo.
# Inputs        : none
# Outputs       : none
# Bidirectional : none
# Returns       : none
# Notes         : -
##############################################################################
def printLogo():
    "Print program logo."
    print "\n" * 3, "=" * 60
    print "       _         _                 _       "
    print "  ___ | |__  ___| | ___   ___ __ _| |_ ___ "
    print " / _ \\| '_ \\/ __| |/ _ \\ / __/ _` | __/ _ \\"
    print "| (_) | |_) \\__ \\ | (_) | (_| (_| | ||  __/"
    print " \\___/|_.__/|___/_|\\___/ \\___\\__,_|\\__\\___|"
    print "=" * 60


##############################################################################
# Name          : askForDropPoint
# Purpose       : Obtain initial drop point information
# Inputs        : none
# Outputs       : none
# Bidirectional : none
# Returns       : List containing (correctionX, correctionY, depth, (lat, lon))
# Notes         : none
##############################################################################
def askForDropPoint():
    """Obtain initial drop point information.
    This function returns a list: (correctionX, correctionY, depth, (lat, lon))
    usage: a = askForDropPoint()
    """
    lat   = float(raw_input("Drop point Latitude (decimal format  [N:+, S:-]): "))
    lon   = float(raw_input("Drop point Longitude (decimal format [E:+, W:-]): "))
    depth = float(raw_input("Drop point depth (meters): "))
    return (0, 0, depth, (lat, lon))

##############################################################################
# Name          : isDataFrom8011M
# Purpose       : Is the data file generated by 8011M or an old deck box?
# Inputs        : none
# Outputs       : none
# Bidirectional : none
# Returns       : True-8011M generated, False-Old Style
# Notes         : none
##############################################################################
def isDataFrom8011M():
    """Ask if the data was generated using an 8011M.
    This function returns: False-Old Type, True-New Type (8011M).
    The old acoustic deck boxes returned the range value as a time (in milliseconds)
    that was the one way travel time.  The new 8011M boxes return the 2 way travel time
    (including release turn around time).  Therefore, to get the one way travel time
    we must apply the following formula:

    oneWayTime = (Time - TurnaroundTime)/2.0

    The turn around time for the acoustic releases is given as 13mSec.
    usage: a = isDataFrom8011M()
    """

    text = raw_input("Was an 8011M acoustic box used to collect the data ([Y]/N): ")
    if text == '' or text[0] in ['Y', 'y']:
        return True
    return False


##############################################################################
# Name          : readDataFile
# Purpose       : To parse a text file and extract lat, lon, and range data
# Inputs        : Filename    - The file to open
#                 device8011M - True=8011M used, False=Old System 
# Outputs       : none
# Bidirectional : none
# Returns       : data     - A list of (lat, lon, travelTimeRange, string)
# Notes         : -
#   
##############################################################################
def readDataFile(filename, device8011M=True):
    """Read a text file containing ranging information and place the data into a list.
    The format of the data will be:

    Rng                Position                      Altitude    Sample Time
    ===  =====================================       ==========  ============================
    304 msec. Lat: 32 38.2945 N  Lon: 117 26.0574 W  Alt: -25.10 Time(UTC): 2004:329:00:58:58
    267 msec. Lat: 32 38.2812 N  Lon: 117 26.0229 W  Alt: -24.79 Time(UTC): 2004:329:00:59:18
    230 msec. Lat: 32 38.2684 N  Lon: 117 25.9877 W  Alt: -25.39 Time(UTC): 2004:329:00:59:38

    Lines not fitting the format will be ignored.
    This function returns a dataList containing: [lat, lon, tTime, line] where
        lat   - Ship's Latitude
        lon   - Ship's Longitude
        tTime - Ranging travel Time (msec).  Multiply by velocity of sound in water to get distance.
        line  - The line of text corresponding to the data point.  This is done so that when a
                modified list is created, it can be save to a file without the need for back conversions 
                (and loss of data since this routine ignores the 'Sample Time').
                
    usage: data = readDataFile(filename)
    """

    dataFormat1 = re.compile (r"""
                ^                       # Go to the beginning of the line
                \s*                     # Grab any leading spaces (or none)
                (\d+)                   # Look for the range value (in msec)
                \s+msec\.\s+Lat:\s+     # One or more space(s) "msec.", space(s), "Lat:", space(s)
                (\d+)\s+(\S+)\s+        # Degrees, space(s), minutes (any character except whitespace), space(s)
                (N|S)                   # North or South marker
                \s+Lon:\s+              # Space(s), "Lon:", follwed by 1 or more spaces
                (\d+)\s+(\S+)\s+        # Degrees, space(s), minutes (any character except whitespace), space(s)
                (W|E)                   # West or East marker
                \s+Alt:\s+              # One or more spaces, "Alt:", followed by 1 ore more spaces
                (\S+)                   # Altitude (any character except whitespace)
                .*                      # Grab the rest of the line
                """,
                re.IGNORECASE | 
                re.VERBOSE        # !!!! NEED !!!! Allows commented string (like above)
                )


    # Read all the lines in a file and extract lines that match template
    lines = open(filename, 'r').readlines()
    dataList = []
    for line in lines:
        match1 = dataFormat1.match(line)
        if (match1):
            # It matches the current data format - compute the values and add it to the list
            d = match1.groups()
            if device8011M:
                tTime = (float(d[0])-13.0)/2.0  #  Turn around Time = 13mSec; Two way travel time
            else:
                tTime = float(d[0])
            lat = (float(d[1]) + float(d[2])/60.0) * [-1,1][d[3] in ['n', 'N']]
            lon = (float(d[4]) + float(d[5])/60.0) * [-1,1][d[6] in ['e', 'E']]
            alt = float(d[7])
            dataList.append([lat, lon, tTime, line.strip()])

    return dataList


##############################################################################
# Name          : convertLatLonToXYDeltas
# Purpose       : To flatten an Lat/Lon coordinate to a xy-location
# Inputs        : deltaLat    - latitude difference between fixed point and
#                               point of interest.
#                 deltaLon    - longitude difference between fixed point and
#                               point of interest.
#                 lat         - latitude of fixed point
# Outputs       : none
# Bidirectional : none
# Returns       : Differential (X, Y) cooridinates corresponding to the 
#                 differential (dlon, dlat) coordinates.
# Notes         : This uses a spherical representation of the earth with
#                 the radius of the earth (at the equator) equal to 111Km.
#
##############################################################################
def convertLatLonToXYDeltas(deltaLat, deltaLon, fixedLat):
    """Convert LAT/LON-deltas to xy-delts (measured in meters).
    This treats the Earth as a sphere - longitude distance is dependent 
    on the cos(latitude)

    Radius of earth at equator: 6378136m  (Merit, 1983)
    1 degree = 60 minutes = (6378136 * 2 * pi / 360) = 111.31947 Km

    The conversion is done on coorinates relative to a fixed
    point.  The earth model breaks if you use the crossing of 
    the equator and meridian line as the origin and should be avoided.

    Only the real (absolute) latitude of the fixed position 
    is necessary to compute the xy-deltas.
    usage:  (deltaX, deltaY) = convertLatLonToXYDeltas(deltaLat, deltaLon, fixedLat)
    """
    degToMeter = 111319.47

    dy = deltaLat * degToMeter
    dx = deltaLon * degToMeter * math.cos(fixedLat * (math.pi / 180.))
    
    return (dx, dy)
    
##############################################################################
# Name          : convertXYToLatLonDeltas
# Purpose       : Take the (dx,dy) data and compute (dLat,dLon) offsets.
# Inputs        : dx       - x positional change from fixed location
#                 dy       - y positional change from fixed location
#                 fixedLat - Latitude of fixed location.
# Outputs       : none
# Bidirectional : none
# Returns       : deltaLat - Latitude difference corresponding to dy.
#                 deltaLon - longitude difference corresponding to dx.
# Notes         : The spherical model of the earth breaks if you use absolute
#                 positions.  Therefore, this routine converts x/y relative
#                 distances to lon/lat relative differences.
##############################################################################
def convertXYToLatLonDeltas(dx, dy, fixedLat):
    """Convert XY-deltas (measured in meters) to Lat/Lon-deltas.
    usage:  (deltaLat, deltaLon) = convertXYToLatLonDeltas(dx, dy, fixedLat)
    """
    degToMeter = 111319.47

    deltaLat   = dy / degToMeter
    deltaLon   = dx / degToMeter / math.cos(fixedLat * (math.pi / 180.))
    return (deltaLat, deltaLon)

##############################################################################
# Name          : convertDeltaXYToFixedLatLon
# Purpose       : To take dx & dy distances in meters and used a fixed
#                 (lat, lon) to determine a new (lat, lon)
# Inputs        : none
# Outputs       : none
# Bidirectional : none
# Returns       : lat, lon  - Latitude, Longitude
# Notes         : -
##############################################################################
def convertDeltaXYToFixedLatLon(dx, dy, fixedLat, fixedLon):
    """Convert xy-deltas to an absolute location.
    usage: (lat, lon) = convertDeltaXYtoFixedLatLon(dx, dy, fixedLat, fixedLon)
    """
    (dLat, dLon) = convertXYToLatLonDeltas(dx, dy, fixedLat)
    lat, lon  = (fixedLat + dLat, fixedLon + dLon)
    return (lat, lon)

##############################################################################
# Name          : dec2deg
# Purpose       : To change decimal degrees to degrees and degree Minutes.
# Inputs        : dec       - decimal latitude or longitude
# Outputs       : none
# Bidirectional : none
# Returns       : deg, min  - degrees and minutes
# Notes         : -
##############################################################################
def dec2deg(dec):
    "Convert decimal degrees to degree minutes."
    deg = int(dec)
    min = (dec - float(deg)) * 60.0
    return (deg, math.fabs(min))
    

##############################################################################
# Name          : computeResidual
# Purpose       : Compute the error between the measured range and the computed
#                 range.  For this problem the computed range is calculated using
#                 the pythagorean theorem (distance between 2 objects in x,y,z space).
# Inputs        : r, x, y, z  - Measured range, x, y, depth
# Outputs       : none
# Bidirectional : none
# Returns       : residual
# Notes         : -
##############################################################################
def computeResidual (r, x, y, z):
    "Compute residual between measured range and calculated range."
    return (r - math.sqrt((x*x) + (y*y) + (z*z)))

##############################################################################
# Name          : findMatrixMin
# Purpose       : Search through a 2D array and determine where in the matrix
#                 the minimum value is located.
# Inputs        : grid    -  2D array
# Outputs       : none
# Bidirectional : none
# Returns       : row, column, and minimum value
# Notes         : There must be a quicker/nicer way to do this?
##############################################################################
def findMatrixMin(grid):
    "Locate the position in a 2-D array where the minimum value occurs."
    row,col = (len(grid), len(grid[0]))
    (minr, minc, minVal) = (0, 0, grid[0][0])

    for r in range(row):
        for c in range(col):
            if minVal > grid[r][c]:
                (minr, minc, minVal) = (r, c, grid[r][c])
    return (minr, minc, minVal)


##############################################################################
# Name          : computeGridResiduals
# Purpose       : (see notes below)
# Inputs        : gridSize      - Determines how many points we will search
#                                 (A 2D array of gridSize x gridSize will be created)
#                 spacing       - Spacing (in meters) between the grid points
#                 dropInfo      - The grid will be centered around this point.
#                 dataList      - Collected data points.
# Outputs       : none
# Bidirectional : none
# Returns       : New drop point
# Notes         : -
##############################################################################
def computeGridResiduals(gridSize, spacing, dropInfo, dataList):
    """Search a 2-D area looking for the smallest error between measured and 
    calculated errors.
    
    A 2D grid is created with given size and spacing.  The grid is centered around an
    initial OBS drop position.  A sweep is done "moving" the OBS drop position to
    a new "assumed" drop position (a point on the 2D grid).  

    Using the measured boat positions and a fixed depth; the range between "assumed" OBS
    position and the boat can be calculated.  This calculated range is subtracted from
    the actual measured range to give an error (residual).  This is done for every
    measured position of the ship/boat.  All the individual errors are squared (giving
    more weigth to larger errors and also removing ambiguity between +/- errors) and then
    summed.  This value, sum(error^2), is used as an indicator of the true position of
    the OBS.  It is assumed (correctly?) that the OBS position having the least sum(error^2)
    corresponds to the actual OBS position.

                                                                        xy-Plane
                                                                     (OBS Position)
       (Side View: xy-plane vs. z)                             +-----+-----+-----+-----+
                                                               |     |     |     |     |
             +--------B  <-- Boat (measured Lat/Lon)           |     |     |     |     |
             |       /                                         +-----+-----+-----+-----+
             |      /                                          |     |     |     |     |
       fixed |     /                                           |     |     |     |     |
       depth |    / measured                                   +-----+-----I-----+-----+
             |   /  range                                      |     |     |     |     |
             |  /                                              |     |     |     |     |
             | /                                               +-----A-----+-----+-----+
             |/                                                |     |     |     |     |
             O  <---- OBS on Ocean Floor (grid point)          |     |     |     |     |
                                                               +-----+-----+-----+-----+

                                        I - Assumed initial drop point  (located at 0,0) 
                                        A - Actual drop point (located at -1,-1)          
                                        The least sum(error^2) will occur at position 'A'.
    usage: ((dx, dy, depth), min) = computeGridResiduals(gridSize, spacing, dropInfo, dataList)
    """

    # Make the gridSize an odd integer
    halfLen = int(gridSize)/2
    gridLen = (halfLen*2) + 1
    (dropX, dropY, depth, fixedPoint) = dropInfo[:]

    # Create a grid centered on the drop point.  Compute the residual (error)
    # assuming a point on the grid is the "new" drop point.  
    grid = []
    obsX = dropX - (spacing * halfLen)
    obsY = dropY + (spacing * halfLen)

    for y in range(halfLen, -halfLen-1, -1):     #  halfLen, ..., 0, ..., -halfLen
        row = []
        for x in range(-halfLen, halfLen+1):     # -halfLen, ..., 0, ...,  halfLen
            sum = 0
            for (boatX, boatY, boatR, text) in dataList:
                residual = computeResidual(boatR, (boatX-obsX), (boatY-obsY), depth)
                sum += (residual * residual)
            row.append(sum)
            obsX += spacing
        # Increment the OBS y location and reset the x location
        grid.append(row)
        obsY -= spacing
        obsX  = dropX - (spacing * halfLen)

    # Find the minimum value (error^2) in the 2-d array.  Its location in the array
    # corresponds to where the real drop point is located.
    (mrow, mcol, mVal) = findMatrixMin(grid)
    (minX, minY) = (mcol-halfLen, halfLen-mrow)
    #print "Minimum Points at: (X: %d, Y:%d, Val:%.5f)" % (minX, minY, mVal)
    return ((dropX + (minX*spacing), dropY + (minY*spacing), depth, fixedPoint), mVal)


##############################################################################
# Name          : printInfo
# Purpose       : Display the Lat/Lon/Range/Depth on the screen
# Inputs        : d  - data point (xLon, yLat, zRange)
#                 printRange - Print range tag or depth tag
# Outputs       : none
# Bidirectional : none
# Returns       : none
# Notes         : -
##############################################################################
def printInfo (d,fwobj, printRange=False):
    "Display the lat, lon, depth (or range) on the screen."
    (dx, dy, depth, (fixedLat, fixedLon)) = d
    lat, lon = convertDeltaXYToFixedLatLon(dx, dy, fixedLat, fixedLon)
    latDeg, latMin = dec2deg(lat)
    lonDeg, lonMin = dec2deg(lon)
    if printRange:
        asdfprint = " Lat: %d %.4f (%.4f),  Lon: %d %.4f (%.4f), range: %.4f" % (latDeg, latMin, lat,lonDeg, lonMin, lon, d[2]) 
        print asdfprint
        fwobj.write(asdfprint)
    else:
        asdfprint = " Lat: %d %.4f (%.4f),  Lon: %d %.4f (%.4f), depth: %.4f" % (latDeg, latMin, lat,lonDeg, lonMin, lon, d[2])
        print asdfprint
        fwobj.write(asdfprint)

##############################################################################
# Name          : processData
# Purpose       : The main routine
# Inputs        : dropInfo    - Initial drop position (latitude, longitude, RangeTravelTime)
#                 data        - A list of all boat locations plus travel time ranges to the OBS
# Outputs       : none
# Bidirectional : none
# Returns       : none
# Notes         : -
##############################################################################
def processData (dropInfo, data,fwobj):
    "Cycle through several grids/spacings to find the location of the OBS."
    print   "=" * 60, "\n"
    fwobj.write("\n")
    fwobj.write("=" * 60)
    fwobj.write("\n")
    
    # Create a grid of possible drop points and compute the "error" between the assumed drop
    # point and the measured ranges/locations.  Find the location where the minimum error occurs
    # and relocate the drop point to that location.  Rerun the calculation using a tighter grid and the new
    # computed drop point.
    gridInfo = ((100, 100), (100, 10), (20, 1), (20, .1))
    newDropInfo = dropInfo
    numDataPoints = len(data)
    for (size, spacing) in gridInfo:
        asdfprint= "\nRunning calculation on %dx%d grid (with %.4f meter spacing)" % (size, size, spacing)
        print asdfprint
        fwobj.write(asdfprint)
        asdfprint= "\nDrop point located at:  \n "
        print asdfprint
        printInfo(newDropInfo,fwobj)
        asdfprint= ""
        print asdfprint
        fwobj.write(asdfprint)

        (newDropInfo, error) = computeGridResiduals(gridSize=size, spacing=spacing, dropInfo=newDropInfo, dataList=data)
        asdfprint =  "\nNew drop point location: \n"
        print asdfprint
        fwobj.write(asdfprint)
        printInfo(newDropInfo,fwobj)
        asdfprint = "\nsum(Residual^2) = %.4f\n" % error
        print asdfprint
        fwobj.write(asdfprint)
        if (numDataPoints > 0):
            asdfprint = "\nsqrt(sum(Residual^2)/%d) = %.4f\n" % (numDataPoints, (math.sqrt(error/float(numDataPoints))))
            print asdfprint
            fwobj.write(asdfprint)
        asdfprint = "\n"
        fwobj.write(asdfprint)

    # ------------------------------------------------------------
    # Print out the individual residuals 
    # ------------------------------------------------------------
    print "\n\n"
    fwobj.write("\n\n")
    asdfprint= "Individual residuals (based at final drop point): "
    print asdfprint
    fwobj.write(asdfprint)
    count = 1
    for (boatX, boatY, boatRange, text) in data:
        asdfprint = "\n%3d) " %count
        fwobj.write(asdfprint)
        print "%3d) " % count, 
        residual = computeResidual(boatRange, (boatX-newDropInfo[0]), (boatY-newDropInfo[1]), newDropInfo[2])
        printInfo((boatX, boatY, boatRange, newDropInfo[3]),fwobj, True)
        print "Individual residual: ", residual
        fwobj.write("Individual residual: ")
        fwobj.write("%f"%residual)
        count += 1


    # ------------------------------------------------------------
    # Print out the initial/final locations and offsets
    # ------------------------------------------------------------
    fwobj.write("\n\nInitial Drop:    ")
    print "\nInitial Drop:  ",
    printInfo(dropInfo,fwobj)

    fwobj.write("\nFinal Drop:    ")
    print "\nFinal Drop:    ",
    printInfo(newDropInfo,fwobj)
    fwobj.write('\n')
    print "\n",
    
    fwobj.write(" " * 16)
    fwobj.write("Number of points (N):")
    fwobj.write("%f\n"%numDataPoints) 
    print " " * 16, "Number of points (N):", numDataPoints
    fwobj.write(" " * 16)
    fwobj.write("sum(Residual^2) = %.4f\n"%error)
    print " " * 16, "sum(Residual^2) = %.4f" % error
    if (numDataPoints > 0):
        print  " " * 16, "sqrt(sum(Residual^2)/N)  : ", math.sqrt(error/float(numDataPoints))
        fwobj.write("\n")
        fwobj.write(" " * 16)
        fwobj.write("sqrt(sum(Residual^2/N)  : ")
        fwobj.write("%f\n"%math.sqrt(error/float(numDataPoints)))


    (xlon, ylat) = (newDropInfo[0] - dropInfo[0], newDropInfo[1] - dropInfo[1])
    asdfprint = "\nOffset Distance: Lat=%.4f meters, Lon=%.4f meters, (r=%.4f meters, angle=%.2f)\n" % (
        ylat, xlon, math.sqrt(ylat*ylat + xlon*xlon), math.atan2(ylat, xlon) * 180 / math.pi
        )
    print asdfprint
    fwobj.write(asdfprint)
    

    print "=" * 60, "\n" * 3
    fwobj.write("\n")
    fwobj.write("=" * 60)
    fwobj.write("\n"*3)



##############################################################################
# Name          : main
# Purpose       : main
# Inputs        : -
# Outputs       : -
# Bidirectional : -
# Returns       : -
# Notes         : Main Routine
##############################################################################
def main():
    "main function called only if routine is run as a script."
    printLogo()


    # ------------------------------------------------------------
    # Get the required data and preprocess the data.  Turn
    # all fixed Lat/Lon locations into relative positions with the
    # obs drop point as the center.
    # ------------------------------------------------------------
    try:
        deviceUsed8011M = isDataFrom8011M()
        soundSpeed = float(raw_input("Speed of sound in water (m/s): "))
        dataFileName = raw_input ("Enter location of file to process: ").strip()
        rawData = readDataFile(dataFileName, deviceUsed8011M)
        dropInfo = askForDropPoint()
    except:
        # Don't return, don't fix anything, don't ask for reentry... Take the lazy way 
        # out and just exit!
        print "\nIncorrect Entry.  Program exiting!\n"
        sys.exit(-1)
    outputFileName = dataFileName.split('.')
    outputFileName = outputFileName[0] + '_Corrected.txt'
    print "\n Output File Name: %s"%outputFileName
    fwobj = open(outputFileName,'w')


    (obsDepth, (obsLat, obsLon)) = (dropInfo[2], dropInfo[3])
    diffData = []
    for dLat, dLon, travelTime, text in rawData:
        (dx,dy) = convertLatLonToXYDeltas(dLat-obsLat, dLon-obsLon, obsLat)
        rng = travelTime * soundSpeed / 1000.0
        diffData.append((dx, dy, rng, text))


    # ------------------------------------------------------------
    # Compute the Residual information on the first data set
    # ------------------------------------------------------------
    residuals = [computeResidual(bRange, bX, bY, obsDepth) for (bX,bY,bRange,bText) in diffData]
    resSumSq = 0
    for i in residuals: resSumSq += (i * i)
    for i in range(len(diffData)):
        fwobj.write("\n%3d:: Boat Location: " % (i+1))
        print "%3d) Boat Location: " % (i+1),
        printInfo((diffData[i][0], diffData[i][1], diffData[i][2], (obsLat, obsLon)),fwobj, True)
        fwobj.write("Residual: ")
        fwobj.write("%f" %residuals[i])
        print "Residual: ", residuals[i]
    print "Sum (Residual^2): ", resSumSq
    fwobj.write("\nSum (Residual^2): ")
    fwobj.write("%f\n"%resSumSq)
    if len(residuals) > 0:
        print "Sqrt(Sum(Residuals^2)/N): ", math.sqrt(resSumSq/float(len(residuals)))
        fwobj.write("Sqrt(Sum(Residuals^2)/N): ")
        fwobj.write("%f\n" % math.sqrt(resSumSq/float(len(residuals))))
    # ------------------------------------------------------------

    newDataSet = []
    for i in range(len(residuals)):
        newDataSet.append(diffData[i])

    # Print out the list of data points using original data lines
    if 0:
        fwobj.write("#" * 70)
        print "#" * 70
        fwobj.write("Number of Data Points:    ")
        fwobj.write(len(newDataSet))
        print "Number of Data Points:   ", len(newDataSet)
        fwobj.write("Drop Point Latitude:    ")
        fwobj.write("obsLat")
        print "Drop Point Latitude:     ", obsLat
        fwobj.write("Drop Point Longitude    ")
        fwobj.write(obslon)
        print "Drop Point Longitude:    ", obsLon
        fwobj.write("Drop Point Depth (m)    ")
        fwobj.write(obsDepth)
        print "Drop Point Depth (m):    ", obsDepth
        fwobj.write("Velocity of Water (m/s)    ")
        fwobj.write(soundSpeed)
        print "Velocity of Water (m/s): ", soundSpeed
        for i in newDataSet:
            print i[3].strip()
            fwobj.write(i[3].strip())
        print "#" * 70
        fwobj.write("#" * 70)
   

    # Process the data and see if we wish to recompute with the new setup
    processData(dropInfo, newDataSet,fwobj)

    
# ----------------------------------------------------------------------------------
# Check for script invocation (run 'main' if the script is not imported as a module)
# ----------------------------------------------------------------------------------
if __name__ ==  '__main__':
    # Get the data and convert the lat/lon into kilometers and travel time into range (in meters)
    main()
