#!/usr/bin/env python2.7

from matplotlib import pyplot as plt
import sys

ifn1 = "XBT562.txt"
ifn2 = "XBT563.txt"

# load data
InFile1 = open(ifn1)
data1 = []
for line in InFile1:
    data1.append(line.split())
InFile1.close()

InFile2 = open(ifn2)
data2 = []
for line in InFile2:
    data2.append(line.split())
InFile2.close()

Depth = []
TDiff = []
T1 = []
T2 = []
VDiff = []
for i in range(len(data1)):
    # check to make sure depth is same
    ZDiff = float(data1[i][0]) - float(data2[i][0])
    if ZDiff != 0:
        print "Depths not the same! Diff is "+str(ZDiff)+"  Exiting..."
        sys.exit()
    Depth.append(float(data1[i][0])*-1)
    TDiff.append(float(data1[i][1]) - float(data2[i][1]))
    T1.append(float(data1[i][1]))
    T2.append(float(data2[i][1]))
    VDiff.append(float(data1[i][2]) - float(data2[i][2]))

Plot = plt.figure()
plt.plot(TDiff,Depth)
plt.title("Difference in XBT Temperature Measurements")
plt.ylabel("Depth (m)")
plt.xlabel("Temp Diff. (degrees C)")
#plt.show()
Plot.savefig("XBT_TDiff.pdf")

Plot = plt.figure()
plt.plot(T1,Depth)
plt.title("XBT Temperature Measurement 1")
plt.ylabel("Depth (m)")
plt.xlabel("Temp. (degrees C)")
#plt.show()
Plot.savefig("XBT_T1.pdf")

Plot = plt.figure()
plt.plot(T2,Depth)
plt.title("XBT Temperature Measurement 2")
plt.ylabel("Depth (m)")
plt.xlabel("Temp. (degrees C)")
#plt.show()
Plot.savefig("XBT_T2.pdf")

Plot = plt.figure()
plt.plot(VDiff,Depth)
plt.title("Difference in XBT Velocity Estimates")
plt.ylabel("Depth (m)")
plt.xlabel("Velocity Diff. (m/s)")
#plt.show()
Plot.savefig("XBT_VDiff.pdf")


