#!/usr/bin/env python
# coding: utf-8

# Very basic script loading data from LabNSN2 data files and plotting a histogram

import sys

if len(sys.argv)<2:
    print('Usage: python LabNSN2DataExample.py file1 file2 ... fileN')
    sys.exit()
    
datafiles = sys.argv[1:]

import csv, numpy
import matplotlib.pyplot as plt

plt.rc('font', family = 'serif')
plt.rcParams["figure.figsize"] = (15,8)
    

# (fake) tdc calibrations, so that time = tdc_value*slope_tdc + offset_tdc
slope_tdc = 4.000			# tdc bin size, ns
offset_tdc = 64			# tdc time offset, ns

# Read the files, append to list only if at least one plane has a valid TDC stop
data = []
for filen in datafiles:
    with open(filen, newline = '') as f:
        reader = csv.reader(f, delimiter = ' ') 
        for row in reader:
            row = row[1:]      #skip the '_'
            p = [int(row[j],16) for j in [2,3,4]]    # I don't care about event# and time here
            if any((item != 4095) for item in p):    # I don't care about empty events here
                data.append(p)
# now data contains the three planes stop values only for non-empty events

# Select data - here only events with P3 TDC stop and nothing else
p3 = []
for ev in data:
    if ev[0] == 4095 and ev[1] == 4095:
        p3.append(ev[2]*slope_tdc+offset_tdc)
        
# Plot it and format
fig, axs = plt.subplots()
y3,x3,h3 = axs.hist(x = p3, color = 'tab:blue', bins = numpy.arange(offset_tdc, 4095*slope_tdc + offset_tdc, 32*slope_tdc))
axs.set_yscale('log')
axs.title.set_fontsize(28)
axs.title.set_text('Type BOTTOM (P3) - (' + str(len(p3)) + ' events)')
axs.set_ylabel('Entries', fontsize = 24)
axs.set_xlabel('Time difference (ns)', fontsize = 24)
plt.show()