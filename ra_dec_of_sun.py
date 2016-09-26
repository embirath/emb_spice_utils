#!/usr/bin/env python2.7

import spiceypy as spice
import numpy as np

# Input: Time (in s/c cblock format)
met = '3/299182112.000'

# Specify target
target = "Sun"

# Load kernels
spice.furnsh('/home/emma/kernels/nh/science_spice/15sci_rhr.tm')

# Convert to ET time
et = spice.scs2e(-98,met)

# Find RA/Dec of Sun. First use spice.spkpos() to get the rectangular
# vector, then use spice.recrad() to get Radious, RA, Dec
target_vector_rec,ltime = spice.spkpos(target,et, 'J2000', 'NONE', 'NEW_HORIZONS')
target_vector_rad = spice.recrad(target_vector_rec)

# Report the result, convert Ra/Dec to degrees
print "Target: " + target
print "RA: {}".format(target_vector_rad[1]*spice.dpr())
print "Dec: {}".format(target_vector_rad[2]*spice.dpr())




