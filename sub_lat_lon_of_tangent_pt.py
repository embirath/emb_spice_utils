#!/usr/bin/env python2.7

# This program reports the Pluto sub-latitude/longitude of 
# the tangent point of a vector "nh-to-star". Input is time
# (in s/c clock format), and Ra/Dec of the star (in degrees).
# The program outputs the Pluto latitude/ longtidue 
# values in degrees. 
#   
#                  
#                         . (star)                
#                        /|\                          
#                         |                               
#                         |                      * * *                
#                         |                    *       *             
#                   .     |/__________________*___ . PlutoCenter                
#                  /|\    |\  pluto_to-tpoint *  /|     *           
#                   |     |                    */      *           
#                   |     |                    / * * *            
#                   |     |                   /                
#         dot-vector|     |                  /                
#                   |     |                 /                
#                   |     |                / nh_to_pluto               
#                   |     |               /                
#                   |     |              /                
#                   |     |             /                   
#                   |     |            /                   
#                   |     |           /                   
#                   |     |          /      
#                   |     |         /      
#                   |     |        /       
#                   |     |       /        
#                   |     |      /         
#                   |     |     /          
# nh_to_star_unit   |     |    /           
#             .     |     |   /            
#            /|\    |     |  / 
#             |     |     | /  
#             |     |     |/   
#                         NH               
#                              

import spiceypy as spice
import numpy as np

# Input Time (in s/c clock format) and the RA/Dec in question
met = '3/299195097.000'
ra = 93.85469568        # deg
dec = 16.14320591       # deg
print ""
print "Input"
print "   MET: " + met
print "   RA/Dec: " + "{0:.2f} / {1:.2f} ".format(ra,dec)


# Load kernels
spice.furnsh('/home/emma/kernels/nh/science_spice/15sci_rhr.tm')

# Convert MET time string to ET
et = spice.scs2e(-98,met)

# Find the NH to Pluto vector, call it nh_to_pluto
nh_to_pluto, ltime = spice.spkpos('Pluto',et, 'J2000', 'NONE', 'NEW_HORIZONS')

# Create a unit vector that points to the specified RA/Dec, 
# call it nh_to_star_unit
nh_to_star_unit = spice.radrec(1.0, ra*spice.rpd(), dec*spice.rpd())

# Now, find vector from Pluto center to Tangent Point on the nh_to_star_unit
# vector. Call it pluto_to_tpoint.
# 
#    We know that: 
#    nh_to_pluto + pluto_to_tpoint = dot_vector 
#    where 
#    dot_vector = (nh_to_pluto DOT nh_to_star_unit ) * nh_to_star_unit
#   
#    Or re-arranged:
#    pluto_to_tpoint = (nh_to_pluto DOT nh_to_star_unit ) * nh_to_star_unit - nh_to_pluto
#    
#    We will use two SPICE functions: 
#    dot_product = spice.vdot( vec1, vec2 )
#    scalar_x_vector = spice.vscl(scalar, unit_vec)
#    
#    or combined: 
#    spice.vscl(spice.vdot(vec1,vec2),unit_vec)
dot_vector = spice.vscl(spice.vdot(nh_to_pluto,nh_to_star_unit),nh_to_star_unit)
pluto_to_tpoint = dot_vector - nh_to_pluto

# Convert the rectangular coordinates to Radius, RA, and Dec
radius, ra_tpoint, dec_tpoint = spice.recrad(pluto_to_tpoint)

# Now we want to rotate this vector from J2000 into Pluto coordinates
# 
# Get the transformation matrix from "J2000" to "IAU_PLUTO"
#   rot_matrix = spice.pxform('from','to',et)
#
# And then rotate pluto_to_tpoint from "J2000" to "IAU_PLUTO"
#   pluto_to_tpoint_in_pluto_coords = numpy.dot(rot_matrix, pluto_to_tpoint)
#
# Note there is a a bug in Python, which gives a warning when you run PXFORM. 
# "ctypes currently produces invalid PEP 3118 type codes, 
# which Numpy notices: http://bugs.python.org/issue10746 http://bugs.python.org/issue10744"
# To hide the warning, we use the warnings module. 
#
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    rot_matrix = spice.pxform("J2000","IAU_PLUTO",et)
pluto_to_tpoint_in_pluto_coords = np.dot(rot_matrix, pluto_to_tpoint)

radius, lat_tpoint, lon_tpoint = spice.recrad(pluto_to_tpoint_in_pluto_coords)
print ""
print "Results: "
print "   Radius: " + "{0:.2f}".format(radius) + " km"
print "   Latitude: " + "{0:.2f}".format(lat_tpoint*spice.dpr()) + " deg"
print "   Longitude: " + "{0:.2f}".format(lon_tpoint*spice.dpr()) + " deg"

print ""

