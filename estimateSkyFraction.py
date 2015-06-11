##
# Given the estimation of the sky magnitude at different bands at Mauna Kea
# and the galaxy brightness at the given radius in the DEIMOS band, estimate the 
# fraction between stellar light and sky brightness

from Deimos_SKiMS_slit__def__ import *

###########################
## System-wide variables ##
###########################

__builtins__.verbosity = False	#If True, keeps the verbosity level high

__builtins__.gal_Reff = 47.86		#Effective radius of galaxy in arcsec

##########
## Main ##
##########

# Distance in R_eff at which the ratio is requested
dist_sel = 3. 

# Extracting measured profile (Noordermeer+08)

realProfilePath = 'photprofs_Rband.txt'

inputData = asciidata.open(realProfilePath)
R_as, mag_R, emag_R = [], [], []
for ii in numpy.arange(len(inputData[0])):
  R_as.append(inputData[1][ii])
  mag_R.append(inputData[2][ii])
  emag_R.append(inputData[3][ii])

R_as_R, mag_R, emag_R = numpy.array(R_as), numpy.array(mag_R), numpy.array(emag_R)

# Extracting measured profile (Debattista+15?)

realProfilePath = 'photprofs_Iband.txt'

inputData = asciidata.open(realProfilePath)
R_as, mag_I = [], []
for ii in numpy.arange(len(inputData[0])):
  R_as.append(inputData[0][ii])
  mag_I.append(inputData[1][ii])

R_as_I, mag_I = numpy.array(R_as), numpy.array(mag_I)

# Interpolating curve and extracting value at defined R_eff
from scipy import interpolate
funct_tmp_R = scipy.interpolate.interp1d(R_as_R, mag_R)
value_at_sel_R = funct_tmp_R(dist_sel*gal_Reff)

funct_tmp_I = scipy.interpolate.interp1d(R_as_I, mag_I)
value_at_sel_I = funct_tmp_I(dist_sel*gal_Reff)


MaunaKea_sky_brightness = {'R': 20.3, 'V':21.1, 'I':19.2} #mag/arcsec
# based on CFHT observatory manual (http://www.cfht.hawaii.edu/Instruments/ObservatoryManual/CFHT_ObservatoryManual_(Sec_2).html)

fraction_R = 10.**((MaunaKea_sky_brightness['R']-value_at_sel_R)/2.5)

fraction_I = 10.**((MaunaKea_sky_brightness['I']-value_at_sel_I)/2.5)

