## DEIMOS mask - SKiMS fill
#v_0.1 - Assuming constant Surface Brightness
#v_0.2 - Sersic SB
#v_0.3 - reboot
#v_0.4 - addition Sky slits
#v_0.5 - addition variable Sersic profile given the mask PA
#v_0.6 - fixed issue with mask saving
#v_0.7 - added radial limit for SuperSKiMS slits
#v_0.8 - fixed issue with SuperSKiMS high SB values

import numpy as np
import pickle

from utils import findBestMask


###########################
## System-wide variables ##
###########################

name = 'n4551'
outdir = './n4551/'

separationSlits = 0.4	#Standard separation between slits in arcsec
minWidthSlits = 3.		#Minimum slit width in arcsec

numberMasks = 2.
# limitRadiusSlits = 5. # Max distance for SS slits in R_eff

gal_Reff = 16.6   #Effective radius of galaxy  (from Atlas3d)
gal_ba = 0.75   #Effective radius of galaxy  (from Atlas3d)
gal_RA = '12:35:37.9' # from NED
gal_Dec = '+12:15:50' # (from NED)
gal_PA = 180 + 70.5 # PA of the galaxy (from Atlas3d), rotated to match ideal alignment stars


##########
## MAIN ##
##########

coneAngle = 180./(numberMasks) #Angle of cone containing the slits

numberIterations = 1e3

listMasks = []
for ii in np.arange(numberMasks):
    mask_PA = gal_PA+ii*180./numberMasks # PA of the mask
    mask_Tmp, maxDist_Tmp = findBestMask(iterations=numberIterations,
                                         separationSlits=separationSlits,
                                         minWidthSlits=minWidthSlits,                                         
                                         gal_Reff=gal_Reff,
                                         gal_ba=gal_ba,
                                         gal_RA=gal_RA,
                                         gal_Dec=gal_Dec,
                                         gal_PA=gal_PA,
                                         mask_PA=mask_PA,
                                         coneAngle=coneAngle)
    listMasks.append([mask_Tmp, maxDist_Tmp, mask_PA])
    #
    mask_Tmp.saveMaskSlits2txt(pathOutput=outdir + 'SS_+' + name + '_' + str(ii) + '.txt')
    #
    mask_Tmp.createOutputDSIM(pathOutput=outdir + 'SS_' + name + '_' + str(round(mask_PA,1)) + '.in')
    #

fileOut = open('SS_design_masks.dat', 'wb')
pickle.dump(listMasks, fileOut, pickle.HIGHEST_PROTOCOL)
fileOut.close()


