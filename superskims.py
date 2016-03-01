## DEIMOS mask - SKiMS fill
#v_0.1 - Assuming constant Surface Brightness
#v_0.2 - Sersic SB
#v_0.3 - reboot
#v_0.4 - addition Sky slits
#v_0.5 - addition variable Sersic profile given the mask PA
#v_0.6 - fixed issue with mask saving
#v_0.7 - added radial limit for SuperSKiMS slits
#v_0.8 - fixed issue with SuperSKiMS high SB values

import os
import pickle
import numpy as np


from utils import findBestMask


###########################
## System-wide variables ##
###########################

# # NGC 4551
# name = 'n4551'
# outdir = '../../gc_selection_march2016/n4551/'
# gal_Reff = 16.6   #Effective radius of galaxy  (from Atlas3d)
# gal_ba = 0.75   #Effective radius of galaxy  (from Atlas3d)
# gal_RA = '12:35:37.9' # from NED
# gal_Dec = '+12:15:50' # (from NED)
# gal_PA = 180 + 70.5 # PA of the galaxy (from Atlas3d), rotated to match ideal alignment stars

# # NGC 4550
# name = 'n4550'
# outdir = '../../gc_selection_march2016/n4550/'
# gal_Reff = 15.5    #Effective radius of galaxy  (from Atlas3d)
# gal_ba = 0.32   #Effective radius of galaxy  (from Atlas3d)
# gal_RA = '12:35:30.6' # from NED
# gal_Dec = '+12:13:15' # (from NED)
# gal_PA = 179 # PA of the galaxy (from Atlas3d), rotated to match ideal alignment stars

# # n4387
# name = "n4387"
# outdir = "../../gc_selection_march2016/n4387/"
# gal_Reff = 15.85   #Effective radius of galaxy in arcsec (from Atlas3d)
# gal_ba = 0.63   # axial ratio of galaxy  (from Atlas3d)
# gal_RA = '12:25:41.7151' # from Atlas3d
# gal_Dec = '+12:48:37.2924' # from Atlas3d
# gal_PA = 143.4 # PA of the galaxy (from Atlas3d), rotated to match ideal alignment stars

# # n3665
# name = 'n3665'
# outdir = '../../gc_selection_march2016/n3665/'
# gal_Reff = 30.9   #Effective radius of galaxy in arcsec (from Atlas3d)
# gal_ba = 0.78   # axial ratio of galaxy  (from Atlas3d)
# gal_RA = '11:24:43.7' # from NED
# gal_Dec = '+38:45:46' # from NED
# gal_PA = 30.9 # PA of the galaxy (from Atlas3d)

# n5475
name = 'n5475'
outdir = '../../gc_selection_march2016/n5475/'
gal_Reff = 16.6   #Effective radius of galaxy in arcsec (from Atlas3d)
gal_ba = 0.3   # axial ratio of galaxy  (from Atlas3d)
gal_RA = '14:05:12.4' # from NED
gal_Dec = '+55:44:31' # from NED
gal_PA = 166.2 # PA of the galaxy (from Atlas3d)



# general settings

separationSlits = 0.4	#Standard separation between slits in arcsec
minWidthSlits = 3.		#Minimum slit width in arcsec

numberMasks = 2.
coneAngle = 180./(numberMasks) #Angle of cone containing the slits
major_angle = 60
minor_angle = 120

numberIterations = 1e3



##########
## MAIN ##
##########

if not os.path.isdir(outdir):
    os.mkdir(outdir)
    
listMasks = []
major_axis = True

for ii in np.arange(numberMasks):
    mask_PA = gal_PA+ii*180./numberMasks # PA of the mask
    # for two masks
    if major_axis:
        mask_cone = major_angle
        major_axis = False
    else:
        mask_cone = minor_angle
    
    mask_Tmp, maxDist_Tmp = findBestMask(iterations=numberIterations,
                                         separationSlits=separationSlits,
                                         minWidthSlits=minWidthSlits,
                                         gal_Reff=gal_Reff,
                                         gal_ba=gal_ba,
                                         gal_RA=gal_RA,
                                         gal_Dec=gal_Dec,
                                         gal_PA=gal_PA,
                                         mask_PA=mask_PA,
                                         coneAngle=mask_cone)
    listMasks.append([mask_Tmp, maxDist_Tmp, mask_PA])
    # mask_Tmp.saveMaskSlits2txt(pathOutput=outdir + 'SS_' + name +
    #                            '_' + str(ii) + '.txt')
    mask_Tmp.createOutputDSIM(pathOutput=outdir + 'SS_' + name
                              + '_' + str(round(mask_PA,1)) + '.in')
    mask_Tmp.write_regions(outdir + 'SS_' + name + '_' +
                           str(round(mask_PA, 1)) + '.reg')

# with open(outdir + 'SS_design_masks.dat', 'wb') as fileOut:
#     pickle.dump(listMasks, fileOut, pickle.HIGHEST_PROTOCOL)



