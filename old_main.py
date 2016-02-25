## DEIMOS mask - SKiMS fill
#v_0.1 - Assuming constant Surface Brightness
#v_0.2 - Sersic SB
#v_0.3 - reboot
#v_0.4 - addition Sky slits
#v_0.5 - addition variable Sersic profile given the mask PA
#v_0.6 - fixed issue with mask saving
#v_0.7 - added radial limit for SuperSKiMS slits
#v_0.8 - fixed issue with SuperSKiMS high SB values

from Deimos_SKiMS_slit__def__v4 import *

###########################
## System-wide variables ##
###########################

'''
NGC 1407 - 2 masks - PA = 20 and 110 degrees
'''

__builtins__.verbosity = False	#If True, keeps the verbosity level high
__builtins__.graphicOutput = False	#If True, shows the final plot at the end of every ScoBen iteration

__builtins__.separationSlits = 0.4	#Standard separation between slits in arcsec
__builtins__.minWidthSlits = 3.		#Minimum slit width in arcsec

__builtins__.numberMasks = 2.
__builtins__.limitRadiusSlits = 5. # Max distance for SS slits in R_eff

__builtins__.gal_Reff = 63.   #Effective radius of galaxy  (from Brodie+14)
__builtins__.gal_ba = 1.-0.03   #Effective radius of galaxy  (from Brodie+14)
__builtins__.gal_RA = '03:40:11.8' # (from NED)
__builtins__.gal_Dec = '-18:34:48' # (from NED)
__builtins__.gal_PA = 35. # PA of the galaxy (from Brodie+14)


''' still to define in which band the mu_0 is'''
__builtins__.mu_0 =  15.# I-band SB of the galaxy in the centre (this is from Spolaor+08)

##########
## MAIN ##
##########

__builtins__.coneAngle = 180./(numberMasks) #Angle of cone containing the slits

numberIterations = 1e2

listMasks = []
for ii in np.arange(numberMasks):
  #__builtins__.mask_PA = gal_PA+ii*180./numberMasks # PA of the mask
  __builtins__.mask_PA = 20.+ii*180./numberMasks # PA of the mask
  mask_Tmp, maxDist_Tmp = findBestMask(iterations = numberIterations, #realProfilePath = profilePath,
  #  sersicParameters = [],
   # bagal=gal_ba, PAgal=gal_PA, PAmask=mask_PA,
    #            mu_0=mu_0, pixelscale=pixelscale
    )
  listMasks.append([mask_Tmp, maxDist_Tmp, mask_PA])
  #
  mask_Tmp.saveMaskSlits2txt(pathOutput='./SS_design_slits_'+str(ii)+'.txt')
  #
  mask_Tmp.createOutputDSIM(pathOutput='./catSS_'+str(round(mask_PA,1))+'.cat')
  #

fileOut = open('SS_design_masks.dat', 'wb')
pickle.dump(listMasks, fileOut, pickle.HIGHEST_PROTOCOL)
fileOut.close()


dummy = raw_input("Best set of slits found. Press a key to continue.")




# PLOTTING
fileIn = open('./SS_design_masks.dat', 'rb')
listMasks = pickle.load(fileIn)
fileIn.close()

## Single mask plot
listMasks[0][0].initializePlot(figSize=(10,5))
listMasks[0][0].plotMask(PAangle=listMasks[0][2])
listMasks[0][0].plotSlits(PAangle=listMasks[0][2])
#mask.plotBoundaries(PAangle=0.)

for ii in np.arange(len(listMasks)-1)+1:
  listMasks[ii][0].ax = listMasks[0][0].ax
  listMasks[ii][0].plotMask(PAangle=listMasks[ii][2])
  listMasks[ii][0].plotSlits(PAangle=listMasks[ii][2])


savefig('AllMasks.pdf', bbox_inches='tight')

