## DEIMOS mask - SKiMS fill
#v_0.1 - Assuming constant Surface Brightness
#v_0.2 - Sersic SB
#v_0.3 - reboot
#v_0.4 - addition Sky slits
#v_0.5 - addition variable Sersic profile given the mask PA

from Deimos_SKiMS_slit__def__v2 import *

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



# PLOTTING

## Single mask plot
listMasks[0][0].initializePlot(figSize=(10,5))
listMasks[0][0].plotMask(PAangle=listMasks[0][2])
listMasks[0][0].plotSlits(PAangle=listMasks[0][2])
#mask.plotBoundaries(PAangle=0.)

listMasks[1][0].ax = listMasks[0][0].ax
listMasks[1][0].plotMask(PAangle=listMasks[1][2])
listMasks[1][0].plotSlits(PAangle=listMasks[1][2])

listMasks[2][0].ax = listMasks[0][0].ax
listMasks[2][0].plotMask(PAangle=listMasks[2][2])
listMasks[2][0].plotSlits(PAangle=listMasks[2][2])

listMasks[3][0].ax = listMasks[0][0].ax
listMasks[3][0].plotMask(PAangle=listMasks[3][2])
listMasks[3][0].plotSlits(PAangle=listMasks[3][2])


savefig('AllMasks.pdf', bbox_inches='tight')










'''
### ISSUE WITH SAVING MASK DESIGN


for ii in listMasks:
  try:
    del ii[0].ax #Issues with pickling of object
  except:
    dummy = True

fileOut = open('MaskObj.dat', 'wb')
pickle.dump(listMasks[0][0], fileOut, pickle.HIGHEST_PROTOCOL)
fileOut.close()

dummy = raw_input("Best set of slits found. Press a key to continue.")

























### MASK Minor

__builtins__.mask_PA = 178.+90. # PA of the mask

mask_Minor, maxDist_Minor = findBestMask(iterations = 1., realProfilePath = 'Profiles/n2549_op_mos_sky_combscale.ell',
                bagal=gal_ba, PAgal=gal_PA, PAmask=mask_PA,
                mu_0=mu_0, pixelscale=1.2233)


#### PLOTS



### MASK Diagonal 1

__builtins__.mask_PA = 178.+45. # PA of the mask

mask_d1, maxDist_d1 = findBestMask(iterations = 100., realProfilePath = 'Profiles/n2549_op_mos_sky_combscale.ell',
                bagal=gal_ba, PAgal=gal_PA, PAmask=mask_PA,
                mu_0=mu_0, pixelscale=1.2233)


#### PLOTS
mask_d1.ax = mask_Major.ax
mask_d1.plotMask(PAangle=mask_PA)
mask_d1.plotSlits(PAangle=mask_PA)



### MASK Diagonal 2

__builtins__.mask_PA = 178.-45. # PA of the mask

mask_d2, maxDist_d2 = findBestMask(iterations = 100., realProfilePath = 'Profiles/n2549_op_mos_sky_combscale.ell',
                bagal=gal_ba, PAgal=gal_PA, PAmask=mask_PA,
                mu_0=mu_0, pixelscale=1.2233)


#### PLOTS
mask_d2.ax = mask_Major.ax
mask_d2.plotMask(PAangle=mask_PA)
mask_d2.plotSlits(PAangle=mask_PA)










mask.plotMask(PAangle=45.)
mask.plotSlits(PAangle=45.)
mask.plotBoundaries(PAangle=45.)



fileOut = open('MaskObj.dat', 'wb')
try:
  del mask.ax	#Issues with pickling of object
except:
  dummy = True

pickle.dump([mask, maxDist], fileOut, pickle.HIGHEST_PROTOCOL)
fileOut.close()

dummy = raw_input("Best set of slits found. Press a key to continue.")


########

fileIn = open('MaskObj.dat', 'rb')
mask, maxDist = pickle.load(fileIn)
fileIn.close()

#### PLOTS
## Single mask plot
mask.initializePlot(figSize=(10,5))
mask.plotMask(PAangle=0.)
mask.plotSlits(PAangle=0.)
mask.plotBoundaries(PAangle=0.)




mask.plotBins()
savefig('SingleSSKiMSmask.pdf', bbox_inches='tight')

# ## Multiple masks plot
# mask.initializePlot(figSize=(10,10))
# mask.plotMask(PAangle=0.)
# mask.plotSlits(PAangle=0.)
# mask.plotBoundaries(PAangle=0.)

# mask.plotMask(PAangle=90.)
# mask.plotSlits(PAangle=90.)
# mask.plotBoundaries(PAangle=90.)

# mask.plotMask(PAangle=135.)
# mask.plotSlits(PAangle=135.)
# mask.plotBoundaries(PAangle=135.)

# mask.plotMask(PAangle=45.)
# mask.plotSlits(PAangle=45.)
# mask.plotBoundaries(PAangle=45.)

# mask.ax.set_ylim([-550,550])

# savefig('MultipleSSKiMSmask.pdf', bbox_inches='tight')

# ## Single mask plot with de Vaucouleurs profile
# mask.initializePlot(figSize=(10,5))
# mask.plotMask(PAangle=0.)
# mask.plotSlits(PAangle=0.)
# mask.plotBoundaries(PAangle=0.)
# mask.plotBins()
# mask.ax.set_ylim([-500,180])


axes2 = mask.fig.add_axes([0.215,0.15,0.595,0.33])
for xlabel_i in axes2.get_xticklabels():
    xlabel_i.set_visible(False)
    xlabel_i.set_fontsize(0.0)

for xlabel_i in axes2.get_yticklabels():
    xlabel_i.set_fontsize(0.0)
    xlabel_i.set_visible(False)

for tick in axes2.get_xticklines():
    tick.set_visible(False)

for tick in axes2.get_yticklines():
    tick.set_visible(False)

xx = numpy.arange(-498.,498.,0.1)
yy = SersicFunct(numpy.abs(xx), 10, gal_Reff)
axes2.set_xlim([-498,498])
axes2.plot(xx, numpy.log10(yy), 'b-')
axes2.set_ylabel(r'$\mu$ [dex]')

savefig('SingleSSKiMSmask_DV.pdf', bbox_inches='tight')
'''
