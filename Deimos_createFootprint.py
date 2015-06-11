## DEIMOS mask - create footprint file for Keck submission
#v_0.1
from Deimos_SKiMS_slit__def__ import *
import glob
###########################
## System-wide variables ##
###########################

__builtins__.verbosity = False	#If True, keeps the verbosity level high
__builtins__.graphicOutput = False	#If True, shows the final plot at the end of every ScoBen iteration

__builtins__.separationSlits = 0.4	#Standard separation between slits in arcsec
__builtins__.minWidthSlits = 3.		#Minimum slit width in arcsec 
__builtins__.coneAngle = 45.		#Angle of cone containing the slits

__builtins__.gal_Reff = 47.86		#Effective radius of galaxy

###############
## FUNCTIONS ##
###############

def convAngCoord(angString):
    dim=len(angString)
    sep=numpy.zeros(dim)
    for i in range(0,dim):  #cerca :
        if (angString[i] == ':'):
            sep[i]=1
    pos=numpy.nonzero(sep)
    #
    arcsec=float(angString[(pos[0][1])+1:])
    arcmin=float(angString[(pos[0][0]+1):(pos[0][1])])
    deg=float(angString[:(pos[0][0])])
    #
    arcsectot=(deg/abs(deg))*(abs(deg*3600)+(arcmin*60)+arcsec)
    deg_tot=(deg/abs(deg))*(abs(deg)+(arcmin/60)+(arcsec/3600))
    #
    return deg,arcmin,arcsec,arcsectot,deg_tot

def Deg2Str(degTot):
  sign = degTot/numpy.abs(degTot)
  dd = int(abs(degTot))
  amm = int(abs(degTot-dd)*60.)
  ass = ((abs(degTot-dd)*60.)-amm)*60
  if sign == 1: 
    return str(dd)+':'+str(amm)+':'+str(round(ass, 2))
  else:
    return str(dd*int(sign))+':'+str(amm)+':'+str(round(ass, 2))


def createReg(listSlits, nameFile='mask.reg', slitAngle = 5., slitWidth = 1., maskCoord = [0,0]):
  #
  maskVertices_Xrel = numpy.array([-498., -498.,360.,420.,460.,498., 498.,-498.,-498.,
                    -460.,-420.,-360., 259.7, 259.7, 259.7, 249.3, 249.3,
                     5.2,   5.2,  -5.2,  -5.2,-249.3,-249.3,-259.7,-259.7]) 
  maskVertices_Yrel = numpy.array([-146.,  146.,  146.,   95.,   52.,   -1., -146., -146.,   -1.,
               52.,   95.,  146.,  146., -146.,  146.,  146., -146., -146.,
               146.,  146., -146., -146.,  146.,  146., -146.])
        #[187., 479.,  479., 428.,385.,332.,187.,187.,332.,385.,428.,
         #      479.,479.,187.,479.,479.,187.,187.,479.,479.,187.,187.,479.,479.,187.]
  #
  maskCoordinates_Xabs = maskCoord[0]-(maskVertices_Xrel/3600.)
  maskCoordinates_Yabs = maskCoord[1]-(maskVertices_Yrel/3600.)
  #xMax, xMin = 498., -498.
  #yMax, yMin = 146., -146.
  fileOut = open(nameFile, 'wb')
  #
  fileOut.write("# "+nameFile+"\n")
  fileOut.write("# ds9 regions file generated\n# by Deimos_createFootprint.py\n\n")
  fileOut.write("# BEGIN SLITS")
  fileOut.write("global colour=red\nj2000\n")
  #
  lines = []
  counter = 0
  for ii in listSlits:
    slitRA_0, slitDec_0 = maskCoord[0]-(ii.x0/3600.), maskCoord[1]-(ii.y/3600.)
#    slitRA_c, slitDec_c = maskCoord[0]-(ii.centralCoords[0]/3600.), maskCoord[1]-(ii.centralCoords[1]/3600.)
    slitLength, slitWidth = ii.length, slitWidth
    slitAngle, slitNum = slitAngle, counter
    line = (["box", str(slitRA_0)+"d", str(slitDec_0)+"d", str(slitLength)+'"', 
#    line = (["box", str(slitRA_c)+"d", str(slitDec_c)+"d", str(slitLength)+'"', 
  	           str(slitWidth)+'"', str(int(slitAngle))+"d # color = red, text = {"+str(slitNum)+"}"])
    for item in line:
      fileOut.write("%s\t"% item)
    fileOut.write("\n")
    counter += 1
  #
  fileOut.write("# END SLITS\n# BEGIN TARGETS\n")
  #
  lines = []
  counter = 0
  for ii in listSlits:
    slitName = "SS_"+str(int(counter))
    slitRA_c, slitDec_c = maskCoord[0]-(ii.centralCoords[0]/3600.), maskCoord[1]-(ii.centralCoords[1]/3600.)
    line = ["diamond point", str(slitRA_c)+"d", str(slitDec_c)+"d # color = green, text = {"+(slitName)+"}"]
    for item in line:
      fileOut.write("%s\t"% item)
    fileOut.write("\n")
    counter += 1
  #
  fileOut.write("# END TARGETS\n# BEGIN MASK\n")
  #
  fileOut.write("polygon\t")
  line = []
  for ii in numpy.arange(len(maskCoordinates_Xabs)):
  	line.append(str(maskCoordinates_Xabs[ii])+"d")
  	line.append(str(maskCoordinates_Yabs[ii])+"d")
  #
  for item in line:
    fileOut.write("%s\t"% item)
  #
  fileOut.write("# color = red\n# END MASK\n# BEGIN TV\n")
  #
  # GUIDE STAR FRAME NOT NECESSARY
  #
  fileOut.write("# END TV\n")
  #
  fileOut.close()
  return True

def createDsimInput(listSlits, nameFile='mask.txt', slitAngle = 5., PAgal = 0., 
	      				maskCoord = [0, 0],
	      				extra = False, #To add extra slits near the edges
                rotAngle = 0, #To rotate the slits respect to the galaxy major axis PA
                offset=False):  #To add an x,y offset to the mask centre coordinates in degrees
  #
  xC, yC = maskCoord[0], maskCoord[1]
  angleRot_rad = numpy.pi/180.*((90.-PAgal)+rotAngle)
  #
  fileOut = open(nameFile, 'wb')
  #
  counter = 1
  minRA, maxRA, minDec, maxDec = 0, 0, 0, 0 
  for ii in listSlits:
#
    slitName = 'SS_'+str(int(counter))
    #
    # Rotation respect to the galaxy centre
    # (in arcsec)
    slitX_rot = (ii.centralCoords[0])*numpy.cos(angleRot_rad) - ii.centralCoords[1]*numpy.sin(angleRot_rad)
    slitY_rot = (ii.centralCoords[0])*numpy.sin(angleRot_rad) + ii.centralCoords[1]*numpy.cos(angleRot_rad)
    #
    absCoords_X = xC-(slitX_rot/3600.)/15.
    absCoords_Y = yC-(slitY_rot/3600.)
    #
    slitRA_0_string  = Deg2Str(absCoords_X) #in hours
    slitDec_0_string = Deg2Str(absCoords_Y) #in degrees
    line = ([slitName, slitRA_0_string, slitDec_0_string, str(2000), 
             str(20), 'r', str(int(counter)), '1', '1', int(slitAngle)+(PAgal)-rotAngle, 
             ]) 
    for item in line:
      fileOut.write("%s\t"% item)
    fileOut.write("\n")
    counter += 1
    #
    # For extra objects
    #
    if (ii.centralCoords[0]) < minRA: minRA = (ii.centralCoords[0])
    elif (ii.centralCoords[0]) > maxRA: maxRA = (ii.centralCoords[0])
    if ii.centralCoords[1] < minDec: minDec = ii.centralCoords[1]
    elif ii.centralCoords[1] > maxDec: maxDec = ii.centralCoords[1]
##    
  if extra:	#Creating extra targets near the edges of the mask (+-0.5arcmin)
    #
    for jj in numpy.arange(minDec-50., maxDec+50,  (10.)):	#Dec (arcsec)
      for ii in numpy.arange(minRA, minRA-150, -150./extra):	#RA (arcsec)
      # First side
        slitName = 'eSS_'+str(int(counter))
        #
        # Rotation respect to the galaxy centre
        #
        slitX_rot = ii*numpy.cos(angleRot_rad) - jj*numpy.sin(angleRot_rad)
        slitY_rot = ii*numpy.sin(angleRot_rad) + jj*numpy.cos(angleRot_rad)
        #
        absCoords_X = xC-(slitX_rot/3600.)/15.
        absCoords_Y = yC-(slitY_rot/3600.)
        #
        slitRA_0_string  = Deg2Str(absCoords_X) #in hours
        slitDec_0_string = Deg2Str(absCoords_Y) #in degrees
        #
        line = ([slitName, slitRA_0_string, slitDec_0_string, str(2000), 
  	               str(20), 'r', str(int(counter)), '2', '0', int(slitAngle)+(PAgal)-rotAngle, 
             ]) 
        for item in line:
          fileOut.write("%s\t"% item)
        fileOut.write("\n")
        counter += 1
      #
      for ii in numpy.arange(maxRA, maxRA+150, 150./extra):  #RA (arcsec)
        # Second side
        slitName = 'eSS_'+str(int(counter))
        #
        # Rotation respect to the galaxy centre
        #
        slitX_rot = ii*numpy.cos(angleRot_rad) - jj*numpy.sin(angleRot_rad)
        slitY_rot = ii*numpy.sin(angleRot_rad) + jj*numpy.cos(angleRot_rad)
        #
        absCoords_X = xC-(slitX_rot/3600.)/15.
        absCoords_Y = yC-(slitY_rot/3600.)
        #
        slitRA_0_string  = Deg2Str(absCoords_X) #in hours
        slitDec_0_string = Deg2Str(absCoords_Y) #in degrees
        #
        line = ([slitName, slitRA_0_string, slitDec_0_string, str(2000), 
  	               str(20), 'r', str(int(counter)), '2', '0', int(slitAngle)+(PAgal)-rotAngle, 
             ]) 
        for item in line:
          fileOut.write("%s\t"% item)
        fileOut.write("\n")
        counter += 1
  fileOut.close()
  #
  return True


##########
## MAIN ##
##########

# Retrieve input 
# Slits name, position in the mask, length and width
listInput = glob.glob('MaskObj*.dat')
DistMax = 1000.

# Retrieve the 5 best masks and choose the one with the 
# smallest maxdist.

for ii in listInput:
  try:
    inputFile = open(ii, 'rb')
    maskTMP, maxDist = pickle.load(inputFile)
    if maxDist < DistMax:
      DistMax = maxDist
      mask = maskTMP 
    inputFile.close()
  except:
  	print "Not able to open "+ii


### OPTION 1 -> Create all the files for submission

## Increase the x distances by amount rat
rat = 1.25  
## Duplicate slits inverting the x 
listSlits = mask.listSlits
listSlits2 = []
for ii in listSlits:
  slitTmp = Slit()
  ii.x0 *= rat
  ii.x1 *= rat
  ii.centralCoords[0] *= rat
  slitTmp.createSlit(-ii.x1, ii.length, y=-ii.y)
  listSlits2.append(slitTmp)

'''
to check why the inverted slits are overlapping (check plotting code)
'''

# Manually change Y of last two slits -> central
listSlits2[-1].y, listSlits2[-1].centralCoords[1] = 0, 0
listSlits[-1].y, listSlits[-1].centralCoords[1] = 0, 0
listSlits += listSlits2

# Create .fits file for submission
'''
'''
# Create .reg file for ds9
centreCoord = ['02:40:24.010','+39:03:47.83']
dummy = createReg(listSlits, nameFile='SS_N1023.reg', slitAngle = 5., 
	maskCoord = [convAngCoord(centreCoord[0])[4]*15.,convAngCoord(centreCoord[1])[4]])




### OPTION 2 -> Create only the catalog and then run dsim (simpler)
dummy = createDsimInput(listSlits, nameFile='SS_N1023.txt', slitAngle = 5., PAgal = 83.3,
	      				maskCoord = [convAngCoord(centreCoord[0])[4],convAngCoord(centreCoord[1])[4]],
	      				extra=100, offset=[0.,0.])



os.system('cp ./SS_N1023.txt ./Dsim/SS_N1023.txt')
os.remove('./Dsim/SSN1023.fits')
os.remove('./Dsim/SSN1023.out')





#### Create footprints for rotated masks
'''
dummy = createDsimInput(listSlits, nameFile='SS_N1023_rot45.txt', slitAngle = 5., PAgal = 83.3,
                maskCoord = [convAngCoord(centreCoord[0])[4],convAngCoord(centreCoord[1])[4]],
                extra=100, offset=[0.,0.], rotAngle = 0.)

 
To solve problem of slits getting further out from the centre with increasing angle in DSIM 
'''
### Just modifying the reg file for ds9
pathOutput = './Dsim/SSN1023_allangles.reg'
fileOut = open(pathOutput, 'wb')
#
fileOut.write("# \n")
fileOut.write("# ds9 regions file generated\n# by Deimos_createFootprint.py\n\n")
fileOut.write("# BEGIN SLITS")
fileOut.write("global colour=red\nj2000\n")

listAngles = [0, 45, 90, 135]
maskCoord = [convAngCoord(centreCoord[0])[4],convAngCoord(centreCoord[1])[4]]


lines = []
counter = 0

cc = ['red', 'green', 'blue', 'yellow']
for ii in listSlits:
 for jj in listAngles:
  if jj == listAngles[0]: colour = cc[0]
  elif jj == listAngles[1]: colour = cc[1]
  elif jj == listAngles[2]: colour = cc[2]
  elif jj == listAngles[3]: colour = cc[3]
  # 
  angleRot_rad = numpy.pi/180.*((90.-PAgal)+jj)
  #
  xRot = (ii.x0/3600.)*numpy.cos(angleRot_rad) - (ii.y/3600.)*numpy.sin(angleRot_rad)
  yRot = (ii.x0/3600.)*numpy.sin(angleRot_rad) + (ii.y/3600.)*numpy.cos(angleRot_rad) 
  #
  slitRA_0, slitDec_0 = (maskCoord[0]*15.-xRot), maskCoord[1]-yRot
  slitLength, slitWidth = ii.length, 1.
  slitAngle, slitNum = -((90.-PAgal)+jj)+5, counter
  line = (["box", str(slitRA_0)+"d", str(slitDec_0)+"d", str(slitLength)+'"', 
             str(slitWidth)+'"', str(int(slitAngle))+"d # "+
             "color = "+colour
 #            +", text = {"+str(slitNum)+
 #            "_"+str(int(jj))+"}"
             ])
  for item in line:
    fileOut.write("%s\t"% item)
  fileOut.write("\n")
 counter += 1

fileOut.write("# END SLITS\n# BEGIN TARGETS\n")

fileOut.close()

  #
  lines = []
  counter = 0
  for ii in listSlits:
    slitName = "SS_"+str(int(counter))
    slitRA_c, slitDec_c = maskCoord[0]-(ii.centralCoords[0]/3600.), maskCoord[1]-(ii.centralCoords[1]/3600.)
    line = ["diamond point", str(slitRA_c)+"d", str(slitDec_c)+"d # color = green"]
    #, text = {"+(slitName)+"}"]
    for item in line:
      fileOut.write("%s\t"% item)
    fileOut.write("\n")
    counter += 1
  #
  fileOut.write("# END TARGETS\n# BEGIN MASK\n")
  #
  fileOut.write("polygon\t")
  line = []
  for ii in numpy.arange(len(maskCoordinates_Xabs)):
    line.append(str(maskCoordinates_Xabs[ii])+"d")
    line.append(str(maskCoordinates_Yabs[ii])+"d")
  #
  for item in line:
    fileOut.write("%s\t"% item)
  #
  fileOut.write("# color = red\n# END MASK\n# BEGIN TV\n")
  #
  # GUIDE STAR FRAME NOT NECESSARY
  #
  fileOut.write("# END TV\n")
  #
  fileOut.close()
  return True