## DEIMOS mask - SKiMS fill
import asciidata, numpy, pylab, random, scipy, time, copy, os, glob
from scipy import optimize
try:
  import cPickle as pickle
except:
  import pickle

from pylab import *
from sys import stdout

# v.2 - 25/09/15  - Addition of extra option for Sersic profile
#                 - Addition 'saving catalog' function

###############
## Functions ##
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

    arcsectot=(deg/abs(deg))*(abs(deg*3600)+(arcmin*60)+arcsec)
    deg_tot=(deg/abs(deg))*(abs(deg)+(arcmin/60)+(arcsec/3600))
#
    return deg,arcmin,arcsec,arcsectot,deg_tot

def Coords2String(coordValue, #Converts degrees into string dd:amm:ass
                  hh=False): # If True, the output is in hh:mm:ss
    if hh:
      coordValue /= 15.
      signDD = ''
    else:
      if coordValue >= 0:
        signDD = '+'
      else:
        signDD = '-'
    #
    dd = int(coordValue)
    float_am = (coordValue-int(coordValue))*60.
    am = int(float_am)
    ass = (float_am-am)*60.
    #
    outString = signDD+str(dd)+':'+str(am)+':'+str(ass)
    return outString



def relativeCoordinates(RA, Dec, RA_gal, Dec_gal): #In degrees
  # From Huchra+91
  DeltaRA = numpy.sin(numpy.radians(RA-RA_gal))*numpy.cos(numpy.radians(Dec))
  DeltaDec = (numpy.sin(numpy.radians(Dec))*numpy.cos(numpy.radians(Dec_gal))-              numpy.cos(numpy.radians(RA-RA_gal))*numpy.cos(numpy.radians(Dec))*numpy.sin(numpy.radians(Dec_gal)))
  #
  return numpy.degrees(DeltaRA), numpy.degrees(DeltaDec)


def SersicFunct(r, I0, Re, m=4, bm = 7.67): #if m == 1 is an exponential function, m == 4 is a de Vaucouleur profile (UPDATED VERSION WITH I0)
  #bm = 7.67  #Cai et al. 2008
  r = numpy.array(r)*1.
  return I0 * (numpy.e**(-bm*((r/float(Re))**(1./m))))

def SersicFunctChi2(initialGuesses, R, I, eI): #Function that returns the chi2 -> For the Sersic Fit
  bm, Re, m = initialGuesses
#  bm = 7.67 #Cai et al. 2008
  R = numpy.array(R)*1.
  FittedI = I[0] * (numpy.e**(-bm*((R/float(Re))**(1./m))))
  #
  # Measure Chi2
  chi2 = numpy.sum(((FittedI-I)**2.))#/eI)
  return chi2

def withinMask(xP, yP):
  xP = numpy.abs(xP)
  if yP > -1 and xP > 360.:
    if (360. <= xP < 420.) and yP < -0.85*xP+452.:
      return True
    elif (420. <= xP < 460.) and yP < -1.075*xP+546.5:
      return True
    elif (460. <= xP < 498.) and yP < -1.9347368421*xP+693.5789473684:
      return True
    else:
      return False
  else:
    return True

def withinCones(xP, yP, q=0):
  xP = numpy.abs(xP)
  yline1 = (numpy.tan((coneAngle/2.)*numpy.pi/180.))*numpy.array(xP)+q
  yline2 = -(numpy.tan((coneAngle/2.)*numpy.pi/180.))*numpy.array(xP)+q
  if yline1 > yP > yline2:
    return True
  else:
    return False



def readProfileEllipse(filepath, pixelscale=1., mu_0 = 0.,
                        bagal=1., PAgal=90., PAmask=90.):
  # Reads SB profile from ELLIPSE output
  #
  # pixelscale is arcsec/pixel
  #
  inputData = np.loadtxt(filepath,
      dtype={ 'names': ('ROW', 'SMA', 'INTENS', 'INT_ERR', 'PIX_VAR', 'RMS', 'ELLIP', 'ELLIP_ERR', 'PA', 'PA_ERR', 'X0', 'X0_ERR', 'Y0', 'Y0_ERR', 'GRAD', 'GRAD_ERR', 'GRAD_R_ERR', 'RSMA', 'MAG', 'MAG_LERR', 'MAG_UERR', 'TFLUX_E', 'TFLUX_C', 'TMAG_E', 'TMAG_C', 'NPIX_E', 'NPIX_C', 'A3', 'A3_ERR', 'B3', 'B3_ERR', 'A4', 'A4_ERR', 'B4', 'B4_ERR', 'NDATA', 'NFLAG', 'NITER', 'STOP', 'A_BIG SAREA'),
              'formats': ('I', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', 'S10', )})
  #
  # Cleaning and floating
  R_arcsec = inputData['SMA']; R_arcsec = R_arcsec.astype('float')*pixelscale
  intens = inputData['INTENS']; intens[np.where(intens == 'INDEF')] = NaN; intens = intens.astype('float')
  errintens = inputData['INT_ERR']; errintens[np.where(errintens == 'INDEF')] = NaN; errintens = errintens.astype('float')
  #
  mag = -2.5*np.log10(intens/pixelscale**2.)
  emagp = (-2.5*np.log10((intens+errintens)/pixelscale**2.)) - mag
  emagm = mag-(-2.5*np.log10((intens-errintens)/pixelscale**2.))
  emag = np.sqrt(emagp**2.+emagm**2.)
  #
  # Scaling
  mag += (mu_0-mag[0])
  #
  theta = np.radians(PAgal-PAmask) #Angle respect to the galaxy major axis
  ## extraction profile
  R_theta = R_arcsec/np.sqrt((np.cos(theta))**2.+(np.sin(theta)*(1./bagal))**2.)  # Distances along the selected PA
  return R_theta, mag, emag

def readProfileTxt(filepath):
  # Reads SB profile from txt output
  #
  try:
    inputData = numpy.array(asciidata.open(realProfilePath))
    R_as, mag_R, emag_R = inputData[1,:], inputData[2,:], inputData[3,:]
  except:
    inputData = asciidata.open(realProfilePath)
    R_as, mag_R, emag_R = [], [], []
    for ii in numpy.arange(len(inputData[0])):
      R_as.append(inputData[1][ii])
      mag_R.append(inputData[2][ii])
      emag_R.append(inputData[3][ii])
      #
    R_as = numpy.array(R_as); mag_R = numpy.array(mag_R); emag_R = numpy.array(emag_R)
  #
  return R_as, mag_R, emag_R




def createGrid(xRange, yRange, resolution):
  [xMin, xMax], [yMin, yMax] = xRange, yRange
  gridX = numpy.empty((int((yMax-yMin)/resolution), int((xMax-xMin)/resolution)))
  gridY = numpy.empty((int((yMax-yMin)/resolution), int((xMax-xMin)/resolution)))
  for ii in numpy.arange(int((xMax-xMin)/resolution)):
    gridY[:, ii] = numpy.arange(yMin, yMax, resolution)
  for ii in numpy.arange(int((yMax-yMin)/resolution)):
    gridX[ii, :] = numpy.arange(xMin, xMax, resolution)
  return gridX, gridY


def findBestMask(iterations=100., maxDist = 1000., realProfilePath=False,
  sersicParameters = False, #sersicParameters = [bm, Re, m, I0]
  **kwargs): # includes 'pixelscale' (in case the profile is from ellipse) and galaxy and mask parameters (PAmask, PAgalaxy, bagalaxy and mu_0)
  import Deimos_SKiMS_slit__def__v2 as tst
  #
  if realProfilePath: #Extracts Sersic profile's parameters from input SB profile
    try:
      try:
        R_as, mag_R, emag_R = readProfileTxt(realProfilePath)
      except:
        R_as, mag_R, emag_R = readProfileEllipse(realProfilePath, pixelscale=kwargs['pixelscale'], PAmask=kwargs['PAmask'], bagal=kwargs['bagal'], mu_0=kwargs['mu_0'], PAgal=kwargs['PAgal'])
      #
      inputPar = [R_as, mag_R, emag_R]
      initialGuesses = [7.67, gal_Reff, 4] #Guesses for minimization, #b=7.67 Cai et al. 2008, Reff=47.9arcsec from Brodie+14, Sersic Index = 4 (de Vaucouleurs)
      #
      bm, Re, m = scipy.optimize.fmin(SersicFunctChi2, initialGuesses,
              args=(R_as, mag_R, emag_R), ftol = 0.1,  disp = True)
      I0 = mag_R[0]
    except:
      print "Error with SB profile. Using a default Sersic profile. "
      realProfilePath = False
  #
  for ii in arange(iterations):
    t1 = time.time()
    print "\n###########"
    print "Iteration "+str(int(ii+1))+"/"+str(int(iterations))+"\n"
    print "Creating mask..."
    tmpObj = tst.Mask()
    print "\r DONE!"
    print "Creating slits..."
    if realProfilePath:
      tmpObj.createSlits(sersicPar=[bm, Re, m, I0])
    elif sersicParameters:
      tmpObj.createSlits(sersicPar=sersicParameters)
    else:
      tmpObj.createSlits()
    print "\r DONE!"
    print "Finding largest empty space between the slits."
    tmpDist = tmpObj.getMaxEmptySpace()
    print "\r DONE!"
    if tmpDist < maxDist:
      maxDist = tmpDist
      ## Adding Sky Slits
      tmpObj.createSkySlits()
      #
      Mask = tmpObj
    t2 = time.time()
    print "Elapsed time: "+str(t2-t1)
    print "###########\n"
    tmpObj.__del__()
  return Mask, maxDist

#############
## Classes ##
#############

class generalClass(object):
#
  def __init__(self):
    if verbosity: print "Object "+str(self)+" created."
#
  def __del__(self):
    if verbosity: print "Object "+str(self)+" deleted."


class Mask(generalClass):
#
  def __init__(self):
    generalClass.__init__(self)
    self.xM = xM = [-498.,-498.,360.,420.,460.,498., 498.,-498.,-498.,
                    -460.,-420.,-360., 259.7, 259.7, 259.7, 249.3, 249.3,
                     5.2,   5.2,  -5.2,  -5.2,-249.3,-249.3,-259.7,-259.7]
    self.yM = [-146.,  146.,  146.,   95.,   52.,   -1., -146., -146.,   -1.,
               52.,   95.,  146.,  146., -146.,  146.,  146., -146., -146.,
               146.,  146., -146., -146.,  146.,  146., -146.]
        #[187., 479.,  479., 428.,385.,332.,187.,187.,332.,385.,428.,
         #      479.,479.,187.,479.,479.,187.,187.,479.,479.,187.,187.,479.,479.,187.]
    self.xMax, self.xMin = 498., -498.
    self.yMax, self.yMin = 146., -146.
    self.posCurrent = 0
    self.spaceSkySlits = 50 # Arcsec from the edges dedicated to sky slits
#
    self.listSlits = []
    self.fig = None
#
  def createSlits(self, sersicPar=False, #[bm, Re, m]
                        SNmin=35.):
    self.posCurrent = separationSlits/2.
    while self.posCurrent < numpy.max(self.xM)-self.spaceSkySlits:
      if sersicPar:
        slitwidth = 1.
        exptime = 120.*60. #(2 hours)
        extraConstant = 4e6 #To convert into real signal -> manually changed...
        FluxDensity = extraConstant*exptime*10.**(SersicFunct(self.posCurrent, sersicPar[3], sersicPar[1], bm=sersicPar[0], m=sersicPar[2])/(-2.5))
        lengthSlitMeasured = (SNmin**2./(slitwidth*FluxDensity))
        lengthSlit = numpy.max([minWidthSlits, lengthSlitMeasured])
      else:
        lengthSlit = numpy.max([minWidthSlits, 0.0075/SersicFunct(self.posCurrent, 10, gal_Reff)])
        FluxDensity = SersicFunct(self.posCurrent, 10, gal_Reff)
      tmpObj = Slit()
      tmpObj.createSlit(self.posCurrent, lengthSlit, SB=FluxDensity)
      if self.posCurrent + lengthSlit < max(self.xM)-self.spaceSkySlits:
        self.listSlits.append(tmpObj)
        self.posCurrent += lengthSlit + separationSlits
      else:
        tmpObj = Slit()
        tmpObj.createSlit(self.posCurrent, max(self.xM)-self.spaceSkySlits-self.posCurrent-separationSlits/2.)
        self.listSlits.append(tmpObj)
        self.posCurrent += lengthSlit + separationSlits
    #
#
  def createSkySlits(self, lengthSkySlit = 3.):
    #
    self.posCurrent = max(self.xM)-self.spaceSkySlits
    #
    while self.posCurrent < numpy.max(self.xM):
      tmpObj = Slit()
      tmpObj.createSlit(self.posCurrent, lengthSkySlit, target='Sky')
      self.listSlits.append(tmpObj)
      self.posCurrent += lengthSkySlit + separationSlits
#
  def getMaxEmptySpace(self, resolution = 1.): # Creates circles and finds the maximum size circle that can be built without touching the slits
    #Define list of points within the cones
    gridPointsX, gridPointsY = createGrid([0., self.xMax], [self.yMin, self.yMax], resolution)
    ListCircles = []
    xtmp, ytmp, rtmp = [], [], []
    for ii in numpy.arange(len(gridPointsX)):
      for jj in numpy.arange(len(gridPointsX[0])):
        if withinCones(gridPointsX[ii, jj], gridPointsY[ii, jj]) and withinMask(gridPointsX[ii, jj], gridPointsY[ii, jj]):  #Check if point within the angle cone and mask edges
          distances = []
          for kk in self.listSlits:
            distances.append(numpy.sqrt((gridPointsX[ii, jj]-kk.centralCoords[0])**2+(gridPointsY[ii, jj]-kk.centralCoords[1])**2))
          ListCircles.append(numpy.nanmin(distances))
          xtmp.append(gridPointsX[0][jj]); ytmp.append(gridPointsY[ii][0]); rtmp.append(numpy.nanmin(distances))
    return numpy.max(ListCircles)
#
  def initializePlot(self, figSize=(10,5)):
    plt.ion()
    self.fig = figure(num=0, figsize=figSize)
    clf()
    self.ax = subplot(111)
    self.ax.set_aspect('equal')
    self.ax.set_xlabel(r'$\Delta$x [arcsec]')
    self.ax.set_ylabel(r'$\Delta$y [arcsec]')
#

  def plotMask(self, PAangle=0.):
    theta = np.radians(PAangle)
    xMrot = numpy.array(self.xM)*numpy.cos(theta) - numpy.array(self.yM)*numpy.sin(theta)
    yMrot = numpy.array(self.xM)*numpy.sin(theta) + numpy.array(self.yM)*numpy.cos(theta)
#    self.ax.plot(self.xM, self.yM, 'k-')
    self.ax.plot(xMrot, yMrot, 'k-')
#
  def plotSlits(self, PAangle=0.):
    theta = np.radians(PAangle)
    # Separating Science slits from Sky slits
    selSKiMS, selSky = [], []
    for ii in self.listSlits:
      if ii.type == 'SKiMS':
        selSKiMS.append(ii)
      elif ii.type == 'Sky':
        selSky.append(ii)
    #
    for ii in selSKiMS:
      xSlitR = numpy.array([ii.x0, ii.x1])*numpy.cos(theta) - numpy.array([ii.y, ii.y])*numpy.sin(theta)
      ySlitR = numpy.array([ii.x0, ii.x1])*numpy.sin(theta) + numpy.array([ii.y, ii.y])*numpy.cos(theta)
      xSlitL = numpy.array([-ii.x0, -ii.x1])*numpy.cos(theta) - numpy.array([-ii.y, -ii.y])*numpy.sin(theta)
      ySlitL = numpy.array([-ii.x0, -ii.x1])*numpy.sin(theta) + numpy.array([-ii.y, -ii.y])*numpy.cos(theta)
      self.ax.plot(xSlitR, ySlitR, '-r')
      self.ax.plot(xSlitL, ySlitL, '-r')
      #
    for ii in selSky:
      xSlitR = numpy.array([ii.x0, ii.x1])*numpy.cos(theta) - numpy.array([ii.y, ii.y])*numpy.sin(theta)
      ySlitR = numpy.array([ii.x0, ii.x1])*numpy.sin(theta) + numpy.array([ii.y, ii.y])*numpy.cos(theta)
      xSlitL = numpy.array([-ii.x0, -ii.x1])*numpy.cos(theta) - numpy.array([ii.y, ii.y])*numpy.sin(theta)
      ySlitL = numpy.array([-ii.x0, -ii.x1])*numpy.sin(theta) + numpy.array([ii.y, ii.y])*numpy.cos(theta)
      self.ax.plot(xSlitR, ySlitR, '-b')
      self.ax.plot(xSlitL, ySlitL, '-b')
      #
#
  def plotBoundaries(self, q=0, PAangle=0.):
    theta = PAangle * numpy.pi/180.
    xlimits = [self.xMin, self.xMax]
    yline1 = (numpy.tan((coneAngle/2.)*numpy.pi/180.+theta))*numpy.array(xlimits)+q
    yline2 = -(numpy.tan((coneAngle/2.)*numpy.pi/180.+theta))*numpy.array(xlimits)+q
    line1 = self.ax.plot(xlimits, yline1, 'g:')
    line2 = self.ax.plot(xlimits, yline2, 'g:')
#
  def plotBins(self): #Only for PA=0
    for ii in self.listSlits:
      self.ax.fill_between([ii.x0, ii.x1], [self.yMin,self.yMin], [self.yMax,self.yMax], color='black', alpha=0.2)
      self.ax.fill_between([-ii.x0, -ii.x1], [self.yMin,self.yMin], [self.yMax,self.yMax], color='black', alpha=0.2)
#
  def saveMaskSlits2txt(self, pathOutput=''):
    if not(pathOutput):
      pathOutput = './mask'+str(round(mask_PA))+'.txt'
    theta = np.radians(mask_PA)
    Xmin, Xmax, Ymin, Ymax = [], [], [], []
    for ii in self.listSlits:
      for jj in [-1., 1.]:
        Xmin.append(jj*ii.x0*np.cos(theta)-jj*ii.y*np.sin(theta))
        Xmax.append(jj*ii.x1*np.cos(theta)-jj*ii.y*np.sin(theta))
        Ymin.append(jj*ii.x0*np.sin(theta)+jj*ii.y*np.cos(theta))
        Ymax.append(jj*ii.x1*np.sin(theta)+jj*ii.y*np.cos(theta))
    #
    np.savetxt(pathOutput, np.transpose((Xmin, Xmax, Ymin, Ymax)),
              delimiter='\t', header='Xmin\tXmax\tYmin\tYmax')
#
  def createOutputDSIM(self, pathOutput='./catSS.txt'):
    # It creates a catalog object for DSIM
    # All the slits are placed/rotated respect to their center and the 2 semi-  lengths are provided
    RAgal_deg, Decgal_deg = convAngCoord(gal_RA)[4]*15., convAngCoord(gal_Dec)[4]
    counter = 1
    listStrings = []
    for ii in self.listSlits:
      for side in ['L', 'R']:
        # Rotate with mask
        rotAngle = np.radians(90.-mask_PA)
        # RIGHT HAND SLIT
        if side == 'R':
          xSlit = numpy.array([ii.x0, ii.x1])*numpy.cos(rotAngle) - numpy.array([  ii.y, ii.y])*numpy.sin(rotAngle)
          ySlit = numpy.array([ii.x0, ii.x1])*numpy.sin(rotAngle) + numpy.array([  ii.y, ii.y])*numpy.cos(rotAngle)
        elif side == 'L':
          xSlit = numpy.array([-ii.x0, -ii.x1])*numpy.cos(rotAngle) - numpy.array  ([-ii.y, -ii.y])*numpy.sin(rotAngle)
          ySlit = numpy.array([-ii.x0, -ii.x1])*numpy.sin(rotAngle) + numpy.array  ([-ii.y, -ii.y])*numpy.cos(rotAngle)
        #
        # Assign coordinates
        posX, posY = np.average([xSlit]), np.average([ySlit])
        RAslit, Decslit = RAgal_deg-posX/3600., Decgal_deg-posY/3600.
        #
        if ii.type == 'SKiMS':
          nameSlit_str = 'SS_'+str(int(mask_PA))+'_'+str(counter)
          mag_str = str(round(ii.SB, 2))
          pcode_str = str(int(1001-counter))
          sample_str = '1'
          select_str = '1' #Pre-selected
        elif ii.type == 'Sky':
          nameSlit_str = 'SS_sky_'+str(int(mask_PA))+'_'+str(counter)
          mag_str = '0'
          pcode_str = str(int(300-counter))
          sample_str = '3'
          select_str = '1' #Pre-selected
        #
        RASlit_str = Coords2String(RAslit, hh=True)
        DecSlit_str = Coords2String(Decslit)
        equinox_str = '2000'
        passband_str = 'R'  #Not true
        slitPA_str = str(round(mask_PA+5., 2))
        len1_str, len2_str = str(round(ii.length/2., 2)), str(round(ii.length/2.,2))
        slitWidth_str = '1'
        #
        #SAVE
        listStrings.append(nameSlit_str+'\t'+RASlit_str+'\t'+DecSlit_str+'\t'+  equinox_str+'\t'+mag_str+'\t'+passband_str+'\t'+pcode_str+'\t'+  sample_str+'\t'+select_str+'\t'+slitPA_str+'\t'+len1_str+'\t'+len2_str+  '\t'+slitWidth_str+'\n')
        counter += 1
        #
    outputTable =  np.column_stack(listStrings)
    np.savetxt(pathOutput, outputTable, delimiter="", fmt="%s")
    return True



class Slit(generalClass):
#
  def __init__(self):
    self.x0, self.x1 = 0, 0
    self.y = 0
    self.SB = 0
    self.length = 0
    self.centralCoords = [0,0]
    self.type = 'NaN'
#
  def createSlit(self, pos0, length, y='random', xRange=[-498.,498.], yRange=[-146.,146.], target='SKiMS', SB=0.):
    self.x0 = pos0
    self.length = length
    self.x1 = self.x0 + self.length
    self.SB = SB
    self.type = target
    if y == 'random':
      import random
      elevation = (self.x0+self.length/2.)*numpy.tan((coneAngle/2.)*numpy.pi/180.)
      flag = True
      while flag:
        y = (random.random()*2.*elevation)-elevation
        if (yRange[0] < y < yRange[1]) and (withinMask(self.x1, y)):
          self.y = y
          flag = False
    else:
      self.y = y
    self.centralCoords = [numpy.mean([self.x0, self.x1]), self.y]










