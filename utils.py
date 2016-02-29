## DEIMOS mask - SKiMS fill
import sys
import time
import numpy as np


# v.2 - 25/09/15  - Addition of extra option for Sersic profile
#                 - Addition 'saving catalog' function
# v.3 - 29/09/15  - Fixed issue with single mask
#                 - Added radial limit for slits
# v.4 - 16/11/15  - Fixed issue with too high SB magnitudes
# v.5 - 28/01/16  - Allow the option for not saving the slit length

###############
## Functions ##
###############

def convAngCoord(angString):
    dim=len(angString)
    sep=np.zeros(dim)
    for i in range(0,dim):  #cerca :
        if (angString[i] == ':'):
            sep[i]=1
    pos=np.nonzero(sep)
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
    outString = signDD+str(np.abs(dd))+':'+str(np.abs(am))+':'+str(np.abs(ass))
    return outString


def relativeCoordinates(RA, Dec, RA_gal, Dec_gal):
    #In degrees
    # From Huchra+91
    DeltaRA = np.sin(np.radians(RA-RA_gal))*np.cos(np.radians(Dec))
    DeltaDec = (np.sin(np.radians(Dec))*np.cos(np.radians(Dec_gal)) -
                np.cos(np.radians(RA-RA_gal))*np.cos(np.radians(Dec)) *
                np.sin(np.radians(Dec_gal)))
    return np.degrees(DeltaRA), np.degrees(DeltaDec)


def SersicFunct(r, I0, Re, m=4, bm = 7.67):
    #if m == 1 is an exponential function, m == 4 is a de Vaucouleur profile (UPDATED VERSION WITH I0)
    #bm = 7.67  #Cai et al. 2008
    r = np.array(r)*1.
    return I0 * (np.e**(-bm*((r/float(Re))**(1./m))))


def withinMask(xP, yP):
    xP = np.abs(xP)
    a = (yP < -1) | (xP < 360)
    b = (360 <= xP) & (xP < 420) & (yP < -0.85*xP+452.)
    c = (420. <= xP) & (xP < 460.) & (yP < -1.075*xP+546.5)
    d = (460. <= xP) & (xP < 498.) & (yP < -1.9347368421*xP+693.5789473684)
    return a | b | c | d


def withinCones(xP, yP, coneAngle, q=0):
    xP = np.abs(xP)
    yline1 = (np.tan((coneAngle/2.)*np.pi/180.))*np.array(xP)+q
    yline2 = -(np.tan((coneAngle/2.)*np.pi/180.))*np.array(xP)+q
    return (yline2 < yP) & (yP < yline1)


def withinSlitRange(xP, yP, maxRadius):
    return np.sqrt(xP**2.+yP**2) <= maxRadius


class Mask:
    
    def __init__(self, gal_Reff, gal_ba, gal_RA, gal_Dec, gal_PA, mask_PA,
                 coneAngle, separationSlits, minWidthSlits,
                 limitRadiusSlits=5):
        self.xM = [-498.,-498.,360.,420.,460.,498., 498.,-498.,-498.,
                   -460.,-420.,-360., 259.7, 259.7, 259.7, 249.3, 249.3,
                   5.2,   5.2,  -5.2,  -5.2,-249.3,-249.3,-259.7,-259.7]
        self.yM = [-146.,  146.,  146.,   95.,   52.,   -1., -146., -146.,   -1.,
                   52.,   95.,  146.,  146., -146.,  146.,  146., -146., -146.,
                   146.,  146., -146., -146.,  146.,  146., -146.]
        #[187., 479.,  479., 428.,385.,332.,187.,187.,332.,385.,428.,
        #      479.,479.,187.,479.,479.,187.,187.,479.,479.,187.,187.,479.,479.,187.]
        self.xMax, self.xMin = 498., -498.
        self.yMax, self.yMin = 146., -146.
        self.gal_Reff = gal_Reff
        self.gal_ba = gal_ba
        self.gal_RA = gal_RA
        self.gal_Dec = gal_Dec
        self.gal_PA = gal_PA
        self.mask_PA = mask_PA
        self.coneAngle = coneAngle
        self.separationSlits = separationSlits
        self.minWidthSlits = minWidthSlits
        self.posCurrent = 0
        self.spaceSkySlits = 50 # Arcsec from the edges dedicated to sky slits
        self.limitRadiusSlits = limitRadiusSlits
        self.max_dist_SKiMS_slits = 0.
        self.listSlits = []
    
        # R_eff along the mask direction    
        Reff_major = np.sqrt((gal_Reff**2.)/gal_ba) # Reff along the major axis
        Reff_minor = np.sqrt((gal_Reff**2.)*gal_ba) # Reff along the minor axis
        DPA_mask_gal = np.radians(gal_PA-mask_PA)
        self.Reff_PA = np.sqrt((Reff_major*np.cos(DPA_mask_gal))**2.+(Reff_minor*np.sin(DPA_mask_gal))**2.)


    def createSlits(self, SNmin=35.):
        self.posCurrent = self.separationSlits/2.
        while self.posCurrent < self.limitRadiusSlits*self.Reff_PA:#np.max(self.xM)-self.spaceSkySlits:
            lengthSlit = np.max([self.minWidthSlits, 0.0075/SersicFunct(self.posCurrent, 10, self.gal_Reff)])
            FluxDensity = SersicFunct(self.posCurrent, 10, self.gal_Reff)
            tmpObj = Slit(self.posCurrent, lengthSlit, self.coneAngle, SB=FluxDensity)
            if self.posCurrent + lengthSlit < self.limitRadiusSlits*self.Reff_PA:#max(self.xM)-self.spaceSkySlits:
                self.listSlits.append(tmpObj)
                self.posCurrent += lengthSlit + self.separationSlits
            else:
                tmpObj = Slit(self.posCurrent,
                              self.limitRadiusSlits * self.Reff_PA - self.posCurrent - self.separationSlits/2.,
                              self.coneAngle)
                self.listSlits.append(tmpObj)
                self.posCurrent += lengthSlit + self.separationSlits
                self.max_dist_SKiMS_slits = self.posCurrent
        self.xcoords = np.array([slit.centralCoords[0] for slit in self.listSlits])
        self.ycoords = np.array([slit.centralCoords[1] for slit in self.listSlits])

                
    def createSkySlits(self, lengthSkySlit = 3.):
        self.posCurrent = max(self.xM)-self.spaceSkySlits
        while self.posCurrent < np.max(self.xM):
            tmpObj = Slit(self.posCurrent, lengthSkySlit, self.coneAngle, target='Sky')
            self.listSlits.append(tmpObj)
            self.posCurrent += lengthSkySlit + self.separationSlits


    def getMaxEmptySpace(self, resolution = 0.1):
        # Creates circles and finds the maximum size circle that can be built without touching the slits
        #Define list of points within the cones

        x_samples = np.arange(0, self.xMax, resolution)
        y_samples = np.arange(self.yMin, self.yMax, resolution)

        x = np.tile(x_samples, y_samples.shape)
        y = np.tile(y_samples, (x_samples.size, 1)).T.flatten()
        
        pick = (withinCones(x, y, self.coneAngle) &
                withinSlitRange(x, y, self.max_dist_SKiMS_slits) &
                withinMask(x, y))
        
        x_slits = np.tile(self.xcoords, (x[pick].size, 1)).T
        y_slits = np.tile(self.ycoords, (y[pick].size, 1)).T

        nslits = len(self.listSlits)
        
        x_points = np.tile(x[pick], (nslits, 1))
        y_points = np.tile(y[pick], (nslits, 1))
        
        distances = np.amin(np.sqrt((x_slits - x_points)**2 +
                                    (y_slits - y_points)**2),
                            axis=0)
        return np.sum(distances)

    
    def saveMaskSlits2txt(self, pathOutput=''):
        if not(pathOutput):
            pathOutput = './mask'+str(round(self.mask_PA))+'.txt'
        theta = np.radians(self.mask_PA)
        Xmin, Xmax, Ymin, Ymax = [], [], [], []
        for ii in self.listSlits:
            for jj in [-1., 1.]:
                Xmin.append(jj*ii.x0*np.cos(theta)-jj*ii.y*np.sin(theta))
                Xmax.append(jj*ii.x1*np.cos(theta)-jj*ii.y*np.sin(theta))
                Ymin.append(jj*ii.x0*np.sin(theta)+jj*ii.y*np.cos(theta))
                Ymax.append(jj*ii.x1*np.sin(theta)+jj*ii.y*np.cos(theta))
        np.savetxt(pathOutput, np.transpose((Xmin, Xmax, Ymin, Ymax)),
                   delimiter='\t', header='Xmin\tXmax\tYmin\tYmax')


    def write_regions(self, output):
        '''
        Save the current mask design to a ds9 regions file.

        Parameters
        ----------
        output: str, filename to save regions
        
        Returns
        -------
        None
        '''
        with open(output, 'w') as f:
            f.write('# Region file format: DS9 version 4.1\n')
            f.write('global color=red move=0 select=0\n')
            f.write('j2000\n')
            RAgal_deg, Decgal_deg = convAngCoord(self.gal_RA)[4]*15., convAngCoord(self.gal_Dec)[4]
            counter = 1
            
            for slit in self.listSlits:
                # x_slit, y_slit = slit.centralCoords
                # rotAngle = np.radians(90. - self.mask_PA)
                # RAslit = RAgal_deg + (-x_slit * np.cos(rotAngle) +
                #                        y_slit * np.sin(rotAngle)) / 3600.
                # Decslit = Decgal_deg + (x_slit * np.sin(rotAngle) +
                #                          y_slit * np.cos(rotAngle)) / 3600.
                for side in ['L', 'R']:
                    # Rotate with mask
                    rotAngle = np.radians(90.-self.mask_PA)
                    # RIGHT HAND SLIT
                    if side == 'R':
                        xSlit = (np.array([slit.x0, slit.x1])*np.cos(rotAngle) -
                                 np.array([  slit.y, slit.y])*np.sin(rotAngle))
                        ySlit = (np.array([slit.x0, slit.x1])*np.sin(rotAngle) +
                                 np.array([  slit.y, slit.y])*np.cos(rotAngle))
                    elif side == 'L':
                        xSlit = (np.array([-slit.x0, -slit.x1])*np.cos(rotAngle) -
                                 np.array  ([-slit.y, -slit.y])*np.sin(rotAngle))
                        ySlit = (np.array([-slit.x0, -slit.x1])*np.sin(rotAngle) +
                                 np.array  ([-slit.y, -slit.y])*np.cos(rotAngle))
                    # Assign coordinates
                    posX, posY = np.average([xSlit]), np.average([ySlit])
                    RAslit, Decslit = RAgal_deg-posX/3600., Decgal_deg-posY/3600.
                    RASlit_str = Coords2String(RAslit, hh=True)
                    DecSlit_str = Coords2String(Decslit)

                    if slit.type == 'SKiMS':
                        nameSlit_str = 'SS_'+str(int(self.mask_PA))+'_'+str(counter)
                    elif slit.type == 'Sky':
                        nameSlit_str = 'SS_sky_'+str(int(self.mask_PA))+'_'+str(counter)

                    slitPA_str = str(round(self.mask_PA+5., 2))
                    len1_str, len2_str = str(round(slit.length/2., 2)), str(round(slit.length/2.,2))

                    name_str = nameSlit_str
                    ra_str = RASlit_str
                    dec_str = DecSlit_str
                    width_str = str(float(len1_str) + float(len2_str)) + '\"'
                    height_str = '1.0\"'
                    angle_str = str(float(slitPA_str) + 90.)
                    f.write('box(' + ', '.join([ra_str, dec_str, width_str, height_str, angle_str]) +
                            ') # text={' + name_str + '}\n')
                    counter += 1
        
        
    def createOutputDSIM(self, pathOutput='./catSS.txt', save_slit_length=False):
        # It creates a catalog object for DSIM
        # All the slits are placed/rotated respect to their center and the 2 semi-  lengths are provided
        RAgal_deg, Decgal_deg = convAngCoord(self.gal_RA)[4]*15., convAngCoord(self.gal_Dec)[4]
        counter = 1
        listStrings = []
        for ii in self.listSlits:
            for side in ['L', 'R']:
                # Rotate with mask
                rotAngle = np.radians(90.-self.mask_PA)
                # RIGHT HAND SLIT
                if side == 'R':
                    xSlit = np.array([ii.x0, ii.x1])*np.cos(rotAngle) - np.array([  ii.y, ii.y])*np.sin(rotAngle)
                    ySlit = np.array([ii.x0, ii.x1])*np.sin(rotAngle) + np.array([  ii.y, ii.y])*np.cos(rotAngle)
                elif side == 'L':
                    xSlit = np.array([-ii.x0, -ii.x1])*np.cos(rotAngle) - np.array  ([-ii.y, -ii.y])*np.sin(rotAngle)
                    ySlit = np.array([-ii.x0, -ii.x1])*np.sin(rotAngle) + np.array  ([-ii.y, -ii.y])*np.cos(rotAngle)
                # Assign coordinates
                posX, posY = np.average([xSlit]), np.average([ySlit])
                RAslit, Decslit = RAgal_deg-posX/3600., Decgal_deg-posY/3600.
                
                if ii.type == 'SKiMS':
                    nameSlit_str = 'SS_'+str(int(self.mask_PA))+'_'+str(counter)
                    if ii.SB > 500:
                        mag_str = str(20)   # Bogus number, just to prevent issues with DSIM
                    else:
                        mag_str = str(round(ii.SB, 2))
                    pcode_str = str(int(1001-counter))
                    sample_str = '1'
                    select_str = '1' #Pre-selected
                elif ii.type == 'Sky':
                    nameSlit_str = 'SS_sky_'+str(int(self.mask_PA))+'_'+str(counter)
                    mag_str = '0'
                    pcode_str = str(int(300-counter))
                    sample_str = '3'
                    select_str = '1' #Pre-selected

                RASlit_str = Coords2String(RAslit, hh=True)
                DecSlit_str = Coords2String(Decslit)
                equinox_str = '2000'
                passband_str = 'R'  #Not true
                slitPA_str = str(round(self.mask_PA+5., 2))
                len1_str, len2_str = str(round(ii.length/2., 2)), str(round(ii.length/2.,2))
                slitWidth_str = '1'

                if save_slit_length:
                    listStrings.append(nameSlit_str+'\t'+RASlit_str+'\t'+DecSlit_str+'\t' +
                                       equinox_str+'\t'+mag_str+'\t'+passband_str+'\t'+pcode_str +
                                       '\t'+  sample_str+'\t'+select_str+'\t'+slitPA_str+'\t' +
                                       len1_str+'\t'+len2_str+  '\t'+slitWidth_str+'\n')
                else:
                    listStrings.append(nameSlit_str+'\t'+RASlit_str+'\t'+DecSlit_str+'\t'+
                                       equinox_str+'\t'+mag_str+'\t'+passband_str+'\t'+pcode_str+
                                       '\t'+  sample_str+'\t'+select_str+'\t'+slitPA_str+'\n')
                counter += 1
        outputTable =  np.column_stack(listStrings)
        np.savetxt(pathOutput, outputTable, delimiter="", fmt="%s")
        return True


class Slit:

    def __init__(self, pos0, length, coneAngle, y='random', xRange=[-498.,498.], yRange=[-146.,146.], target='SKiMS', SB=0.):
        self.x0 = pos0
        self.length = length
        self.x1 = self.x0 + self.length
        self.SB = SB
        self.type = target
        if y == 'random':
            if coneAngle < 180:
                elevation = (self.x0+self.length/2.)*np.tan((coneAngle/2.)*np.pi/180.)
            else:
                elevation = yRange[1]
            flag = True
            while flag:
                y = (np.random.rand()*2.*elevation)-elevation
                if (yRange[0] < y < yRange[1]) and (withinMask(self.x1, y)):
                    self.y = y
                    flag = False
        else:
            self.y = y
        self.centralCoords = [np.mean([self.x0, self.x1]), self.y]

    
def findBestMask(gal_Reff, gal_ba, gal_RA, gal_Dec, gal_PA, mask_PA, coneAngle,
                 iterations=100., maxDist=np.inf, separationSlits=0.4, minWidthSlits=3):
    t1 = time.time()
    for ii in np.arange(iterations):
        tmpObj = Mask(gal_Reff, gal_ba, gal_RA, gal_Dec, gal_PA, mask_PA,
                      coneAngle, separationSlits, minWidthSlits)
        tmpObj.createSlits()
        tmpDist = tmpObj.getMaxEmptySpace()
        if tmpDist < maxDist:
            maxDist = tmpDist
            ## Adding Sky Slits
            tmpObj.createSkySlits()
            m = tmpObj
    t2 = time.time()
    print str((t2 - t1) / iterations), 's per iteration'
    return m, maxDist



    
    










