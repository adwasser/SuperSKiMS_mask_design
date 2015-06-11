## DEIMOS mask - SKiMS fill
#v_0.1 - Assuming constant Surface Brightness
#v_0.2 - Sersic SB 
#v_0.3 - reboot

from Deimos_SKiMS_slit__def__ import *

###########################
## System-wide variables ##
###########################

__builtins__.verbosity = False	#If True, keeps the verbosity level high
__builtins__.graphicOutput = False	#If True, shows the final plot at the end of every ScoBen iteration

__builtins__.separationSlits = 0.5	#Standard separation between slits in arcsec
__builtins__.minWidthSlits = 3.		#Minimum slit width in arcsec 
__builtins__.coneAngle = 45.		#Angle of cone containing the slits

__builtins__.gal_Reff = 47.86		#Effective radius of galaxy


###############
## Functions ##
###############

import multiprocessing as mp
from multiprocessing import Process
import platform

def status(proc):
  """Check for processes status"""
  if proc.is_alive==True:
    return 'alive'
  elif proc.is_alive==False:
    return 'dead'
  else:
    return proc.is_alive()

def worker(inputqueuenumber):
  """Running the different 100 iterations"""
  mask, maxDist = findBestMask(iterations = 100.,
                realProfilePath='photprofs_Rband.txt') #It uses the real profile to define the slit length
  fileOut = open('MaskObj'+str(inputqueuenumber)+'.dat', 'wb')
  try:
    del mask.ax	#Issues with pickling of object
  except:
    dummy = True
  pickle.dump([mask, maxDist], fileOut, pickle.HIGHEST_PROTOCOL)
  fileOut.close()

##########
## MAIN ##
##########

if platform.system() == 'Linux':
  #Let's parallel!
#  nproc = mp.cpu_count()
  nproc = 5
#
  procs = []    # processes container
  # Start the worker processes
  for ii in numpy.arange(5):
    print "Process loop ", ii
    procs.append(mp.Process(target=worker, args=([ii+1])))
#
  for ii in procs:
    ii.start()
#
  for ii in procs:
    print "Process ", ii," @ " , ii.pid, " is ", status(ii)
# Wait processes to finish
  while ((os.path.exists('./MaskObj1.dat')) and (os.path.exists('./MaskObj2.dat')) and (os.path.exists('./MaskObj3.dat')) and
  (os.path.exists('./MaskObj4.dat')) and (os.path.exists('./MaskObj5.dat'))) != True:
    time.sleep(10) # loose 10 seconds
#



###
# Find best configuration # 
import glob
maxDist, bestMask = 1e10, ''
for ii in glob.glob('MaskObj*dat'):
  fileIn = open(ii, 'rb')
  mask, maxDistTmp = pickle.load(fileIn)
  fileIn.close()
  #
  if maxDistTmp < maxDist:
    maxDist = maxDistTmp
    bestMask = ii

os.system('cp '+bestMask+' BestMask.dat')

print "DONE"


'''

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

## Multiple masks plot
mask.initializePlot(figSize=(10,10))
mask.plotMask(PAangle=0.)
mask.plotSlits(PAangle=0.)
mask.plotBoundaries(PAangle=0.)

mask.plotMask(PAangle=90.)
mask.plotSlits(PAangle=90.)
mask.plotBoundaries(PAangle=90.)

mask.plotMask(PAangle=135.)
mask.plotSlits(PAangle=135.)
mask.plotBoundaries(PAangle=135.)

mask.plotMask(PAangle=45.)
mask.plotSlits(PAangle=45.)
mask.plotBoundaries(PAangle=45.)

mask.ax.set_ylim([-550,550])

savefig('MultipleSSKiMSmask.pdf', bbox_inches='tight')

## Single mask plot with de Vaucouleurs profile
mask.initializePlot(figSize=(10,5))
mask.plotMask(PAangle=0.)
mask.plotSlits(PAangle=0.)
mask.plotBoundaries(PAangle=0.)
mask.plotBins()
mask.ax.set_ylim([-500,180])


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
