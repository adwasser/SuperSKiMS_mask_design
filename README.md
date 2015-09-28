#### SuperSKiMS_mask_design
###### Python code to retrieve an optimal mask design for SKiMS on DEIMOS. 

###### Latest version: Deimos_SKiMS_slit.V0.5.py - 28th September 2015

Copyright: Nicola Pastorello (nicola.pastorello@gmail.com, 2015)


_To run: 
1. Set initial galaxy parameters in the header section of 'Deimos_SKiMS_slit.V0.5.py'
2. run: "python Deimos_SKiMS_slit.V0.5.py"_


The code finds the optimal SuperSKiMS slit distribution, given the galaxy surface 
brightness profile and the orientation of the DEIMOS mask. 
The output is a '.cat' file containing the SuperSKiMS and sky slits that is meant to be 
merged with a standard DEIMOS catalog. Such a final catalog has to be used in 
'dsimulator' in order to create the DEIMOS mask design file that will be submitted to 
the DEIMOS milling website. 

The slit distribution is symmetrical. The best configuration is found via a 
MonteCarlo approach. This configuration is the one that minimize the largest 
empty contiguous area within the mask. 
The length of the slits depends on the galaxy surface brightness **(still to be properly implemented**). 
This can be:
1. provided as an input txt file
2. provided as a Sersic profile, given the profile's parameters
3. assumed as a de Vaucouleurs profile

Still to fix:
-[] decide which band to use to photometric profile in.
-[] retrieving galaxy parameters from dictionary

