#### SuperSKiMS_mask_design
###### Python code to retrieve an optimal mask design for SKiMS on DEIMOS. 

###### Latest version: Deimos_SKiMS_slit.V0.5.py - 28th September 2015

Copyright: Nicola Pastorello (nicola.pastorello@gmail.com, 2015)


   To run: 
   1) Set initial galaxy parameters in the header section of 'Deimos_SKiMS_slit.V0.5.py'
   2) run: "python Deimos_SKiMS_slit.V0.5.py"


The code finds the optimal SuperSKiMS slit distribution, given the galaxy surface 
brightness profile and the orientation of the DEIMOS mask. 
The output is a '.cat' file containing the SuperSKiMS and sky slits that is meant to be 
merged with a standard DEIMOS catalog. Such a final catalog has to be used in 
'dsimulator' in order to create the DEIMOS mask design file that will be submitted to 
the DEIMOS milling website. 
