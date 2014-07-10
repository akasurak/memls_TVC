README TO THE IAP SNOW EMISSION MODEL MEMLS, Version 3, by Christian Matzler, IAP, April 2007, incl. rev. May 2012, July 2014

For further information, see the MEMLS DOCUMENTATION, and the MEMLS reports and publications:

MŠtzler C., "Notes on microwave radiation from snow samples and from a layered snowpack", IAP Research Report 96-9, University of Bern, Switzerland, 1996, updated 2004.

A. Wiesmann, C. MŠtzler, "Microwave emission model of layered snowpacks", Remote Sensing of Environment, Vol. 70, No. 3, pp. 307-316 (1999).

C. MŠtzler and A. Wiesmann, "Extension of the Microwave Emission Model of Layered Snowpacks to Coarse-Grained Snow", Remote Sensing of Environment, Vol. 70, No. 3, pp. 317-325 (1999).

A. Wiesmann, C. Fierz and C. MŠtzler, "Simulation of microwave emission from physically modeled snowpacks", Annals of Glaciology, Vol. 31, pp. 397-405 (2000).

C. MŠtzler and A. Wiesmann, Microwave Emission Model of Layered Snowpacks, Documentation for MEMLS, Version 3, Institute of Applied Physics - University of Bern, Sidlerstrasse 5, 3012 Bern, Switzerland (2007), with minor revisions May (2012).

Pdf files of these articles can be obtained from the IAP Web site, http://www.iap.unibe.ch/

A revision of MEMLS became adequate when it was discovered (in 2006) that certain errors occurred in the option of the Improved Born Approximation. These errors should be eliminated in the present Version 3.

Version 3, described here, uses updated formulas for the dielectric constants of ice and water, including now also slightly saline ice, to account for salt-contaminated snow on sea ice (Chapter 5.3 in MŠtzler et al. 2006). For pure ice the update is a very small change in the imaginary part of the dielectric constant, since the formulas used in Versions 1 and 2 were mainly confirmed with new measurements. However, the salinity also affects the imaginary part even at the level of a few ppm; the increase is significant if the temperature approaches the melting point of ice. The applicable salinity range covers 0 to 0.1 ppt, but may still work up to 1 ppt. The dielectric constant of water is used in case of wet snow. Here we use the two-Debye relaxation formula of Liebe et al. (1991). Although more accurate formulations exist (Ellison, Chapter 5.2, in MŠtzler et al. 2006; Ellison, 2007), the simplicity of Liebe's expression is an advantage, and the accuracy is sufficient in view of the knowledge on the liquid-water content in wet snow.In addition the distinction on scattering between different types of snow can now be specified by a graintype parameter with graintype=1 for natural dry snow as observed by MŠtzler (1996a), graintype=2 for spherical scatterers, and graintype=3 for thin spherical shells (see routine descriptions). The current option is graintype=1.
Chapter 5 in C. MŠtzler (Ed.), P.W. Rosenkranz, A. Battaglia and J.P. Wigneron (Co-Eds.), "Thermal Microwave Radiation - Applications for Remote Sensing", IET Electromagnetic Waves Series 52, London, UK (2006)..

The option of 'graintype', is set to 1 for 'natural snow'. Other options are 2 for spheres and 3 for thin spherical shells. This feature was introduced by Ludovic Brucker, and at IAP we implemented it consistently in this MEMLS version. It meant that the real part of the dielectric constant of dry snow had to be concentrated into a routine together with the K^2 factor used in the Born approximation (see vK2eps.m).
 
Furthermore the computation can provide both the emissivity and the brightness temperature. This may be relevant due to the fact that snowpacks are often not isothermal. 

In 2013-2014 MEMLS was extended by an active version to compute backscattering.

In the enclosed archive you find Version 3 (including the active part as well) of MATLAB routines needed to run MEMLS. 

The main routines memlsmain (and amemlsmain) 
are the basic programs to compute the brightness temperature of a snowpack (and including the backscatter coefficient). The main input parameters are contained in the input file (test-input files are included) where each line contains the data of one snow layer, with the top line representing the bottom layer just above the ground, and the bottom line representing the top snow layer!

Each line (layer) contains the following 7 input data:
Layer no., 
Layer temperature (K),
volumetric liquid-water content (Vol. fraction range: 0-1)
snow density (kg/m3)
Layer thickness (cm)
Snow salinity (parts per thousand)
Exponential correlation length (mm)

lmain2 and fmain2 are similar to memlsmain, but with computing
the spectral (lmain2) or angular (fmain2) dependence.
These routines provide two plots each, and a datafile containing the numeric data.
these are the emissivities eh,v at h and v pol., the emitted brightness temperatures at h and v pol (i,e for Tsky=0) and the effective temperature Teffh,v=Teh,v/eh,v (K).

All other routines are called by the main programs.
Some routine names and variable names have changed from earlier versions. But the structure is still very similar.

try:
memlsmain(21, 50, 0, 0, 'test1layer.txt', 0, 273, 11)
amemlsmain(21, 50, 0, 0, 0, 0, 'test1layer.txt', 0, 273, 12, 0.1, 0.1)
lmain2(10, 100, 50, 0, 0, 'test1layer.txt', 273, 11)
fmain2(35, 0, 80, 0, 0, 'test1layer.txt', 273, 11)

All subroutines are commented. For questions contact matzler@iap.unibe.ch

