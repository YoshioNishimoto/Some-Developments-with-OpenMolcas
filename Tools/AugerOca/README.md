This is the readme file for auger-OCA

Contact person: Bruno Tenorio
b.nunes.bruno@gmail.com

# USE "pip install auger-oca" TO INSTALL the auger-oca packages
#=====================================

 OCA stands for One-Center Approximation. It is a simple way to get Auger intensities bypassing the
 calculation of continuum wave function.
 
 Auger-OCA is a collection of python scripts (programs) to post-process a collection
 of 2-el reduced transition density matrices (also called Auger densities) generated 
 by the RASSI module of OpenMolcas. It also requires the \$project.rassi.h5 file in the background.

 The keyword TDYS must be used together with DYSON in your RASSI input calculation (see example below
 for C 1s edge).

#\&RASSI 
# Nrof JobIphs =  5 all 
# CIprint 
# Ejob 
# TDYS 
# 1 
# C 1s 
# Dyson 

 Check auger_test/ folder for complete examples. 

It simply runs as:

 for RAES (resonant Auger):
 
 \$python3 \$MOLCAS/Tools/auger_oca/auger_main.py -d \$WorkDir --raes --spec \&

 If non-resonant Auger (AES) of singlet final states:

 \$python3 \$MOLCAS/Tools/auger_oca/auger_main.py -d \$WorkDir --aes --s --spec \&

 If non-resonant Auger (AES) of triplet final states:

 \$python3 \$MOLCAS/Tools/auger_oca/auger_main.py -d \$WorkDir --aes --t --spec \&

++ See also

 \$python3 \$MOLCAS/Tools/auger_oca/auger_main.py -h 

 for more info and usage. 

+Important: 1) auger_main.py reads a file \$project.rassi.h5 in the background. You must execute auger_main.py 
 from any directory where you have your \$project.rassi.h5 file (either \$CurrDir or \$WorkDir work).
 2) auge_main.py will create a directory called 'auger_outputs/' containing a collection of
 [Auger_OCA.r2TM_XXX_YYY_ZZZ.out,...] files. The binding energy and Auger intensities (plus other info)
 are writen on each of these files.
 3) '--spec' is optional. One can also use a bash script ('oca.spectrum.sh' on auger_test folder) to collect the Auger {energy,intensity}
 from those [Auger_OCA.r2TM_XXX_YYY_ZZZ.out,...] files -> in a file called spectrum.out.
 4) Again, don't forget to install the Auger OCA packages with 'pip install auger-oca'
  