#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./sepsunocc2.py dir obsid det

""" This seperates the 01 and 02 _cl.evt files for both det A and B into an alternate directory 
    (event_sep_cl) with sub directories for both occulted and non-occulted times. Inside each of those 
    (if both exist) the times are further separated based on the SUNSHINE flag in the event file header.
    The output is the GTI and full event file for both A and B in each instance, do note though, that while
    it looks for these, if the telescope was never occulted you will not have that file folder or if it was 
    never in the sun and so on"""

import sys
import os
import string

dir=sys.argv[1]
obsid=sys.argv[2]
det=sys.argv[3]
dir=dir+"/"+obsid+"/"


#for i in range(2):
#print('This is currently set to only process the 02s')
#blah = raw_input('Hit enter if that is ok with you')
for i in [0]: # to force only 01 processing
#for i in [1]: # to force only 02 processing
    
    if i==0:
        otfile="NOOCC"
        k=str(1)
    elif i==1 and os.path.isfile(dir+"event_cl/nu"+obsid+det+"02_cl.evt"):
        otfile="OCC"
        k=str(2)
    for j in range(1): #currently only have to redo the nosun seps
        if j==0:
            fi="NOSUN"
        elif j==1:
            fi="SUN"

        print "nuscreen "+ \
            " infile="+dir+"event_cl/nu"+obsid+det+"0"+k+"_cl.evt "+ \
            " hkfile="+dir+"event_cl/nu"+obsid+det+"_fpm.hk "+ \
            " mkffile="+dir+"event_cl/nu"+obsid+det+".mkf "+ \
            " outdir="+dir+"event_sep_cl/"+otfile+ \
            " gtiscreen=yes "+ \
            " evtscreen=no "+ \
            " gtiexpr=SUNSHINE=="+str(j)+" "+ \
            " outfile=nu"+obsid+det+"0"+k+"_fullevts_"+fi+".fits "+ \
            " clobber=yes"


        os.system("nuscreen "+ \
              " infile="+dir+"event_cl/nu"+obsid+det+"0"+k+"_cl.evt "+ \
              " hkfile="+dir+"event_cl/nu"+obsid+det+"_fpm.hk "+ \
              " mkffile="+dir+"event_cl/nu"+obsid+det+".mkf "+ \
              " outdir="+dir+"event_sep_cl/"+otfile+ \
              " gtiscreen=yes "+ \
              " evtscreen=no "+ \
              " gtiexpr=SUNSHINE=="+str(j)+" "+ \
              " outfile=nu"+obsid+det+"0"+k+"_fullevts_"+fi+".fits "+ \
              " clobber=yes")


        os.system("mv "+dir+"event_sep_cl/"+otfile+"/nu"+obsid+det+"0"+k+"01_gti.fits "+ \
              dir+"event_sep_cl/"+otfile+"/nu"+obsid+det+"0"+k+"01_gti_"+fi+".fits")
