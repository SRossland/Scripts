#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./sepsunocc.py dir obsid outdir det

import sys
import os
import string

dir=sys.argv[1]
obsid=sys.argv[2]
outdir=sys.argv[3]
det=sys.argv[4]
dir=dir+"/"+obsid+"/"


for i in range(4):
    if i==0:
        os.system("nuscreen "+ \
              " infile="+dir+"event_cl/nu"+obsid+det+"01_cl.evt "+ \
              " hkfile="+dir+"event_cl/nu"+obsid+det+"_fpm.hk "+ \
              " mkffile="+dir+"event_cl/nu"+obsid+det+".mkf "+ \
              " outdir="+dir+outdir+"/NOSUN" \
              " gtiscreen=yes "+ \
              " evtscreen=no "+ \
              " gtiexpr=SUNSHINE==0 "+ \
              " outfile=nu"+obsid+det+"01_fullevts_NOSUN.fits "+ \
              " clobber=yes")

        print "nuscreen "+ \
            " infile="+dir+"event_cl/nu"+obsid+det+"01_cl.evt "+ \
            " hkfile="+dir+"event_cl/nu"+obsid+det+"_fpm.hk "+ \
            " mkffile="+dir+"event_cl/nu"+obsid+det+".mkf "+ \
            " outdir="+dir+outdir+"/NOSUN" \
            " gtiscreen=yes "+ \
            " evtscreen=no "+ \
            " gtiexpr=SUNSHINE==0 "+ \
            " outfile=nu"+obsid+det+"01_fullevts_NOSUN.fits "+ \
            " clobber=yes"

        os.system("mv "+dir+outdir+"/NOSUN/nu"+obsid+det+"0101_gti.fits "+ \
              dir+outdir+"/NOSUN/nu"+obsid+det+"0101_gti_NOSUN.fits")

    elif i==1:
        os.system("nuscreen "+ \
                  " infile="+dir+"event_cl/nu"+obsid+det+"01_cl.evt "+ \
                  " hkfile="+dir+"event_cl/nu"+obsid+det+"_fpm.hk "+ \
                  " mkffile="+dir+"event_cl/nu"+obsid+det+".mkf "+ \
                  " outdir="+dir+outdir+"/SUN" \
                  " gtiscreen=yes "+ \
                  " evtscreen=no "+ \
                  " gtiexpr=SUNSHINE==1 "+ \
                  " outfile=nu"+obsid+det+"01_fullevts_SUN.fits "+ \
                  " clobber=yes")
            
        print "nuscreen "+ \
                " infile="+dir+"event_cl/nu"+obsid+det+"01_cl.evt "+ \
                " hkfile="+dir+"event_cl/nu"+obsid+det+"_fpm.hk "+ \
                " mkffile="+dir+"event_cl/nu"+obsid+det+".mkf "+ \
                " outdir="+dir+outdir+"SUN" \
                " gtiscreen=yes "+ \
                " evtscreen=no "+ \
                " gtiexpr=SUNSHINE==1 "+ \
                " outfile=nu"+obsid+det+"01_fullevts_SUN.fits "+ \
                " clobber=yes"

        os.system("mv "+dir+outdir+"/SUN/nu"+obsid+det+"0101_gti.fits "+ \
          dir+outdir+"/SUN/nu"+obsid+det+"0101_gti_SUN.fits")

    elif i==2:
        os.system("nuscreen "+ \
                  " infile="+dir+"event_cl/nu"+obsid+det+"01_cl.evt "+ \
                  " hkfile="+dir+"event_cl/nu"+obsid+det+"_fpm.hk "+ \
                  " mkffile="+dir+"event_cl/nu"+obsid+det+".mkf "+ \
                  " outdir="+dir+outdir+"/NOOCC" \
                  " gtiscreen=yes "+ \
                  " evtscreen=no "+ \
                  " gtiexpr=OCCULTED==0 "+ \
                  " outfile=nu"+obsid+det+"01_fullevts_NOOCC.fits "+ \
                  " clobber=yes")
            
        print "nuscreen "+ \
                " infile="+dir+"event_cl/nu"+obsid+det+"01_cl.evt "+ \
                " hkfile="+dir+"event_cl/nu"+obsid+det+"_fpm.hk "+ \
                " mkffile="+dir+"event_cl/nu"+obsid+det+".mkf "+ \
                " outdir="+dir+outdir+"/NOOCC" \
                " gtiscreen=yes "+ \
                " evtscreen=no "+ \
                " gtiexpr=OCCULTED==0 "+ \
                " outfile=nu"+obsid+det+"01_fullevts_NOOCC.fits "+ \
                " clobber=yes"
    
        os.system("mv "+dir+outdir+"/NOOCC/nu"+obsid+det+"0101_gti.fits "+ \
          dir+outdir+"/NOOCC/nu"+obsid+det+"0101_gti_NOOCC.fits")

    elif i==3:
        """This is a check to see if the file even exists and goes on if it doesn't"""
        if os.path.isfile(dir+"event_cl/nu"+obsid+det+"02_cl.evt"):
            os.system("nuscreen "+ \
                  " infile="+dir+"event_cl/nu"+obsid+det+"02_cl.evt "+ \
                  " hkfile="+dir+"event_cl/nu"+obsid+det+"_fpm.hk "+ \
                  " mkffile="+dir+"event_cl/nu"+obsid+det+".mkf "+ \
                  " outdir="+dir+outdir+"/OCC" \
                  " gtiscreen=yes "+ \
                  " evtscreen=no "+ \
                  " gtiexpr=OCCULTED==1 "+ \
                  " outfile=nu"+obsid+det+"02_fullevts_OCC.fits "+ \
                  " clobber=yes")
    
            print "nuscreen "+ \
                " infile="+dir+"event_cl/nu"+obsid+det+"02_cl.evt "+ \
                " hkfile="+dir+"event_cl/nu"+obsid+det+"_fpm.hk "+ \
                " mkffile="+dir+"event_cl/nu"+obsid+det+".mkf "+ \
                " outdir="+dir+outdir+"/OCC" \
                " gtiscreen=yes "+ \
                " evtscreen=no "+ \
                " gtiexpr=OCCULTED==1 "+ \
                " outfile=nu"+obsid+det+"02_fullevts_OCC.fits "+ \
                " clobber=yes"

            os.system("mv "+dir+outdir+"/OCC/nu"+obsid+det+"0201_gti.fits "+ \
              dir+outdir+"/OCC/nu"+obsid+det+"0201_gti_OCC.fits")






