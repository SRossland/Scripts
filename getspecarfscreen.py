#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./getspecarfscreen.py dir obsid regcore subdir outdir det rmf

import sys
import os
import string

dir=sys.argv[1]
obsid=sys.argv[2]
regcore=sys.argv[3]
outdir=sys.argv[5]
det=sys.argv[6]
subdir=sys.argv[4]
dir=dir+"/"+obsid+"/"
arf=sys.argv[7]
rmfarf=" runmkarf="+arf+" runmkrmf=yes "

i=0
params = ["SUN", "NOSUN", "OCC", "NOOCC"]
while i < len(params):
    os.system("nuproducts "+ \
              "infile="+dir+"event_cl/"+subdir+"/nu"+obsid+det+"01_fullevts_"+params[i]+".fits "+ \
              " srcregionfile="+dir+"event_cl/"+regcore+".reg "+ \
              " indir="+dir+"event_cl "+ \
              " outdir="+dir+"event_cl/"+outdir+rmfarf+ \
              " bkgextract=no "+ \
              " lcfile=NONE "+ \
              " instrument=FPM"+det+" steminputs=nu"+obsid+ \
              " stemout="+regcore+" boxsize=20 "+ \
              " clobber=yes")

    print "nuproducts "+ \
            "infile="+dir+"event_cl/"+subdir+"/nu"+obsid+det+"01_fullevts_"+params[i]+".fits "+ \
            " srcregionfile="+dir+"event_cl/"+regcore+".reg "+ \
            " indir="+dir+"event_cl "+ \
            " outdir="+dir+"event_cl/"+outdir+rmfarf+ \
            " bkgextract=no "+ \
            " lcfile=NONE "+ \
            " instrument=FPM"+det+" steminputs=nu"+obsid+ \
            " stemout="+regcore+" boxsize=20 "+ \
            " clobber=yes"
              
    os.system("mv "+dir+"event_cl/"+outdir+"/"+regcore+"_sr.rmf "+ \
            dir+"event_cl/"+outdir+"/"+regcore+"_sr_orig.rmf")
    os.system("cmprmf "+dir+"event_cl/"+outdir+"/"+regcore+"_sr_orig.rmf "+ \
            dir+"event_cl/"+outdir+"/"+regcore+"_sr.rmf 1e-6")
    os.system("rm -f "+dir+"event_cl/"+outdir+"/"+regcore+"_sr_g30.pha")
    os.system("grppha "+dir+"event_cl/"+outdir+"/"+regcore+"_sr.pha "+ \
            dir+"event_cl/"+outdir+"/"+regcore+ \
            "_sr_g30.pha 'chkey RESPFILE "+dir+"event_cl/"+outdir+"/"+ \
            regcore+"_sr.rmf & group min 30 & exit'")
    os.system("mv "+dir+"event_cl/"+outdir+"/"+regcore+"_im.gif "+ \
            dir+"event_cl/"+outdir+"/"+regcore+"_im_screen"+params[i]+".gif")
    os.system("mv "+dir+"event_cl/"+outdir+"/"+regcore+"_ph.gif "+ \
            dir+"event_cl/"+outdir+"/"+regcore+"_ph_screen"+params[i]+".gif")
    os.system("mv "+dir+"event_cl/"+outdir+"/"+regcore+"_sr_orig.rmf "+ \
            dir+"event_cl/"+outdir+"/"+regcore+"_sr_orig_screen"+params[i]+".rmf")
    os.system("mv "+dir+"event_cl/"+outdir+"/"+regcore+"_sk.img "+ \
            dir+"event_cl/"+outdir+"/"+regcore+"_sk_screen"+params[i]+".img")
    os.system("mv "+dir+"event_cl/"+outdir+"/"+regcore+"_sr_g30.pha "+ \
            dir+"event_cl/"+outdir+"/"+regcore+"_sr_g30_screen"+params[i]+".pha")
    os.system("mv "+dir+"event_cl/"+outdir+"/"+regcore+"_sr.pha "+ \
            dir+"event_cl/"+outdir+"/"+regcore+"_sr_screen"+params[i]+".pha")
    os.system("mv "+dir+"event_cl/"+outdir+"/"+regcore+"_sr.rmf "+ \
            dir+"event_cl/"+outdir+"/"+regcore+"_sr_g30_screen"+params[i]+".rmf")
    if sys.argv[7]=='yes' or 'Yes' or 'YES' or 'YeS' or 'yEs' or 'yES' or 'yeS':
        os.system("mv "+dir+"event_cl/"+outdir+"/"+regcore+"_sr.arf "+ \
                  dir+"event_cl/"+outdir+"/"+regcore+"_sr_screen"+params[i]+".arf")
    i += 1
