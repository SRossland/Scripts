#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./elevcorr.py dir obsid ELV

# Created by:  Steven P. Rossland, 2018

""" The given 02 files from heasarc are filtered by an elev<=3 parameter.  This allows some light leak in, so to correct for this we rescreen the evt file for an ELV<=0.  This was found to be acceptable empiraclly by looking at a bright cluster (80202014006). Note, the guide mentions the flag is ELEV, but the file has a colmun header of ELV.  """

import sys
import os
import string

dir=sys.argv[1]
obsid=sys.argv[2][0:11]
ELEV=sys.argv[3]
dir=dir+"/"+obsid+"/"


for det in ['A','B']:
#for det in ['B']:    

	if os.path.isfile(dir+"event_cl/nu"+obsid+det+"02_cl.evt"):
        	otfile="test"
	elif det == 'A':
		print("No A event found")
		continue
	else:
		print("Can not find file 02 for "+obsid)
		sys.exit()

	print "nuscreen "+ \
		" infile="+dir+"event_cl/nu"+obsid+det+"02_cl.evt\n"+ \
        	" gtiscreen=yes\n"+ \
        	" evtscreen=no\n"+ \
        	" gtiexpr='ELV<="+ELEV+"'\n"+ \
        	" gradeexpr=NONE\n" + \
        	" statusexpr=NONE\n" + \
        	" outdir="+dir+otfile+"\n"+ \
        	" obsmode=OCCULTATION\n" + \
        	" hkfile="+dir+"event_cl/nu"+obsid+det+"_fpm.hk\n"+ \
        	" mkffile="+dir+"event_cl/nu"+obsid+det+".mkf\n"+ \
        	" outfile=nu"+obsid+det+"02_cl.evt\n" + \
        	" clobber=yes"

	os.system("nuscreen "+ \
        	" infile="+dir+"event_cl/nu"+obsid+det+"02_cl.evt "+ \
		" gtiscreen=yes"+ \
		" evtscreen=no"+ \
		" gtiexpr='ELV<="+ELEV+"'"+ \
		" gradeexpr=NONE" + \
		" statusexpr=NONE" + \
		" outdir="+dir+otfile+ \
		" obsmode=OCCULTATION" + \
        	" hkfile="+dir+"event_cl/nu"+obsid+det+"_fpm.hk "+ \
        	" mkffile="+dir+"event_cl/nu"+obsid+det+".mkf "+ \
#	" usrgtifile=NONE" + \
		" outfile=nu"+obsid+det+"02_cl.evt" + \
		" clobber=yes")
	
	os.system("rm "+dir+"event_cl/nu"+obsid+det+"02_gti.fits")
	os.system("mv "+dir+otfile+"/*.fits "+dir+"event_cl/nu"+obsid+det+"02_gti.fits")
	os.system("rm "+dir+"event_cl/nu"+obsid+det+"0202_gti.fits")
	os.system("mv "+dir+otfile+"/*.evt "+dir+"event_cl/")
	with open('/uufs/astro.utah.edu/common/home/u1019304/temp/batchcomplete.txt','a+') as O:
                        O.write(obsid+' '+det+'\n')
#	os.system("rm -f -r "+dir+"/test")
