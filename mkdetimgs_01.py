#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./mkdetimgs_01.py dir obsid elowlist ehighlist 01or02 sunocc sunno

# This code should replace mkdetimgs.py and mkdetimags_02.py by allowing the user to 
# pass sci or occ call and sun/nosun call to make the det images as desired.

import sys
import os
import string

try:
    fobs = open(sys.argv[1]+'/'+sys.argv[2]+'/'+sys.argv[2]+'.dat', 'r')
except IOError:
    obsids=[sys.argv[2]]
    if len(sys.argv) == 6:
        print('len 6 call')
        dir=sys.argv[1]+"/"+obsids[0]+"/event_cl/"
    else:
        print('There are ',len(sys.argv),' arguments')
        dir=sys.argv[1]+"/"+obsids[0]+"/event_sep_cl/"+sys.argv[6]+"/"
else:
    obsids = fobs.readlines()
    for i in range(len(obsids)):
        obsids[i]=obsids[i].rstrip()
    fobs.close()
    dir=sys.argv[1]+'/'+sys.argv[2]+'/'

slow=sys.argv[3].split(',')
shigh=sys.argv[4].split(',')
usrtype = sys.argv[5]
elow=[int(i) for i in slow]
ehigh=[int(i) for i in shigh]
clow=[str(int((i-1.6)/0.04+1)) for i in elow]
chigh=[str(int((i-1.6)/0.04)) for i in ehigh]

for det in ["A","B"]:
#for det in ["B"]:
    xsel=open(dir+"xsel.xco","w")
    xsel.write("session1\n")
    iobs=0
    for obsid in obsids:
        edir=sys.argv[1]+'/'+obsid+'/'
#	edir = './'+obsid
        if len(sys.argv) == 6:
            if os.path.isfile(edir+"event_cl/nu"+obsid+det+usrtype+"_cl.evt"):
                xsel.write("read event "+edir+"event_cl/nu"+obsid+det+usrtype+"_cl.evt\n")
            else:
                xsel.write("exit\n")
                xsel.write("no\n")
                xsel.close()
#            else:
#                if os.path.isfile(edir+"event_sep_cl/"+sys.argv[6]+"/nu"+obsid+det+usrtype+"_fullevts_"+sys.argv[7]+".fits") == True:
#                    xsel.write("read event "+obsid+"/event_sep_cl/"+sys.argv[6]+"/nu"+obsid+det+usrtype+"_fullevts_"+sys.argv[7]+".fits\n")
#                else:
#                    xsel.write("exit\n")
#                    xsel.write("no\n")
#                    xsel.close()

            if iobs == 0:
                xsel.write("./\n")
                xsel.write("yes\n")
            iobs=iobs+1

    for i in range(len(slow)):
        try:
            fobs=open(dir+"im"+det+slow[i]+"to"+shigh[i]+"keVDET.fits",'r')
        except IOError:
            blah=1
        else:
            os.system("rm -f -r "+dir+"im"+det+slow[i]+"to"+ \
                      shigh[i]+"keV.fits")
        xsel.write('filter pha_cutoff '+clow[i]+' '+chigh[i]+'\n')
        xsel.write('set xyname det1x det1y\n')
        xsel.write("extract image\n")
        xsel.write("save image\n")
        if len(sys.argv) == 6:
            xsel.write(dir+"im"+det+slow[i]+"to"+shigh[i]+"keVDET.fits\n")
        else:
            xsel.write(dir+"im"+det+slow[i]+"to"+shigh[i]+"keV"+sys.argv[7]+"DET.fits\n")
        xsel.write("clear\n")
        xsel.write("pha_cutoff\n")

    xsel.write("exit\n")
    xsel.write("no\n")
    xsel.close()
    os.system("xselect @"+dir+"xsel.xco")
    os.system("rm -r -f "+dir+"xsel.xco")

