import os, sys, numpy as np


obs = sys.argv[1]

obsdir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/'

event = os.path.join(obsdir,obs)
reg = os.path.join(event,'event_defcl/excl.reg')
event = os.path.join(event,'event_defcl/im3to30keV.fits')

os.system('ds9 {} -cmap Hsv -smooth yes -scale limits 0 10 -smooth radius 5 -region {}'.format(event,reg))
#os.system('ds9 -region {}'.format(reg))




