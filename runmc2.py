#!/uufs/astro.utah.edu/common/home/u1019304/VENV3.7.9/bin/python3

# Syntax: runmc2.py path/to/directory

import os, sys, numpy as np
import time

work_dir = sys.argv[1]

dir_list = os.listdir(work_dir)
os.chdir(work_dir)
# need list of energy values
elow, ehigh = [],[]

for i,li in enumerate(dir_list):
    if 'longformat' in li:
        if 'A_nosun' in li:
            elow.append(li.split('_')[2])
            ehigh.append(li.split('_')[3])

if (len(elow) == 0) or (len(elow) != len(ehigh)):
	print('You suck again!')
	sys.exit()

energy = sorted(zip(elow,ehigh), key = lambda x: float(x[0]))
elow, ehigh = zip(*energy)

elow = list(elow)
ehigh = list(ehigh)

tele = ['A','B']
for det in ['A','B']:
    detmapfile = '/uufs/astro.utah.edu/common/home/u1019304/temp/detmap'+det+'.fits'
    gradfile = '/uufs/astro.utah.edu/common/home/u1019304/my_idl/nuskybgd/auxil/det'+det+'_det1.img'
    for i,lo in enumerate(elow):
        start_time = time.time()
        datafile = work_dir+'/fits_file/Data'+det+'_01_nosun_'+lo+'_'+ehigh[i]+'keV.fits'
        binfile = datafile.replace('Data','BINS')
        expfile = datafile.replace('Data','Exp')
        normfile = work_dir+'/'+det+'_nosun_'+lo+'_'+ehigh[i]+'_keVparams_longformat.txt'
        os.system('mcmcNormspy3.py '+datafile+' '+binfile+' '+detmapfile+' '+expfile+' '+gradfile+' '+normfile)
        print('----------%s seconds-----------' %(time.time()-start_time))

 
