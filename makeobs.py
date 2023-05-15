#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# This script takes a list and makes files to pass to my slurm script

import numpy as np, sys, os

dirl = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs/'
list = np.loadtxt(dirl+sys.argv[1])
dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs/OBS_ID'

os.system('rm -f '+dir+'/*.obs')

for i in range(len(list)):
    if os.path.isfile(dir+'/'+str(int(list[i]))+'.obs') == False:
	f = open(dir+'/'+str(int(list[i]))+'.obs','a')
    else: 
	continue
