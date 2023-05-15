#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

import os, sys, numpy as np



log = sys.argv[1]

obs = np.genfromtxt(log,dtype='|U11')

evt_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/'

for i,ob in enumerate(obs):
  for det in ['A']:
    os.system('sepsunocc2.py '+evt_dir+' '+ob+' '+det)
