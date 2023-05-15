import os, sys, numpy as np
from astropy.io import fits

obslist = sys.argv[1]

obs = np.genfromtxt('/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs/{}'.format(obslist), dtype=int)



obs_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/'

Aexp,Bexp = 0,0

for ob in obs:
    eventA = '{}{}/event_cl/nu{}A02_cl.evt'.format(obs_dir,ob,ob)
    eventB = eventA.replace('A02_cl','B02_cl')
    if os.path.isfile(eventA):
        with fits.open(eventA) as hdul:
            Aexp += hdul[0].header['EXPOSURE']
    if os.path.isfile(eventB):
        with fits.open(eventB) as hdul:
            Bexp += hdul[0].header['EXPOSURE']

print('FPMA Exposure: {}'.format(Aexp))
print('FPMB Exposure: {}'.format(Bexp))





