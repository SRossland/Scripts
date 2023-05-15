import os, sys, numpy as np
from astropy.io import fits


# This program is purely to stack 02 data files that have not been
# filtered yet. 

obsdir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/'
obslist = sys.argv[1]
elow = sys.argv[2]
ehigh = sys.argv[3]

mask_dir = '/uufs/astro.utah.edu/common/home/u1019304/temp'
with fits.open(mask_dir+'/fullmaskA_final.fits') as mAhdu:
    edgemaskA = mAhdu[0].data

with fits.open(mask_dir+'/fullmaskB_final.fits') as mBhdu:
    edgemaskB = mBhdu[0].data

minbin = (float(elow)-1.6)/0.04
maxbin = (float(ehigh)-1.6)/0.04

Adata, Bdata = np.zeros((360,360)), np.zeros((360,360))
zeroElevDirectory = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OneElevations'

#obs_orig = np.genfromtxt(obslist)
with open(obslist) as f:
    lines = f.readlines()

obs_orig = []
for i in range(2,len(lines)):
    li = lines[i].split()
    obs_orig.append(li[0])


for det in ['A','B']:
    for i,ob in enumerate(obs_orig):
        ob = int(ob)
        event = '{}/nu{}{}02_1ELV.fits'.format(zeroElevDirectory,ob,det)
        #event = '{}{}/event_defcl/nu{}{}02_cl.evt'.format(obsdir,ob,ob,det)
        if os.path.isfile(event):
            if det == 'A': em = np.copy(edgemaskA)
            else: em = np.copy(edgemaskB)
            em = em.astype('float64')
            with fits.open(event) as hdul:
                PI = hdul[1].data['PI']
                X = hdul[1].data['X']
                Y = hdul[1].data['Y']
                DET1_X = hdul[1].data['DET1X']
                DET1_Y = hdul[1].data['DET1Y']

            idx_good = (X > 0)*(Y > 0)

            for m in range(len(PI[idx_good])):
                elow_int = int(round((float(elow)-1.6)/0.04))
                ehigh_int = int(round((float(ehigh)-1.6)/0.04))
                if elow_int <= PI[idx_good][m] < ehigh_int:
                    if (em[DET1_Y[idx_good][m]-1,DET1_X[idx_good][m]-1] == 0):
                        continue
                    else:
                        if det == 'A':
                            Adata[DET1_Y[idx_good][m]-1,DET1_X[idx_good][m]-1] += 1
                        else:
                            Bdata[DET1_Y[idx_good][m]-1,DET1_X[idx_good][m]-1] += 1

fits.writeto('/uufs/astro.utah.edu/common/home/u1019304/temp/A_02_1ELV.fits',Adata)
fits.writeto('/uufs/astro.utah.edu/common/home/u1019304/temp/B_02_1ELV.fits',Bdata)








