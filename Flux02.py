import os, sys, numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


obslist = sys.argv[1]

if len(sys.argv) > 2:
    lowE = float(sys.argv[2])
    highE = float(sys.argv[3])

    le = int(round((lowE-1.6)/0.04))
    he = int(round((highE-1.6)/0.04))
else:
    lowE, highE = 3.0, 12.0
    le = int(round((3.0-1.6)/0.04))
    he = int(round((12.0-1.6)/0.04))

tag = obslist.split('/')[-1]
det = tag[4]
tag = tag.split('.')[0]
tag = tag.split('sun')[-1]
# get the data from all observations
# It will be assumed that the data will already be divided into the angles

OBS = np.genfromtxt(obslist, dtype='|U11')
# Need to go through all the PI values that are available like in the image build

# So, I need the mask, and I need the read in the PI, X, Y, DET1X, DET1Y, EXP, and GRADE

with fits.open('/uufs/astro.utah.edu/common/home/u1019304/temp/fullmaskA_final.fits') as hdu:
    maskA = hdu[0].data
with fits.open('/uufs/astro.utah.edu/common/home/u1019304/temp/fullmaskB_final.fits') as hdul:
    maskB = hdul[0].data

avg_flux_A = []
avg_flux_B = []

for ob in OBS:
        
    if det == 'A':
        mask = np.copy(maskA)
    else:
        mask = np.copy(maskB)

    event_file = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/{}/event_sep_cl/OCC/nu{}{}02_fullevts_SUN.fits'.format(ob,ob,det)

    if not os.path.isfile(event_file):
        continue

    with fits.open(event_file) as hdul:
        PI = hdul[1].data['PI']
        EXP = hdul[1].header['EXPOSURE']
        X = hdul[1].data['X']
        Y = hdul[1].data['Y']
        DET1_X = hdul[1].data['DET1X']
        DET1_Y = hdul[1].data['DET1Y']
        GRADE = hdul[1].data['GRADE']

    flux = 0
    idx_good = (X > 0)*(Y > 0)*(GRADE <= 26)

    for m in range(len(PI[idx_good])):
        if le <= PI[idx_good][m] < he:
            if (mask[DET1_Y[idx_good][m]-1,DET1_X[idx_good][m]-1] == 0):
                flux += 1
    flux /= EXP

    if det == 'A':
        avg_flux_A.append(flux)
    else:
        avg_flux_B.append(flux)


if det == 'A':
    checklist = [] 

    Astd = np.std(avg_flux_A)
    Amean = np.mean(avg_flux_A)
    Ax = np.arange(len(avg_flux_A))
    
    for i,check in enumerate(avg_flux_A):
        if check > 2*Astd:
            checklist.append(OBS[i])

    plt.plot(Ax, avg_flux_A, 'bo')
    plt.axhline(y=Amean+2*Astd,color='k',linestyle='--')
    plt.axhline(y=Amean+Astd,color='r',linestyle='--')
    plt.axhline(y=Amean-Astd,color='r',linestyle='--')
    plt.axhline(y=Amean, color = 'b', linestyle='-')
    plt.suptitle('A norm photon count with error bounds')
    plt.title('Energies from {} to {} at angles of {}'.format(lowE,highE,tag))
    plt.xlabel('OBSID Ref')
    plt.ylabel('Photon Count Value (adj for exposure)')
    plt.savefig('Acountavg_{}.png'.format(tag))

    plt.close()

if det == 'B':
    checklist = []

    Bstd = np.std(avg_flux_B)
    Bmean = np.mean(avg_flux_B)
    Bx = np.arange(len(avg_flux_B))

    for i,check in enumerate(avg_flux_B):
        if check > 2*Bstd:
            checklist.append(OBS[i])

    plt.plot(Bx, avg_flux_B, 'bo')
    plt.axhline(y=Bmean+2*Bstd,color='k',linestyle='--')
    plt.axhline(y=Bmean+Bstd,color='r',linestyle='--')
    plt.axhline(y=Bmean-Bstd,color='r',linestyle='--')
    plt.axhline(y=Bmean,color='b',linestyle='-')
    plt.suptitle('B norm photon count with error bounds')
    plt.title('Energies from {} to {} at angles of {}'.format(lowE,highE,tag))
    plt.xlabel('OBSID Ref')
    plt.ylabel('Photon Count Value (adj for exposure)')
    plt.savefig('Bcountavg_{}.png'.format(tag))


    plt.close()

if len(checklist) > 0:
    with open('highflux{}_{}.txt'.format(det,tag),'a+') as O:
        for ob in checklist:
            O.write('{}\n'.format(ob))

