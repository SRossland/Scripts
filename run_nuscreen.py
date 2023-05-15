import os, sys, numpy as np

#obsids = np.genfromtxt(sys.argv[1])
#obsids = [int(i) for i in obsids]

with open(sys.argv[1]) as f:
    lines = f.readlines()

obsids = []
for i in range(2,len(lines)):
    li = lines[i].split()
    obsids.append(li[0])

obsDirectory = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS'
checkDirectory = obsDirectory.replace('OBS','OneElevations')

for i,ob in enumerate(obsids):
    # Check if the file exists
    for det in ['A','B']:
        if os.path.isfile('{}/nu{}{}02_1ELV.fits'.format(checkDirectory,ob,det)):
            continue
        os.system('nuscreen_pipe.sh {}/{} {}'.format(obsDirectory,ob,det))







