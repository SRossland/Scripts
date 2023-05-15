#!/uufs/astro.utah.edu/common/home/u1019304/VENV3.7.9/bin/python3

# Syntax: BoreSightCheck.py obslist lower_expected higher_expected

# Purpose:  To verify the rotation results are in line with expectation by 
#           comparing those values with the boresight to sun angle in the 
#           NuSTAR attorb files.


import sys, os, numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

obslist = sys.argv[1]
#lower_expected_angle = float(sys.argv[2])
#higher_expected_angle = float(sys.argv[3])


obsdir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS'
obs = np.genfromtxt(obslist,dtype='|U11')

boresight_angle, PA_Roll_angle = [],[]

for i,ob in enumerate(obs):
    eventA = '{}/{}/event_cl/nu{}A.attorb'.format(obsdir,ob,ob)
    eventB = eventA.replace('A.attorb','B.attorb')
    if os.path.isfile(eventA): 
        event = eventA
    else:
        event = eventB
    with fits.open(event) as hdul:
        BSA = hdul[1].data['SUN_ANGLE']
        Roll = hdul[1].data['ROLL']
    Roll_avg = np.average(Roll)
    BSA_avg = np.average(BSA)
    boresight_angle.append(BSA_avg)
    PA_Roll_angle.append(Roll_avg)
#    if lower_expected_angle <= BSA_avg <= higher_expected_angle:
#        print('{} is good'.format(ob))
#    else:
#        print('{} is bad: {} '.format(ob, BSA_avg))

#fi_name_end = obslist.split('_')[3]
fi_name = '02_OBS_BSA_Master.txt'#.format(fi_name_end)

for i,ob in enumerate(obs):
    with open(fi_name,'a+') as O:
        O.write('{} {} {} {}'.format(ob,boresight_angle[i],PA_Roll_angle[i],'\n'))



