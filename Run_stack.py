#!/uufs/astro.utah.edu/common/home/u1019304/VENV3.7.9/bin/python3

import os, sys
import numpy as np

strip_list, sl_rr, sl_lr = [],[],[]

dir_list = os.listdir('/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/SolarData')

for li in dir_list:
    if '02_OBS_sun_BSA_' in li:
        if 'RR.txt' in li:
            sl_rr.append(li)
            continue
        if 'RL.txt' in li:
            sl_lr.append(li)
            continue
        else:
            strip_list.append(li)

obsdir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS'

def StackList(lis,directory,lang,hang):
    os.system('Sun_image_stack.py '+obsdir+' '+lis+' 3 5 BOTH '+directory)
    for det in ['A','B']:
        os.system('mv '+directory+'/Data'+det+'_02_nosun_3.0_5.0keV.fits '+directory+'/Data'+det+'_02_nosun_3_5_'+lang+'_'+hang+'.fits')
        os.system('mv '+directory+'/Norm'+det+'_02_nosun_3.0_5.0keV.fits '+directory+'/Norm'+det+'_02_nosun_3_5_'+lang+'_'+hang+'.fits')
        os.system('rm '+directory+'/Exp*.fits')


home_dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/SolarData'
strip_dir = os.path.join(home_dir,'BSA_strip_fits/nosun_fits')
split_dir = os.path.join(home_dir,'BSA_split_strip_fits/nosun_fits')

for li in strip_list:
    stuff = li.split('_')
    low_angle = stuff[4].split('to')[0]
    high_angle = stuff[4].split('to')[1].rstrip('.txt')
    StackList(li,strip_dir,low_angle,high_angle)

for li in sl_rr:
    stuff = li.split('_')
    low_angle = stuff[4].split('to')[0]
    high_angle = stuff[4].split('to')[1].rstrip('.txt')
    StackList(li,split_dir,low_angle,high_angle)

for li in sl_lr:
    stuff = li.split('_')
    low_angle = stuff[4].split('to')[0]
    high_angle = stuff[4].split('to')[1].rstrip('.txt')
    StackList(li,split_dir,low_angle,high_angle)

print('done')

