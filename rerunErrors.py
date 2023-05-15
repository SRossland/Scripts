#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntax: rerunErrors.py 

# Written by: Steven Rossland 03/21

# Purpose: To rerun count_stat_new.py and runDetErrors.py to ensure proper error finding

import os, sys, numpy as np
from astropy.io import fits

# To run count_stat_new.py, we'll need to use the telescope letter [A,B], Data image, 
#      low energy, high energy, seperation type [full, sun, nosun], and working directory

# to get those, i'll assume this program is running from the directory where those files are stored, i.e. /uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/CXB/[01,02]/[full,sun,nosun],****steps/

def runDiffEnergy():
  for i,low in enumerate(elowA_sort):
    high = ehighA_sort
    os.system('count_stat_new.py A fits_file/DataA_01_nosun_'+str(low)+'_'+str(high)+'keV.fits '+str(low)+' '+str(high)+' nosun .')
  for i,low in enumerate(elowB_sort):
    high = ehighB_sort
    os.system('count_stat_new.py B fits_file/DataB_01_nosun_'+str(low)+'_'+str(high)+'keV.fits '+str(low)+' '+str(high)+' nosun .')
  os.system('runDetErrors.py')

listroot = os.listdir('./old_params_used_1st_newmask')

elowA,ehighA,elowB,ehighB = [],[],[],[]
for i,li in enumerate(listroot):
  if 'A_nosun' in li:
    if 'longformat' in li:
      crap = li.split('_')
      elowA.append(crap[2])
      ehighA.append(crap[3])
  if 'B_nosun' in li:
    if 'longformat' in li:
      crap = li.split('_')
      elowB.append(crap[2])
      ehighB.append(crap[3])

if len(elowA) != len(elowB):
  difference = True
else: 
  difference = False

elowA = [float(i) for i in elowA]
elowB = [float(i) for i in elowB]
ehighA = [float(i) for i in ehighA]
ehighB = [float(i) for i in ehighB]

energyA = sorted(zip(elowA,ehighA), key = lambda x: x[0])
energyB = sorted(zip(elowB,ehighB), key = lambda x: x[0])

elowA_sort, ehighA_sort = zip(*energyA)
elowB_sort, ehighB_sort = zip(*energyB)

homedir = os.getcwd()

if difference: runDiffEnergy()
else:
  for i,low in enumerate(elowA_sort):
    high = ehighA_sort[i]
    for tele in ['A','B']:
      os.system('count_stat_new.py '+tele+' fits_file/Data'+tele+'_01_nosun_'+str(low)+'_'+str(high)+'keV.fits '+str(low)+' '+str(high)+' nosun '+homedir)
  os.system('runDetErrors.py')

