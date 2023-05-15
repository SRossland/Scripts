#!/bin/python3

import os, sys, numpy as np

working_dir = sys.argv[1]

# Det_errors.py Data Bins Exp gradfile normilization_value_file

list_root = os.listdir(working_dir)

Anorms, Bnorms = [],[]
elowA, ehighA = [],[]
elowB, ehighB = [],[]
for i, li in enumerate(list_root):
  if 'A_nosun' in li:
    if 'longformat.txt' in li:
      Anorms.append(li)
      crap = li.split('_')
      elowA.append(crap[2])
      ehighA.append(crap[3])
  if 'B_nosun' in li:
    if 'longformat.txt' in li:
      Bnorms.append(li)
      crap = li.split('_')
      elowB.append(crap[2])
      ehighB.append(crap[3])

# May want to sort the energy lists so it goes in order, but it's not needed

elowA = [float(i) for i in elowA]
#ehighA = [float(i) for i in ehighA]
elowB = [float(i) for i in elowB]
#ehighB = [float(i) for i in ehighB]

energyA = sorted(zip(elowA,ehighA,Anorms), key = lambda x: x[0])
energyB = sorted(zip(elowB,ehighB,Bnorms), key = lambda x: x[0])

elowA, ehighA, Anorms = zip(*energyA)
elowB, ehighB, Bnorms = zip(*energyB)

# comes out as a tuple, need to make them a list
elowA = list(elowA)
ehighA = list(ehighA)
elowB = list(elowB)
ehighB = list(ehighB)
Anorms = list(Anorms)
Bnorms = list(Bnorms)

# make them str elements *this could all be compressed to a single line per variable.....but no
elowA = [str(i) for i in elowA]
#ehighA = [str(i) for i in ehighA]
elowB = [str(i) for i in elowB]
#ehighB = [str(i) for i in ehighB] 

# Checks:
if (len(elowA) != len(ehighA)) or (len(elowB) != len(ehighB)): print('different energy list lengths'); sys.exit()
if (len(elowA) != len(Anorms)) or (len(elowB) != len(Bnorms)): print('different file list lengths to energy list lenghts'); sys.exit()
if len(elowA) != len(Bnorms): print('different B file list lengths: '+str(len(elowA))+' '+str(len(Bnorms)))


tele = ['A','B']
norms = [Anorms, Bnorms]
for det in [0,1]:
  detfile = '/uufs/astro.utah.edu/common/home/u1019304/temp/detmap'+tele[det]+'.fits'
  gradfile = '/uufs/astro.utah.edu/common/home/u1019304/my_idl/nuskybgd/auxil/det'+tele[det]+'_det1.img'
  for i,norm in enumerate(norms[det]):
    if tele[det] == 'A':
      el = elowA[i]
      eh = ehighA[i]
    if tele[det] == 'B':
      el = elowB[i]
      eh = ehighB[i]
    datafile = 'fits_file/Data'+tele[det]+'_01_nosun_'+el+'_'+eh+'keV.fits'
    binfile = 'fits_file/BINS'+tele[det]+'_01_nosun_'+el+'_'+eh+'keV.fits'
    expfile = 'fits_file/Exp'+tele[det]+'_01_nosun_'+el+'_'+eh+'keV.fits'
    os.system('Det_errors.py '+datafile+' '+binfile+' '+detfile+' '+expfile+' '+gradfile+' '+norm)

print('full list finished')
 
# think of passing the det and energy values since i'm gathering them here 
