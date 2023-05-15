#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

import os, sys, numpy as np


#obsids = ['10110001002','10110003002','10110004002','10110005001']
obsids = ['10110005001']

# need to pass on is the telescope letter, data, low energy, high energy, obsid, and aCXB norm value

homedir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/CRAB/'

el, eh, ACXB, BCXB = [],[],[],[]

with open(homedir+'A_params.txt','r+') as f:
  lines = f.readlines()
for line in lines:
  el.append(line.split()[0])
  eh.append(line.split()[1])
  ACXB.append(line.split()[2])

expA = float(lines[0].split()[5])
with open(homedir+'B_params.txt','r+') as f:
  lines = f.readlines()
for line in lines:
  BCXB.append(line.split()[2])
expB = float(lines[0].split()[5])


#for i,ob in enumerate(obsids):
for k,ob in enumerate(obsids):
#  i = 1
#  i = 2
  i = 3
  if i == 0: dets = ['A','B']
  if (i == 1) or (i == 3):  dets = 'B'
  if i == 2: dets = 'A'
  
  for det in dets:
    for j in range(len(el)):
      data = homedir+ob+'/'+det+'_'+el[j]+'_'+eh[j]+'keV_crab.fits'
      if det == 'A': aCXBnorm = float(ACXB[j])/expA
      if det == 'B': aCXBnorm = float(BCXB[j])/expB
      os.system('count_stat_crab.py '+det+' '+data+' '+el[j]+' '+eh[j]+' '+ob+' '+str(aCXBnorm))
  



