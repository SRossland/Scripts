#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# syntax: run_det_errors.py telescope[A/B] working_directory

# This program is a wrapper to run Det_errors.py

import os, sys, numpy as np

detp = sys.argv[1]
root_file = sys.argv[2]



le, he = [],[]
dir_list = os.listdir(root_file)
for li in dir_list:
  if 'longformat' in li:
    if detp+'_nosun' in li:
      le.append(float(li.split('_')[2]))
      he.append(float(li.split('_')[3]))

le,he = list(zip(*sorted(zip(le,he))))

le = [str(i) for i in le]
he = [str(i) for i in he]

for det in detp:
  detmap = '/uufs/astro.utah.edu/common/home/u1019304/temp/detmap'+det+'.fits'
  aCXBgrad = '/uufs/astro.utah.edu/common/home/u1019304/my_idl/nuskybgd/auxil/det'+det+'_det1.img'
  for i,low in enumerate(le):
    data_file = root_file+'/fits_file/Data'+det+'_01_nosun_'+low+'_'+he[i]+'keV.fits'
    bins_file = data_file.replace('Data','BINS')
    exp_file = data_file.replace('Data','Exp')
    norm_file = root_file+'/'+det+'_nosun_'+low+'_'+he[i]+'_keVparams_longformat.txt'
    os.system('Det_errors.py '+data_file+' '+bins_file+' '+detmap+' '+exp_file+' '+aCXBgrad+' '+norm_file)







