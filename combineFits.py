#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntax: *.py directory A/B low_energy high_energy output_file 

import os, sys, numpy as np
from astropy.io import fits

########################

def GetData(fil):
  with fits.open(fil) as hdul:
    Dat = hdul[0].data
  return Dat

########################


workDir = sys.argv[1]
det = sys.argv[2]
typ = sys.argv[3]
le = float(sys.argv[4])
he = float(sys.argv[5])
save_name = sys.argv[6]

DirList = os.listdir(workDir)
workList = []

for fi in DirList:
  if typ+det in fi:
    if (float(fi.split('_')[3]) >= le) and (float(fi.split('_')[4].split('keV')[0]) <= he):
      workList.append(fi)

data_grid = np.zeros((360,360))
for fi in workList:
  f = os.path.join(workDir,fi)
  dat = GetData(f)
  data_grid += dat

fits.writeto(save_name,data_grid)


