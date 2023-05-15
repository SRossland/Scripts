#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# This is just a quick program to check the images of obs from the data 
# to see if there is anything going on with the data, like a shadow, source
# etc.

# syntax: check_image.py path/to/list det low_energy high_energy idx=None

import os, sys, numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

obs_list = sys.argv[1]
de = sys.argv[2]
elow = sys.argv[3]
ehigh = sys.argv[4]
if len(sys.argv) == 6:
  idx = sys.argv[5]
else: idx = None
################################################
def getevt(fi):
    with fits.open(fi) as hdul:
      PI = hdul[1].data['PI']
      X = hdul[1].data['DET1X']
      Y = hdul[1].data['DET1Y']
    return PI, X, Y

def getmask(det):
  mask_dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/fullmask'+det+'.fits'
  with fits.open(mask_dir) as hdul:
    em = hdul[0].data
  return em

def getexcl(excl_dir):
  srcfile = open(excl_dir,'r+')
  x_vec,y_vec,rad,data = [],[],[],[]
  for line in srcfile.readlines():
    data.append(line.rstrip())
  srcfile.close()
  for q in range(len(data[3:])):
    xv, yv, radv = data[3+q].split(',')
    x_vec.append(xv[7:])
    y_vec.append(yv)
    rad.append(radv[:-1])
  [float(t) for t in x_vec]
  [float(t) for t in y_vec]
  [float(t) for t in rad]
  return x_vec, y_vec, rad
###############################################


obsid = []
saved_obs = []

f = open(obs_list,'r+')
for line in f.readlines():
  obsid.append(line.rstrip())
f.close()

if idx == None:
  idx_start = 0
else:
  idx_start = obsid.index(idx)

img = np.zeros((360,360))

for i in range(idx_start,len(obsid)):
  ex_path = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/'+obsid[i]+'/event_defcl/excl.reg'
  fi = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/'+obsid[i]+'/event_cl/nu'+obsid[i]+de+'01_cl.evt'
#  if os.path.isfile(ex_path) == False: continue
  if os.path.isfile(fi):
    PI,X,Y = getevt(fi)
  else: continue
  mask = getmask(de)
  for m in range(len(PI)):
    if float(elow) <= PI[m]*0.04+1.6 <= float(ehigh):
      if (X[m] <= -1) or (mask[Y[m]-1,X[m]-1] == 0):
        continue
      else:
        img[Y[m]-1,X[m]-1] += 1
  if np.sum(img) != 0.0: 
    fig, ax = plt.subplots(1)
    ax.imshow(img)
    plt.title(str(obsid[i]))
    plt.show()
#  if os.path.isfile(ex_path):
#      x_vec, y_vec, rad = getexcl(ex_path)
#      if len(x_vec) >= 2:
#        os.chdir('/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/'+obsid[i]+'/event_defcl')
#        os.system('ds9 /uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/'+obsid[i]+'/event_cl/nu'+obsid[i]+de+'01_cl.evt')
    #plt.show()
    #sa = raw_input("Save obsid?  [y,n]:  ")
    #if sa == 'y':
      #saved_obs.append(obsid[i])
  img = None; img = np.zeros((360,360))
  #print(str(i)+" done")

print('End of list')
'''
if len(saved_obs) > 0:
  print(str(len(saved_obs)))
  for i in range(len(saved_obs)):
    with open('OBSIDS_TO_CHECK.txt','a+') as O:
      O.write(saved_obs[i]+'\n')
'''

