#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntax: 20minus38.py obsid le1 he1 le2 he2 PLOTorSAVE

# Purpose: The purpose of this script is take energy bands of images from a single obs and subtract them from the another to give a residual that may show an underlying issue or signature in the image itself that is energy dependent

# Inputs: obs = Obsid, this is assumed to be an obs that is found in the OBS directory of the wik file structure
#	  le1 = the lower energy limit of the first energy band image
#  	  he1 = the high energy limit of the first energy band image
#	  le2 = the lower energy limit of the 2nd image
#  	  he2 = the high energy limit of the 2nd image

import os, sys, numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

obsid = sys.argv[1]
le1 = sys.argv[2]
he1 = sys.argv[3]
le2 = sys.argv[4]
he2 = sys.argv[5]
SorP = sys.argv[6]

elow1 = int(round((float(le1)-1.6)/0.04))
ehigh1 = int(round((float(he1)-1.6)/0.04))
elow2 = int(round((float(le2)-1.6)/0.04))
ehigh2 = int(round((float(he2)-1.6)/0.04))

save_dir = '/uufs/astro.utah.edu/common/home/u1019304/temp/flux_images/'

# To see if the obsid is a file:
obs = []
if os.path.isfile(obsid):
  with open(obsid,'r+') as f:
    lines = f.readlines()
  for line in lines:
    obs.append(line.split()[0])

###############################

def GetExcl(fi):
  f = open(fi,'r+')
  blah = []
  for line in f.readlines():
    blah.append(line)
  x,y,rad = [],[],[]
  for i in range(3,len(blah)):
    pos = blah[i][7:-2].rstrip().split(',')
    x.append(float(pos[0]))
    y.append(float(pos[1]))
    rad.append(float(pos[2]))
  return x,y,rad

def create_circ_mask(h,w,centerx,centery,radius):
  Y,X = np.ogrid[:h,:w]
  dist_from_center = np.sqrt((X-centerx)**2 + (Y-centery)**2)
  mask = dist_from_center <= radius
  return mask

def GetExpoMap(fi):
  with fits.open(fi) as hdul:
    expmap = hdul[0].data
  return expmap

###############################

obsdir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/'
for i,ob in enumerate(obs):
  imageA1 = np.zeros((360,360))
  imageB1 = np.zeros((360,360))
  imageA2 = np.zeros((360,360))
  imageB2 = np.zeros((360,360))
  for det in ['A','B']:
    event = obsdir+ob+'/event_sep_cl/NOOCC/nu'+ob+det+'01_fullevts_NOSUN.fits'
    if not os.path.isfile(event): continue
    regcheck = obsdir+ob+'/event_defcl/excl.reg'
    if os.path.isfile(regcheck):
      srcfile = open(regcheck,'r+')
      x_vec,y_vec,rad = GetExcl(regcheck)
      grid_pi = np.ones((1000,1000))
      for q,vec in enumerate(x_vec):
        mask = create_circ_mask(1000,1000,float(vec)-1,float(y_vec[q])-1,float(rad[q]))
        grid_pi[mask]=0
      nustar_cxb = obsdir.replace('OBS/','CXB/01/expmap/nu'+ob+det+'_nosun')
      exp_map = GetExpoMap(nustar_cxb)
    else:
      grid_pi = np.ones((1000,1000))
      exp_map = np.ones((360,360))
    
    with fits.open('/uufs/astro.utah.edu/common/home/u1019304/temp/fullmask'+det+'_final.fits') as hdul:
      edgemask = hdul[0].data
    em = np.copy(edgemask)
    em = em.astype('float64')
    with fits.open(event) as hdul:
      PI = hdul[1].data['PI']
      X = hdul[1].data['X']
      Y = hdul[1].data['Y']
      DET1_X = hdul[1].data['DET1X']
      DET1_Y = hdul[1].data['DET1Y']
    idx_good1 = (X>0)*(Y>0)*(PI>=elow1)*(PI<ehigh1)
    idx_good2 = (X>0)*(Y>0)*(PI>=elow2)*(PI<ehigh2)

    for m in range(len(PI[idx_good1])):
      if (em[DET1_Y[idx_good1][m]-1,DET1_X[idx_good1][m]-1] == 0) or (grid_pi[Y[idx_good1][m]-1,X[idx_good1][m]-1] == 0):
        continue
      else:
        if det == 'A':
          imageA1[DET1_Y[idx_good1][m]-1,DET1_X[idx_good1][m]-1] += 1
        else:
          imageB1[DET1_Y[idx_good1][m]-1,DET1_X[idx_good1][m]-1] += 1
  
    for m in range(len(PI[idx_good2])):
      if (em[DET1_Y[idx_good2][m]-1,DET1_X[idx_good2][m]-1] == 0) or (grid_pi[Y[idx_good2][m]-1,X[idx_good2][m]-1] == 0):
        continue
      else:
        if det == 'A':
          imageA2[DET1_Y[idx_good2][m]-1,DET1_X[idx_good2][m]-1] += 1
        else:
          imageB2[DET1_Y[idx_good2][m]-1,DET1_X[idx_good2][m]-1] += 1

# Here is the image part of the program, it is assumed that the le1-he1 will have higher counts due to the nature of Xrays
  
  
  finalA = imageA1-imageA2
  finalB = imageB1-imageB2
  if SorP.lower() == 'plot':
    f,(ax1,ax2) = plt.subplots(1,2)
    f.suptitle(ob+' '+le1+' to '+he1+' sub '+le2+' to '+he2)
    ax1.matshow(finalA)
    ax2.matshow(finalB)
    plt.show()

  if SorP.lower() == 'save':
    fits.writeto(save_dir+ob+'_A_'+le1+'to'+he1+'sub'+le2+'to'+he2+'.fits',finalA)
    fits.writeto(save_dir+ob+'_B_'+le1+'to'+he1+'sub'+le2+'to'+he2+'.fits',finalB)
    fits.writeto(save_dir+ob+'_A_'+le1+'to'+he1+'.fits',imageA1)
    fits.writeto(save_dir+ob+'_A_'+le2+'to'+he2+'.fits',imageA2)
    fits.writeto(save_dir+ob+'_B_'+le1+'to'+he1+'.fits',imageB1)
    fits.writeto(save_dir+ob+'_B_'+le2+'to'+he2+'.fits',imageB2)




