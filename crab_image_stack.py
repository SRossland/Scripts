#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# syntax: ./crab_image.py 

import numpy as np, os, sys, datetime, math
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord as sc
import matplotlib.pyplot as plt


mask_dir = '/uufs/astro.utah.edu/common/home/u1019304/temp'
with fits.open(mask_dir+'/fullmaskA_final.fits') as mAhdu:
    edgemaskA = mAhdu[0].data

with fits.open(mask_dir+'/fullmaskB_final.fits') as mBhdu:
    edgemaskB = mBhdu[0].data

Afull = np.zeros((360,360))
Bfull = np.zeros((360,360))

def create_circ_mask(h,w,centerx,centery,radius):
        # Here is am again using X and Y, but they are local and should be treated as such by python
        Y,X = np.ogrid[:h,:w]
        dist_from_center = np.sqrt((X-centerx)**2 + (Y-centery)**2)
        mask = dist_from_center <= radius
        return mask

def getMask(obsid,detector):
  with fits.open(home_dir+obsid+'/event_cl/mask'+detector+'.fits') as hdul:
    mask = hdul[0].data
  return mask

home_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/CRAB/'

obslist = ['10110001002','10110003002','10110004002','10110005001']



el, eh, ACXB, BCXB = [],[],[],[]
with open('A_params.txt','r+') as f:
  lines = f.readlines()
for line in lines:
  el.append(line.split()[0])
  eh.append(line.split()[1])

for j in range(len(el)):
  for i,ob in enumerate(obslist):
    if i == 0:  
      det = ['A','B']
      for de in det:
        fi = home_dir+ob+'/'+de+'_'+el[j]+'_'+eh[j]+'keV_crab.fits'
        fits.writeto(fi,np.zeros((360,360)))
    if (i == 1) or (i == 3):
      fi = home_dir+ob+'/B_'+el[j]+'_'+eh[j]+'keV_crab.fits'
      fits.writeto(fi,np.zeros((360,360)))
    if i == 2:
      fi = home_dir+ob+'/A_'+el[j]+'_'+eh[j]+'keV_crab.fits'
      fits.writeto(fi,np.zeros((360,360)))

try:
  for i,ob in enumerate(obslist):
    if i == 0:
      maskA = getMask(ob,'A')
      maskB = getMask(ob,'B')
    if (i == 1) or (i == 3):
      maskB = getMask(ob,'B')
    if i == 2:
      maskA = getMask(ob,'A')
    
    for j in range(len(el)):
      if i == 0: dets = ['A','B']
      if (i == 1) or (i == 3): dets = 'B'
      if i == 2:  dets = 'A'
      for det in dets:
        fi = home_dir+ob+'/'+det+'_'+el[j]+'_'+eh[j]+'keV_crab.fits'
        event = home_dir+ob+'/event_cl/nu'+ob+det+'01_cl.evt'
        with fits.open(event) as hdul:
          PI = hdul[1].data['PI']
          EXP = hdul[1].header['EXPOSURE']
          X = hdul[1].data['X']
          Y = hdul[1].data['Y']
          DET1_X = hdul[1].data['DET1X']
          DET1_Y = hdul[1].data['DET1Y']


        idx_good = (X > 0)*(Y > 0)

        if det == 'A': mask = np.copy(maskA); em = np.copy(edgemaskA)
        if det == 'B': mask = np.copy(maskB); em = np.copy(edgemaskB)

        for m, pim in enumerate(PI[idx_good]): 
          elow_int = int(round((float(el[j])-1.6)/0.04))
          ehigh_int = int(round((float(eh[j])-1.6)/0.04))
          if elow_int <= pim < ehigh_int:
            if (mask[DET1_Y[idx_good][m]-1,DET1_X[idx_good][m]-1] == 0) or (em[DET1_Y[idx_good][m]-1,DET1_X[idx_good][m]-1] == 0):
              continue
            else:
              if det == 'A':
                Afull[DET1_Y[idx_good][m]-1,DET1_X[idx_good][m]-1] += 1
              else:
                Bfull[DET1_Y[idx_good][m]-1,DET1_X[idx_good][m]-1] += 1
       
         
        with fits.open(fi, mode='update') as hdul:
          if det == 'A': hdul[0].data = Afull
          if det == 'B': hdul[0].data = Bfull
          hdul.flush() 
        
        Afull = np.zeros((360,360))
        Bfull = np.zeros((360,360))
        mask = None; elow_int = None; ehigh_int = None; em = None

except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(str(exc_type) + '\n' + str(e) +'\n' +str(exc_tb.tb_lineno)+ '\n')

