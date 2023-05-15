#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntax: ExcessFlux.py obslist telescope[A/B] lowenergy highenergy

# Purpose of this code is to normalize the flux to the exposure time to get a rate
# and compare it to all other obs in the list. 


import sys, os, numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import time
from astropy import units as u
from astropy.coordinates import SkyCoord as sc

obslist = sys.argv[1]
det = sys.argv[2]
lowe = float(sys.argv[3])
highe = float(sys.argv[4])

######### Definitions ##############
def GetData(fi):
  with fits.open(fi) as hdul:
    data = hdul[1].data
  return data

def GetExpoMap(fi):
  with fits.open(fi) as hdul:
    expmap = hdul[0].data
  return expmap

def create_circ_mask(h,w,centerx,centery,radius):
  Y,X = np.ogrid[:h,:w]
  dist_from_center = np.sqrt((X-centerx)**2 + (Y-centery)**2)
  mask = dist_from_center <= radius
  return mask

def GetExcl(fi):
  f = open(fi,'r+')
  blah = []
  for line in f.readlines():
    blah.append(line)
  x = [];y=[];rad=[]
  for i in range(3,len(blah)):
    pos = blah[i][7:-2].rstrip().split(',')
    x.append(float(pos[0]))
    y.append(float(pos[1]))
    rad.append(float(pos[2]))
  return x,y,rad

def maxFlux(ar, limit):
  values = []; posn = []
  limit = int(limit*0.2)
  if limit > 50: limit = 50
  ranks = sorted( [(x,i) for (i,x) in enumerate(ar)], reverse=True )[:limit]
  for i in range(len(ranks)):
    values.append(ranks[i][0])
    posn.append(ranks[i][1])
  max_flux = np.max(values)
  avg_max = np.average(values)
  stddev = np.std(values)
  return max_flux, avg_max, stddev 

#######################################

try:
  obsid = np.genfromtxt(obslist,'|U11')
  obsdir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/'
  with fits.open('/uufs/astro.utah.edu/common/home/u1019304/temp/fullmask'+det+'_final.fits') as hdul:
    edgemask = hdul[0].data  

  oblist = []; fluxrate = []; avgflux=[]; stderr=[]
  writefilename = 'max_flux_rate_adj_'+det+'_'+str(lowe)+'_'+str(highe)+'_'+time.strftime("%Y%m%d-%H%M%S")+'.txt'

# Loop over obs
  for i,ob in enumerate(obsid):
    print(str(i)+' of '+str(len(obsid)))
    image = np.zeros((360,360))
    event = obsdir+ob+'/event_sep_cl/NOOCC/nu'+ob+det+'01_fullevts_NOSUN.fits'
    if not os.path.isfile(event): continue
    regcheck = obsdir+ob+'/event_defcl/excl.reg'
    if os.path.isfile(regcheck):
      srcfile = open(regcheck,'r+')
      x_vec,y_vec,rad = GetExcl(regcheck)
      grid_pi = np.ones((1000,1000))
      for q,vec in enumerate(x_vec):
        mask = create_circ_mask(1000,1000,float(vec)-1,float(y_vec[q])-1,float(rad[q]))
        grid_pi[mask] = 0
      nustar_cxb = obsdir.replace('OBS/','CXB/01/expmap/nu'+ob+det+'_nosun')
      exp_map = GetExpoMap(nustar_cxb)
    else:
      grid_pi = np.ones((1000,1000))
      exp_map = np.ones((360,360))
    
    em = np.copy(edgemask)
    em = em.astype('float64')
    with fits.open(event) as hdul:
      PI = hdul[1].data['PI']
      EXP = hdul[1].header['EXPOSURE']
      RA = hdul[0].header['RA_PNT']
      DEC = hdul[0].header['DEC_PNT']
      X = hdul[1].data['X']
      Y = hdul[1].data['Y']
      DET1_X = hdul[1].data['DET1X']
      DET1_Y = hdul[1].data['DET1Y']
    c = sc(ra = RA*u.degree, dec = DEC*u.degree, frame='fk5')
    if abs(c.galactic.b.deg) < 10.0: continue
    elow = int(round((lowe-1.6)/0.04))
    ehigh = int(round((highe-1.6)/0.04))
    idx_good = (X>0)*(Y>0)*(PI>=elow)*(PI<ehigh)
    ii = 0
    for m in range(len(PI[idx_good])):
      if (em[DET1_Y[idx_good][m]-1,DET1_X[idx_good][m]-1] == 0) or (grid_pi[Y[idx_good][m]-1,X[idx_good][m]-1] == 0):
        continue
      else:
        image[DET1_Y[idx_good][m]-1,DET1_X[idx_good][m]-1] += 1
        ii += 1
    if ii < 20:
      print('not enought counts: '+str(ii)) 
      continue
    eM = np.copy(exp_map)
    eM /= EXP
    im = np.copy(image)
    im = (im/EXP)*eM
    im = np.ravel(im)
    max_flux_rate, avg_flux_rate, stddev = maxFlux(im, ii) 
    #if max_flux_rate > (avg_flux_rate+stddev):
    #  max_flux_rate = avg_flux_rate+stddev
    #print(str(max_flux_rate)) 
    with open('/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs/'+writefilename,'a+') as O:
       O.write(ob+' '+str(max_flux_rate)+' '+str(avg_flux_rate)+' '+str(stddev)+' '+str(elow)+' '+str(ehigh)+'\n')
    oblist.append(ob) 
    fluxrate.append(max_flux_rate)
    avgflux.append(avg_flux_rate)
    stderr.append(stddev)
  # This is where we will plot the results:
  plt.scatter(oblist,fluxrate,s=1,c="blue",alpha=0.3)
  plt.scatter(oblist,avgflux,s=1,c="red",alpha=0.3)
  plt.errorbar(oblist,avgflux,yerr=stderr,linestyle='None',color="red")
  plt.show()

except Exception as e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)


