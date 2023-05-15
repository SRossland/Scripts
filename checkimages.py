#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

import os, sys, numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

obs_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/CXB/01/nosun/'
blankdir = '3.0_41.28_36steps_BlankSky_newmask'
excl18stepdir = '3.0_10.52_18steps_excl_only'
bin10countdir = '3.0_10.08_20steps_10countbin_exclonly'
bin50countdir = '3.0_10.08_20steps_10countbin_exclonly'
savedir = obs_dir+'countrateimages/'
savefit = obs_dir+'countratefits/'
#
#####defs ##########

def getEnergies(directory):
  li = os.listdir(directory)
  low, high = [],[]
  for i,blah in enumerate(li):
    if 'longformat' in blah:
      low.append(blah.split('_')[2])
      high.append(blah.split('_')[3])
  return low, high

def orderList(direct):
  lowener, highener = getEnergies(direct)
  lo = list(set(lowener))
  hi = list(set(highener))
  lo = [float(i) for i in lo]
  hi = [float(i) for i in hi]
  lo = list(np.sort(lo))
  hi = list(np.sort(hi))
  lo = [str(i) for i in lo]
  hi = [str(i) for i in hi]
  return lo, hi

def getFits(fi):
  with fits.open(fi) as hdul:
    data = hdul[0].data
  return data

def saveFits(fi, data):
  fits.writeto(fi, data)

#########################

blanklo, blankhi = orderList(obs_dir+blankdir)
lo18step, hi18step = orderList(obs_dir+excl18stepdir)
lo10count, hi10count = orderList(obs_dir+bin10countdir)
lo50count, hi50count = orderList(obs_dir+bin50countdir)

##### Same bin width compared #######
tot18A, tot18B, totblankA, totblankB = np.zeros((360,360)), np.zeros((360,360)), np.zeros((360,360)), np.zeros((360,360))
totnormA, totnormB = np.zeros((360,360)), np.zeros((360,360))
for i in range(len(lo18step)):

  for det in ['A','B']:
    blankdata = getFits(obs_dir+blankdir+'/fits_file/Data'+det+'_01_nosun_'+blanklo[i]+'_'+blankhi[i]+'keV.fits')
    blankexp = getFits(obs_dir+blankdir+'/fits_file/Exp'+det+'_01_nosun_'+blanklo[i]+'_'+blankhi[i]+'keV.fits')
    excl18data = getFits(obs_dir+excl18stepdir+'/fits_file/Data'+det+'_01_nosun_'+blanklo[i]+'_'+blankhi[i]+'keV.fits')
    excl18exp = getFits(obs_dir+excl18stepdir+'/fits_file/Exp'+det+'_01_nosun_'+blanklo[i]+'_'+blankhi[i]+'keV.fits')
    blnorm = np.divide(blankdata, blankexp, out=np.zeros_like(blankdata), where=blankexp!=0)
    ex18norm = np.divide(excl18data, excl18exp, out=np.zeros_like(excl18data), where=excl18exp!=0)
    saveFits(savefit+'Blanksky'+det+blanklo[i]+'_'+blankhi[i]+'.fits', blnorm)
    saveFits(savefit+'18step'+det+blanklo[i]+'_'+blankhi[i]+'.fits', ex18norm)
#    plt.matshow(blnorm, origin='lower')
#    plt.colorbar()
#    plt.savefig(savedir+'Blanksky'+det+blanklo[i]+'_'+blankhi[i]+'.png')
#    plt.matshow(ex18norm, origin='lower')
#    plt.savefig(savedir+'18step'+det+blanklo[i]+'_'+blankhi[i]+'.png')
## Do the bin width calc for blank HERE!!!
    binwidth = float(blankhi[i])-float(blanklo[i])
    blnorm_binwidnorm = blnorm/binwidth
#    plt.matshow(blnorm_binwidnorm, origin='lower')
#    plt.savefig(savedir+'Blanksky_binadjusted'+det+blanklo[i]+'_'+blankhi[i]+'.png')
#    plt.close('all')
    if det == 'A':
      tot18A += ex18norm
      totblankA += blnorm
      totnormA += blnorm_binwidnorm
    else:
      tot18B += ex18norm
      totblankB += blnorm
      totnormB += blnorm_binwidnorm
####### Now we have to fix by bin width ##########

saveFits(savefit+'18stepAtot.fits', tot18A)
saveFits(savefit+'18stepBtot.fits', tot18B)
saveFits(savefit+'BlankskytotA.fits', totblankA)
saveFits(savefit+'BlankskytotB.fits', totblankB)
saveFits(savefit+'Blanksky_binadjtotA.fits', totnormA)
saveFits(savefit+'Blanksky_binadjtotB.fits', totnormB)
subA = totblankA-tot18A; subB = totblankB-tot18B
saveFits(savefit+'Blank18subA.fits', subA)
saveFits(savefit+'Blank18subB.fits', subB)

#plt.matshow(tot18A, origin='lower')
#plt.savefig(savedir+'18stepAtot.png')
#plt.matshow(tot18B, origin='lower')
#plt.savefig(savedir+'18stepBtot.png')
#plt.matshow(totblankA, origin='lower')
#plt.savefig(savedir+'BlankskytotA.png')
#plt.matshow(totblankB, origin='lower')
#plt.savefig(savedir+'BlankskytotB.png')
#plt.matshow(totnormA, origin='lower')
#plt.savefig(savedir+'Blanksky_binadjtotA.png')
#plt.matshow(totnormB, origin='lower')
#plt.savefig(savedir+'Blanksky_binadjtotB.png')

#plt.matshow(totblankA-tot18A, origin='lower')
#plt.savefig(savedir+'Blank18subA.png')
#plt.matshow(totblankB-tot18B, origin='lower')
#plt.savefig(savedir+'Blank18subB.png')

#plt.close('all')

tot10A, tot10B, tot50A, tot50B = np.zeros((360,360)), np.zeros((360,360)), np.zeros((360,360)), np.zeros((360,360))

for i in range(len(lo10count)):
  for det in ['A','B']:
    bin10data = getFits(obs_dir+bin10countdir+'/fits_file/Data'+det+'_01_nosun_'+lo10count[i]+'_'+hi10count[i]+'keV.fits')
    bin50data = getFits(obs_dir+bin50countdir+'/fits_file/Data'+det+'_01_nosun_'+lo50count[i]+'_'+hi50count[i]+'keV.fits')
    bin10exp = getFits(obs_dir+bin10countdir+'/fits_file/Exp'+det+'_01_nosun_'+lo10count[i]+'_'+hi10count[i]+'keV.fits')
    bin50exp = getFits(obs_dir+bin50countdir+'/fits_file/Exp'+det+'_01_nosun_'+lo50count[i]+'_'+hi50count[i]+'keV.fits')
    norm10 = np.divide(bin10data, bin10exp, out=np.zeros_like(bin10data), where=bin10exp!=0)
    norm50 = np.divide(bin50data, bin50exp, out=np.zeros_like(bin50data), where=bin50exp!=0)
    bin10 = float(hi10count[i])-float(lo10count[i])
    bin50 = float(hi50count[i])-float(lo50count[i])
    normbin10 = norm10/bin10
    normbin50 = norm50/bin50
    saveFits(savefit+'10count_binadjusted'+det+lo10count[i]+'_'+hi10count[i]+'.fits', normbin10)
    saveFits(savefit+'50count_binadjusted'+det+lo50count[i]+'_'+hi50count[i]+'.fits', normbin50)

#    plt.matshow(normbin10, origin='lower')
#    plt.savefig(savedir+'10count_binadjusted'+det+lo10count[i]+'_'+hi10count[i]+'.png')
#    plt.matshow(normbin50, origin='lower')
#    plt.savefig(savedir+'50count_binadjusted'+det+lo50count[i]+'_'+hi50count[i]+'.png')
#    plt.close('all')
    if det == 'A':
      tot10A += normbin10
      tot50A += normbin50
    else:
      tot10B += normbin10
      tot50B += normbin50

saveFits(savefit+'10count_binadjustedtotA.fits', tot10A)
saveFits(savefit+'10count_binadjustedtotB.fits', tot10B)
saveFits(savefit+'50count_binadjustedtotA.fits', tot50A)
saveFits(savefit+'50count_binadjustedtotB.fits', tot50B)
subA10 = totnormA-tot10A; subB10 = totnormB-tot10B
subA50 = totnormA-tot50A; subB50 = totnormB-tot50B
saveFits(savefit+'Blank10countsubA.fits',subA10)
saveFits(savefit+'Blank10countsubB.fits',subB10)
saveFits(savefit+'Blank50countsubA.fits',subA50)
saveFits(savefit+'Blank50countsubB.fits',subB50)

#plt.matshow(tot10A, origin='lower')
#plt.savefig(savedir+'10count_binadjustedtotA.png')
#plt.matshow(tot10B, origin='lower')
#plt.savefig(savedir+'10count_binadjustedtotB.png')
#plt.matshow(tot50A, origin='lower')
#plt.savefig(savedir+'50count_binadjustedtotA.png')
#plt.matshow(tot50B, origin='lower')
#plt.savefig(savedir+'50count_binadjustedtotB.png')

#plt.matshow(totnormA-tot10A, origin='lower')
#plt.colorbar()
#plt.savefig(savedir+'Blank10countsubA.png')
#plt.matshow(totnormB-tot10B, origin='lower')
#plt.colorbar()
#plt.savefig(savedir+'Blank10countsubB.png')
#plt.matshow(totnormA-tot50A, origin='lower')
#plt.colorbar()
#plt.savefig(savedir+'Blank50countsubA.png')
#plt.matshow(totnormB-tot50B, origin='lower')
#plt.colorbar()
#plt.savefig(savedir+'Blank50countsubB.png')
#
#plt.close('all')
  


