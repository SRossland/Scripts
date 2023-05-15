#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntax: GetCstat.py path/to/directory/with/files

# Purpose: This script takes the values of the given norms from both the fit statistic program 
#          and the mcmc program and compares them to ensure that the mcmc is working correctly.

import os, sys, numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# sys.argv assignments:

data_fileA_aCXB = sys.argv[1]+'/A_params_mcmc.txt'
data_fileB_aCXB = sys.argv[1]+'/B_params_mcmc.txt'
data_fileA_det = sys.argv[1]+'/Adetnormerrors_mcmc.txt'
data_fileB_det = sys.argv[1]+'/Bdetnormerrors_mcmc.txt'


# This path is hard coded in since all files are found in this directory
xray_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/CXB/01/nosun/'
home_dir = '/uufs/astro.utah.edu/common/home/u1019304/'
mask_fileA = home_dir+'temp/fullmaskA_final.fits'
mask_fileB = home_dir+'temp/fullmaskB_final.fits'
grA = home_dir+'my_idl/nuskybgd/auxil/detA_det1.img'
grB = home_dir+'my_idl/nuskybgd/auxil/detB_det1.img'
detA = home_dir+'temp/detmapA.fits'
detB = home_dir+'temp/detmapB.fits'


############ FUNCTIONS	##################

def getData(fi):
  with fits.open(fi) as hdul:
    var = hdul[0].data
  return var

def Cstat(p,A,det0,det1,det2,det3,count):
  x_c, y_c, z_c, w_c, t_c = p
  model = ((x_c*A)+(y_c*det0)+(z_c*det1)+(w_c*det2)+(t_c*det3))
  numerA = model - count + count*(np.log(count) - np.log(model))
  lp = 2*np.sum(numerA)
  return lp

def floatList(li):
  ret = [float(i) for i in li]
  return ret

##########################################

maskA = getData(mask_fileA)
maskB = getData(mask_fileB)
gradA = getData(grA)
gradB = getData(grB)
detmapA = getData(detA)
detmapB = getData(detB)

# I can pull the Cstat values from the brute force method from the longformat files

# I need to recreate the mcmc Cstat values

lowe, highe = [],[]
aCXBA,aCXBB = [],[] # mcmc values for aCXB telescope A and B


with open(data_fileA_aCXB,'r+') as f:
  lines = f.readlines()
for line in lines:
  lowe.append(line.split()[0])
  highe.append(line.split()[1])
  aCXBA.append(line.split()[2])
with open(data_fileB_aCXB,'r+') as f:
  lines = f.readlines()
for line in lines:
  aCXBB.append(line.split()[2])

det0A,det1A,det2A,det3A = [],[],[],[] # detector values for the mcmc run
det0B,det1B,det2B,det3B = [],[],[],[]

with open(data_fileA_det,'r+') as f:
  lines = f.readlines()
for line in lines:
  det0A.append(line.split()[2])
  det1A.append(line.split()[5])
  det2A.append(line.split()[8])
  det3A.append(line.split()[11])

with open(data_fileB_det,'r+') as f:
  lines = f.readlines()
for line in lines:
  det0B.append(line.split()[2])
  det1B.append(line.split()[5])
  det2B.append(line.split()[8])
  det3B.append(line.split()[11])

aCXBA = floatList(aCXBA)
aCXBB = floatList(aCXBB)
det0A = floatList(det0A)
det1A = floatList(det1A)
det2A = floatList(det2A)
det3A = floatList(det3A)
det0B = floatList(det0B)
det1B = floatList(det1B)
det2B = floatList(det2B)
det3B = floatList(det3B)


# Now i need the bins, gradiant, and count

CstatA_mcmc, CstatB_mcmc = [],[]

# out of loop files
# the expmaps will be found in the fits_file directory and will be the same for all energy bins. Find one.

expmapA = sys.argv[1]+'/fits_file/ExpA_01_nosun_'+lowe[0]+'_'+highe[0]+'keV.fits' 
expmapB = sys.argv[1]+'/fits_file/ExpB_01_nosun_'+lowe[0]+'_'+highe[0]+'keV.fits'

with fits.open(expmapA) as hdul:
  exposureA = hdul[0].header['EXPOSURE']
with fits.open(expmapB) as hdul:
  exposureB = hdul[0].header['EXPOSURE']

expA = getData(expmapA)
expB = getData(expmapB)

detmapA = detmapA * maskA
detmapB = detmapB * maskB

expA = expA * maskA
expB = expB * maskB

Aij = gradA[int(len(gradA)/2) - 180:int(len(gradA)/2) + 180, int(len(gradA)/2) - 180:int(len(gradA)/2) + 180]*maskA
Bij = gradB[int(len(gradB)/2) - 180:int(len(gradB)/2) + 180, int(len(gradB)/2) - 180:int(len(gradB)/2) + 180]*maskB

Aij *= 0.0001
Bij *= 0.0001

expmaskA = np.copy(expA)
expmaskB = np.copy(expB)
expmaskA /= exposureA
expmaskB /= exposureB

a0 = np.zeros((360,360)); a1 = np.zeros((360,360)); a2 = np.zeros((360,360)); a3 = np.zeros((360,360))
b0 = np.zeros((360,360)); b1 = np.zeros((360,360)); b2 = np.zeros((360,360)); b3 = np.zeros((360,360))

a0[detmapA == 0] = 1
a1[detmapA == 1] = 1
a2[detmapA == 2] = 1
a3[detmapA == 3] = 1

b0[detmapB == 0] = 1
b1[detmapB == 1] = 1
b2[detmapB == 2] = 1
b3[detmapB == 3] = 1

Aij *= expmaskA; a0 *= expmaskA; a1 *= expmaskA; a3 *= expmaskA; #a2 *= expmaskA
Bij *= expmaskB; b0 *= expmaskB; b1 *= expmaskB; b3 *= expmaskB; #b2 *= expmaskB

for i in range(len(lowe)):

  data_file_A = sys.argv[1]+'/fits_file/DataA_01_nosun_'+lowe[i]+'_'+highe[i]+'keV.fits'      #fits_file directory with the name DataA_01_nosun_blahblahblah.fits  
  data_file_B = sys.argv[1]+'/fits_file/DataB_01_nosun_'+lowe[i]+'_'+highe[i]+'keV.fits'

  dataA = getData(data_file_A)
  dataB = getData(data_file_B)

  bins_file_A = data_file_A.replace('Data','BINS')
  bins_file_B = data_file_B.replace('Data','BINS')

  binsA = getData(bins_file_A)
  binsB = getData(bins_file_B)

  A_bins,det0_bins,det1_bins,det2_bins,det3_bins,e_m,count_bins = [],[],[],[],[],[],[]
  for j in range(int(np.max(binsA))+1): 
    idx = np.where(binsA == j)
    A_bins.append(np.sum(Aij[idx]))
    e_m.append(np.sum(expmaskA[idx]))
    count_bins.append(np.sum(dataA[idx]))
    det0_bins.append(np.sum(a0[idx]))
    det1_bins.append(np.sum(a1[idx]))
    det2_bins.append(np.sum(a2[idx]))
    det3_bins.append(np.sum(a3[idx]))

  p0 = [aCXBA[i],det0A[i],det1A[i],det2A[i],det3A[i]]

  if count_bins[-1] == 0:
    for matrix in [A_bins,det0_bins,det1_bins,det2_bins,det3_bins,e_m,count_bins]:
      var1 = matrix.pop()
      var2 = matrix.pop()
      matrix.append(var1+var2)

  A_bins = np.asarray(A_bins)
  det0_bins = np.asarray(det0_bins)
  det1_bins = np.asarray(det1_bins)
  det2_bins = np.asarray(det2_bins)
  det3_bins = np.asarray(det3_bins)
  e_m = np.asarray(e_m)
  count_bins = np.asarray(count_bins)

  CstatA_mcmc.append(Cstat(p0,A_bins,det0_bins,det1_bins,det2_bins,det3_bins,count_bins))
  
  B_bins,dB0,dB1,dB2,dB3,emB,countB = [],[],[],[],[],[],[]
  for j in range(int(np.max(binsB))+1):
    idx = np.where(binsB == j)
    B_bins.append(np.sum(Bij[idx]))
    emB.append(np.sum(expmaskB[idx]))
    countB.append(np.sum(dataB[idx]))
    dB0.append(np.sum(b0[idx]))
    dB1.append(np.sum(b1[idx]))
    dB2.append(np.sum(b2[idx]))
    dB3.append(np.sum(b3[idx]))

  if countB[-1] == 0:
    for matrix in [B_bins,dB0,dB1,dB2,dB3,emB,countB]:
      var1 = matrix.pop()
      var2 = matrix.pop()
      matrix.append(var1+var2)

  B_bins = np.asarray(B_bins)
  dB0 = np.asarray(dB0)
  dB1 = np.asarray(dB1)
  dB2 = np.asarray(dB2)
  dB3 = np.asarray(dB3)
  emB = np.asarray(emB)
  countB = np.asarray(countB)

  p0B = [aCXBB[i],det0B[i],det1B[i],det2B[i],det3B[i]]
  CstatB_mcmc.append(Cstat(p0B,B_bins,dB0,dB1,dB2,dB3,countB))


# get the brute force Cstat values:

CstatA_brute, CstatB_brute = [],[]
Xvals = []

for i in range(len(lowe)):
  Xvals.append((float(lowe[i])+float(highe[i]))/2.)
  for det in ['A','B']:
    fi = sys.argv[1]+'/'+det+'_nosun_'+lowe[i]+'_'+highe[i]+'_keVparams_longformat.txt'
    with open(fi,'r+') as f:
      lines = f.readlines()
    for line in lines:
      if 'Min Cstat value:' in line:
        if det == 'A':
          CstatA_brute.append(float(line.split()[3]))
        if det == 'B':
          CstatB_brute.append(float(line.split()[3]))


plt.scatter(Xvals,CstatA_mcmc,c='blue',s=4)
plt.scatter(Xvals,CstatA_brute,c='red',s=4)
plt.show()

plt.scatter(Xvals,CstatB_mcmc,c='blue',s=4)
plt.scatter(Xvals,CstatB_brute,c='red',s=4)
plt.show()




