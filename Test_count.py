#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntax: Test_count.py

# Purpose:  This code creates a false image for NuSTAR to test 
# minimization techniques with poissonian noise.  All numbers 
# will be hard coded in and placed so one can change easily

# Types of images created:
#	1) Flat
# 	2) Flat ap, with DET variance
#	3) Gradiant Ap
#	4) Gradiant Ap with DET variance


import numpy as np
from astropy.io import fits

################
# Constants
Aij_con = 0.000002
Bij_con = 0.000003
A0_con = 1.0
A1_con = 2.1
A2_con = 3.2
A3_con = 4.6
B0_con = 1.0
B1_con = 1.0
B2_con = 1.0
B3_con = 2.0

################
#  Masks

M_dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil'

# Edge masks

edgeA_file = M_dir+'/edgemaskA.fits'
edgeB_file = M_dir+'/edgemaskB.fits'

# det masks
# A
det0Amask_file = M_dir+'/det0mask.fits'
det1Amask_file = M_dir+'/det1mask.fits'
det2Amask_file = M_dir+'/det2mask.fits'
det3Amask_file = M_dir+'/det3mask.fits'

# B
det0Bmask_file = M_dir+'/detB0mask.fits'
det1Bmask_file = M_dir+'/detB1mask.fits'
det2Bmask_file = M_dir+'/detB2mask.fits'
det3Bmask_file = M_dir+'/detB3mask.fits'

# Ap 
apAmask_file = '/uufs/astro.utah.edu/common/home/u1019304/my_idl/nuskybgd/auxil/detA_det1.img'
apBmask_file = '/uufs/astro.utah.edu/common/home/u1019304/my_idl/nuskybgd/auxil/detB_det1.img'

################
# Functions

def getfits(file, i):
    with fits.open(file) as hdul:
        var = hdul[i].data
    return var

def noise(img):
    nois = np.random.poisson(img).astype(float)
    return nois 

def writefits(array, filename, dire):
    hdufile = fits.PrimaryHDU(array)
    hdufile.writeto(dire+'/'+filename+'.fits')
    
################

edgeA = getfits(edgeA_file,0)
edgeB = getfits(edgeB_file,0) 
det0A = getfits(det0Amask_file,0)
det1A = getfits(det1Amask_file,0)
det2A = getfits(det2Amask_file,0)
det3A = getfits(det3Amask_file,0)
det0B = getfits(det0Bmask_file,0)
det1B = getfits(det1Bmask_file,0)
det2B = getfits(det2Bmask_file,0)
det3B = getfits(det3Bmask_file,0)
apAma = getfits(apAmask_file,0)
apBma = getfits(apBmask_file,0)

################
dataA_flat_raw = np.ones((360,360))
dataB_flat_raw = np.ones((360,360))

dataA_flat_detvar_raw = np.ones((360,360))
dataB_flat_detvar_raw = np.ones((360,360))

dataA_apgrad_raw = np.ones((360,360))
dataB_apgrad_raw = np.ones((360,360))

dataA_apgrad_detvar_raw = np.ones((360,360))
dataB_apgrad_detvar_raw = np.ones((360,360))

Aij = apAma[len(apAma)/2 - 180:len(apAma)/2 + 180, len(apAma)/2 - 180:len(apAma)/2 + 180]*edgeA
Bij = apBma[len(apBma)/2 - 180:len(apBma)/2 + 180, len(apBma)/2 - 180:len(apBma)/2 + 180]*edgeB
#################

dataA_flat = (dataA_flat_raw*edgeA*Aij_con)
dataB_flat = (dataB_flat_raw*edgeB*Bij_con)

dataA_flat_detvar = noise((dataA_flat_detvar_raw+det0A*A0_con+det1A*A1_con+det2A*A2_con+det3A*A3_con)*edgeA)
dataB_flat_detvar = noise((dataB_flat_detvar_raw+det0B*B0_con+det1B*B1_con+det2B*B2_con+det3B*B3_con)*edgeB)

dataA_apgrad = noise(Aij*edgeA*Aij_con)
dataB_apgrad = noise(Bij*edgeB*Bij_con)

dataA_apgrad_detvar = noise((Aij*Aij_con+det0A*A0_con+det1A*A1_con+det2A*A2_con+det3A*A3_con)*edgeA)
dataB_apgrad_detvar = noise((Bij*Bij_con+det0B*B0_con+det1B*B1_con+det2B*B2_con+det3B*B3_con)*edgeB)

writefits(dataA_flat,'dataA_flat',M_dir)
writefits(dataB_flat,'dataB_flat',M_dir)
writefits(dataA_flat_detvar,'dataA_flat_detvar',M_dir)
writefits(dataB_flat_detvar,'dataB_flat_detvar',M_dir)
writefits(dataA_apgrad,'dataA_apgrad',M_dir)
writefits(dataB_apgrad,'dataB_apgrad',M_dir)
writefits(dataA_apgrad_detvar,'dataA_apgrad_detvar',M_dir)
writefits(dataB_apgrad_detvar,'dataB_apgrad_detvar',M_dir)
