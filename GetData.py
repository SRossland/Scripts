#!/uufs/astro.utah.edu/common/home/u1019304/VENV3.7.9/bin/python3

# Syntax: GetData.py

# Purpose:  To compare SUN images to NOSUN images and get the rate of difference
#           Then to be used w.r.t. the boresight angle and possibly roll....

import os, sys, numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

###############
# Definitions:

def get_data(fi):
    with fits.open(fi) as hdul:
        dat = hdul[0].data
        exp = hdul[0].header['EXPOSURE']
    return dat, exp

def make_fits_file(fi,ary,ex,BSA_low,BSA_high,error):
    fits.writeto(fi,np.zeros((360,360)))
    with fits.open(fi,mode='update') as hdul:
        hdul[0].data = ary
        hdul[0].header['EXPOSURE'] = ex
        hdul[0].header['COMMENT'] = 'Residual image'
        hdul[0].header['COMMENT'] = 'Boresight angle between: {} {}'.format(BSA_low,BSA_high)
        hdul[0].header['ERROR'] = error
        hdul.flush()

def get_exposure(fi):
    with fits.open(fi) as hdul:
        exp = hdul[0].header['EXPOSURE']
    return exp


###############

# main

# assignements

newmaskA = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/SolarData/fullmaskA_final.fits'
newmaskB = newmaskA.replace('maskA','maskB')
with fits.open(newmaskA) as hdul:
    maskA = hdul[0].data
with fits.open(newmaskB) as hdul:
    maskB = hdul[0].data

# hard coding the boresight angles


angles = ['0','60','75','90','105','120','135','150','180']

user_dir = '/uufs/astro.utah.edu/common/home/u1019304/'
#Sun_BSA_fits_dir = os.path.join(user_dir,'NuSTAR/SolarData/BSA_strip_fits')
#Nosun_BSA_fits_dir = os.path.joini(Sun_BSA_fits_dir,'nosun_fits')
#Sun_split_fits_dir = Sun_BSA_fits_dir.replace('strip','split_stript')
#Nosun_split_fits_dir = os.path.join(Sun_split_fits_dir,'nosun_fits')

# Assuming this is ran from the sun file directory
directory = os.getcwd()

# so I am going to grab the Sun data, then look for the cooresponding nosun data
file_list = os.listdir(directory)

sun_fits_files = []
rl_sun_fits, rr_sun_fits = [],[]
for li in file_list:
    if 'Data' in li:
      if '3_5' in li:
        if 'RL.fits' in li:
            rl_sun_fits.append(li)
            continue
        if 'RR.fits' in li:
            rr_sun_fits.append(li)
            continue
        else:
            sun_fits_files.append(li)

nosun_dir = os.path.join(directory,'nosun_fits')

lists_of_fits = [sun_fits_files,rr_sun_fits,rl_sun_fits]

for fits_list in lists_of_fits:
    if len(fits_list) == 0: 
        continue
    for sun_file in fits_list:
        det = sun_file.split('_')[0][4]
        if det == 'A':
            mask = newmaskA
        else:
            mask = newmaskB
        nosun_file = sun_file.replace('_sun_','_nosun_')
        sun_data, sun_exp = get_data(sun_file)
        nosun_data, nosun_exp = get_data('./nosun_fits/{}'.format(nosun_file))
        
        with fits.open(mask) as hdul:
            new_mask = hdul[0].data

        sun_exp_array = np.ones((360,360))*sun_exp*new_mask
        nosun_exp_array = np.ones((360,360))*nosun_exp*new_mask

        sun_data *= new_mask
        nosun_data *= new_mask

        total_count_error_sun = np.sqrt(np.sum(sun_data))
        total_count_error_nosun = np.sqrt(np.sum(nosun_data))

        normalized_sun = np.copy(sun_data)
        normalized_sun = np.divide(normalized_sun,sun_exp_array,out=np.zeros_like(normalized_sun),where=sun_exp_array!=0)


        normalized_nosun = np.copy(nosun_data)
        normalized_nosun = np.divide(normalized_nosun,nosun_exp_array,out=np.zeros_like(normalized_nosun),where=nosun_exp_array!=0)

        norm_count_error_sun = np.sqrt(np.sum(normalized_sun))
        norm_count_error_nosun = np.sqrt(np.sum(normalized_nosun))

        residual_image = normalized_sun - normalized_nosun
        write_directory = os.path.join(directory,'residual_fits')
        file_write_name = sun_file.replace('Data','Resid')
        file_write_name = file_write_name.replace('_sun_','_sun_less_nosun_')
        fi = '{}/{}'.format(write_directory,file_write_name)

        fits.writeto(fi, residual_image)
        with fits.open(fi, mode='update') as hdul:
            hdul[0].header['TOTERRSUN'] = total_count_error_sun
            hdul[0].header['TOTERNOSU'] = total_count_error_nosun
            hdul[0].header['NORMERSUN'] = norm_count_error_sun
            hdul[0].header['NORMERNOS'] = norm_count_error_nosun
            hdul.flush()


















