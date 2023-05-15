#!/uufs/astro.utah.edu/common/home/u1019304/VENV3.7.9/bin/python3

import os, sys, numpy as np
from astropy.io import fits

Solar_dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/SolarData/'

write_dir = os.path.join(Solar_dir,'analysis/')

def get_data(fi,lvl):
    with fits.open(fi) as hdul:
        dat = hdul[lvl].data
    return dat

def get_exp(fi):
    with fits.open(fi) as hdul:
        exposure = hdul[0].header['EXPOSURE']
    return exposure

def grid_space_value(arr,grid):
    values = np.zeros(9)
    #avg_grid_val = get_average_val(grid)
    for i in range(9):
        idx = np.where(grid == i+1)
        pix_count = len(idx[0])
        values[i] = np.sum(arr[idx])/pix_count
        #values[i] = val*avg_grid_val
    return values

def nosun_grid_value(arr,grid,exp):
    values = np.zeros(9)
    errors = np.zeros(9)
    for i in range(9):
        idx = np.where(grid == i+1)
        pix_count = len(idx[0])
        val = np.sum(arr[idx])/exp
        err = np.sqrt(np.abs(np.sum(arr[idx])))/exp
        val /= pix_count
        err /= pix_count
        values[i] = val
        errors[i] = err
    return values, errors

def grid_errors(sun_img,nosun_img,grid,exp_sun,exp_nosun):
    errors = np.zeros(9)
    #avg_grid_val = get_average_val(grid)
    for i in range(9):
        idx = np.where(grid == i+1)
        pix_count = len(idx[0])
        sun_val = np.sqrt(np.abs(np.sum(sun_img[idx])))/exp_sun
        sun_val /= pix_count
        #sun_val *= avg_grid_val
        if exp_nosun == 0.:
            nosun_val = 0.
        else:
            nosun_val = np.sqrt(np.abs(np.sum(nosun_img[idx])))/exp_nosun
            nosun_val /= pix_count
            #nosun_val *= avg_grid_val
        errors[i] = np.sqrt(sun_val**2+nosun_val**2)
    return errors

def get_average_val(arr):
    vals = []
    for i in range(9):
        idx = np.where(arr == i+1)
        vals.append(len(idx[0]))
    return np.average(vals)

angles = ['0','60','75','90','105','120','135','150','180']

gridA_file = os.path.join(Solar_dir,'grid6x6A.fits')
gridB_file = gridA_file.replace('6x6A','6x6B')

maskA_file = os.path.join(Solar_dir,'fullmaskA_final.fits')
maskB_file = maskA_file.replace('A_final','B_final')

detnumA_file = os.path.join(Solar_dir,'detnumA.fits')
detnumB_file = detnumA_file.replace('A.fits','B.fits')

gridA = get_data(gridA_file,0)
gridB = get_data(gridB_file,0)
mA = get_data(maskA_file,0)
mB = get_data(maskB_file,0)
detA = get_data(detnumA_file,1)
detB = get_data(detnumB_file,1)

BSA_dir = os.path.join(Solar_dir,'BSA_strip_fits/residual_fits/')
exp_dir = os.path.join(Solar_dir,'BSA_strip_fits/')
nosun_dir = os.path.join(exp_dir,'nosun_fits/')

for det in ['A','B']:

    if det == 'A':
        grids = gridA
        mask = mA
        detnum = detA
    else:
        grids = gridB
        mask = mB
        detnum = detB

    for k in range(8):
        event = '{}Resid{}_02_sun_less_nosun_3_5_{}_{}.fits'.format(BSA_dir,det,angles[k],angles[k+1])
        event_data = get_data(event,0)
        exp_event = '{}Data{}_02_sun_3_5_{}_{}.fits'.format(exp_dir,det,angles[k],angles[k+1])        
        exp_val = get_exp(exp_event)
        sun_data = get_data(exp_event,0)
        sun_data *= mask
        nosun_event = '{}Data{}_02_nosun_3_5_{}_{}.fits'.format(nosun_dir,det,angles[k],angles[k+1]) 
        nosun_data = get_data(nosun_event,0)
        nosun_exp = get_exp(nosun_event)
        nosun_data *= mask
        for i in range(4):
            a0 = np.zeros((360,360))
            a0[detnum==i] = 1
            subgrid = np.copy(grids)
            subgrid *= a0
            subgrid *= mask
            if nosun_exp == 0.:
                nosun_norm_values = np.zeros(9)
                nosun_norm_values = np.zeros(9)
            else:
                nosun_norm_values, nosun_norm_error = nosun_grid_value(nosun_data, subgrid, nosun_exp)
            vals = grid_space_value(event_data,subgrid)
            errs = grid_errors(sun_data,nosun_data,subgrid,exp_val,nosun_exp)
# TEMP TO WRITE NOSUN VALUE
#            with open(write_dir+'Solar_det_values_BSA.txt','a+') as O:
#                O.write('DET{} {} {} {} {} '.format(i,angles[k],angles[k+1],exp_val,nosun_exp))
#                for i in range(9):
#                    O.write('{} {} '.format(vals[i],errs[i]))
#                O.write('{}'.format('\n'))
            with open(write_dir+'NOSUN_det_values_BSA.txt','a+') as O:
                O.write('DET{} {} {} {}'.format(i,angles[k], angles[k+1], nosun_exp))
                for j in range(9):
                    O.write(' {} {} '.format(nosun_norm_values[j], nosun_norm_error[j]))
                O.write('\n')
                #O.write('{}{}'.format(vals,'\n'))
                #O.write('{}{}'.format(errs,'\n'))


#BSA_dir = BSA_dir.replace('strip_','split_strip_')
#exp_dir = exp_dir.replace('strip_','split_strip_')
#nosun_dir = nosun_dir.replace('strip_','split_strip_')

#for rol in ['RR','RL']:
#  for det in ['A','B']:

#    if det == 'A':
#        grids = gridA
#        mask = mA
#        detnum = detA
#    else:
#        grids = gridB
#        mask = mB
#        detnum = detB

#    for k in range(8):
#        event = '{}Resid{}_02_sun_less_nosun_3_5_{}_{}{}.fits'.format(BSA_dir,det,angles[k],angles[k+1],rol)
#        event_data = get_data(event,0)
#        exp_event = '{}Data{}_02_sun_3_5_{}_{}{}.fits'.format(exp_dir,det,angles[k],angles[k+1],rol)
#        exp_val = get_exp(exp_event)
#        sun_data = get_data(exp_event,0)
#        sun_data *= mask
#        nosun_event = '{}Data{}_02_nosun_3_5_{}_{}{}.fits'.format(nosun_dir,det,angles[k],angles[k+1],rol)
#        nosun_data = get_data(nosun_event,0)
#        nosun_exp = get_exp(nosun_event)
#        nosun_data *= mask
#        for i in range(4):
#            a0 = np.zeros((360,360))
#            a0[detnum==i] = 1
#            subgrid = np.copy(grids)
#            subgrid *= a0
#            subgrid *= mask
#            vals = grid_space_value(event_data,subgrid)
#            errs = grid_errors(sun_data,nosun_data,subgrid,exp_val,nosun_exp)
#            with open(write_dir+'Solar_det_values_BSA_split_'+rol+'.txt','a+') as O:
#                O.write('DET{} {} {} {} {} '.format(i,angles[k],angles[k+1],exp_val,nosun_exp))
#                for i in range(9):
#                    O.write('{} {} '.format(vals[i],errs[i]))
#                O.write('{}'.format('\n'))
                #O.write('{}{}'.format(vals,'\n'))
                #O.write('{}{}'.format(errs,'\n'))










