#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntax: count_stat_raw.py det listname

# This script is similar to count_stat.py, but it looks at the nu#02_cl.evt files instead 
# of an image and uses the data itself for the energy that you specify

import os, sys, numpy as np
from astropy.io import fits
import scipy.optimize as opt

##################
# Directories

edgeA = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/fullmaskA.fits'
edgeB = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/fullmaskB.fits'

#detnumA = '/uufs/astro.utah.edu/common/home/u1019304/my_idl/nuskybgd/auxil/detnumA.fits'
#detnumB = '/uufs/astro.utah.edu/common/home/u1019304/my_idl/nuskybgd/auxil/detnumB.fits'
#pixmapA = '/uufs/astro.utah.edu/common/home/u1019304/my_idl/nuskybgd/auxil/pixmapA.fits'
#pixmapB = '/uufs/astro.utah.edu/common/home/u1019304/my_idl/nuskybgd/auxil/pixmapB.fits'

xray_dir = '/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/xray/NuSTAR/OBS'
list_dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs'

writedir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/spec'
##################

det = sys.argv[1]
lst = sys.argv[2]

#lst = raw_input('List: ')

#el = float(sys.argv[3])
#hl = float(sys.argv[4])

el = eval(raw_input('Low Energy Limit:  '))
hl = eval(raw_input('High Energy Limit: '))

##################
# Functions:
def energy(Channel):
    return Channel*0.04 + 1.6

def chan(Ener):
    return (Ener-1.6)/0.04

def getfits(fi): # need to write this to do the full file, sun, nosun
    with fits.open(fi) as hdul:
        PI = hdul[1].data['PI']
        DETX = hdul[1].data['DET1X']
        DETY = hdul[1].data['DET1Y']
    return PI,DETX,DETY
    
def writefits(data,filename):
    global writedir
    hdu = fits.PrimaryHDU(data)
    hdu.writeto(writedir+'/'+filename)
    

##################
# Get obsid's


f = open(list_dir+'/'+lst,'r+')
obs = []
for line in f.readlines():
    obs.append(line.rstrip())

    
##################
# Get masks and other standard arrays to be used

with fits.open(edgeA) as hdul_ea:
    ed_maskA = hdul_ea[0].data
#with fits.open(edgeB) as hdul_eb:
#    ed_maskB = hdul_eb[0].data
dataA = np.zeros((360,360)) # this will be the stacked data for A
dataA_SUN = np.zeros((360,360))
dataA_NOSUN = np.zeros((360,360))
h = 0
###################
el_chan = chan(el)
hl_chan = chan(hl)
print(el_chan,hl_chan)

for i in range(len(obs)):
#    for detector in ['A']:
#    for detector in ['A','B']:
        for event in ['full','SUN','NOSUN']:

            if event == 'full':
                file_events = xray_dir+'/'+obs[i]+'/event_cl/nu'+obs[i]+det+'02_cl.evt'
		maps = dataA
            elif event == 'SUN':
                file_events = xray_dir+'/'+obs[i]+'/event_sep_cl/OCC/nu'+obs[i]+det+'02_fullevts_'+event+'.fits'
		maps = dataA_SUN
            elif event == 'NOSUN':
                file_events = xray_dir+'/'+obs[i]+'/event_sep_cl/OCC/nu'+obs[i]+det+'02_fullevts_'+event+'.fits'
            	maps = dataA_NOSUN

            if os.path.isfile(file_events) == True:
                PI,DETX,DETY = getfits(file_events)
            else:
                continue
                
#            cldir = xray_dir+'/'+obs[i]+'/event_cl'
#            instrmap = cldir+'/bgd/newinstrmap'+detector+'.fits'
#        
#            with fits.open(instrmap) as hdul_I:
#                intmap = hdul_I[1].data
#            intmap[intmap>0] = 1
        
        # Here I will create a full mask for each obs 
        # This includes the edge mask, inturment mask, setup detmasks
        
            #if detector == 'A':
             #   if event == 'full':
             #       maps = dataA
             #   elif event == 'SUN':
             #       maps = dataA_SUN
             #   elif event == 'NOSUN':
             #       maps = dataA_NOSUN
        # Now to use the masks with the data
        # I think for the expmap, it will be easiest to do a full map, then subtract each det exp at the end
            #maps[1] += exp
            for j in range(len(PI)):
                # see if the PI is in the energy range specified and if the it is not masked
                #PI_en = energy(PI[j])
		if PI[j] <= -1:
		    continue
                elif (DETY[j] == -4095) or (DETX[j] == -4095):
                    continue
                elif (PI[j] >= int(el_chan)) and (PI[j] <= int(hl_chan)) and (ed_maskA[DETY[j]-1,DETX[j]-1] == 1):
                    maps[DETY[j]-1,DETX[j]-1] += 1
	#	    dataA_NOSUN[DETY[j]-1,DETX[j]-1] += 1
                else:
                    continue

datapack = [dataA,dataA_SUN,dataA_NOSUN]
#            dataB_NOSUN,expmapA,expmapA_SUN,expmapA_NOSUN\
#            ,expmapB,expmapB_SUN,expmapB_NOSUN]
names = ['data'+det,'data'+det+'_SUN','data'+det+'_NOSUN']
#            'dataB_NOSUN','expmapA','expmapA_SUN','expmapA_NOSUN'\
#            ,'expmapB','expmapB_SUN','expmapB_NOSUN']
#datapack = [dataA,dataA_SUN,dataA_NOSUN]
#names = ['dataA','dataA_SUN','dataA_NOSUN']
#writefits(dataA_NOSUN,'dataA_test.fits')

#for fil in range(12):
for fil in range(3):
    writefits(datapack[fil],names[fil]+'.fits')
    
