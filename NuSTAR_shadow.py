import numpy as np
from Sec_Scripts import Get_RA_DEC_5deg as GET
from Sec_Scripts import fits_header as fh
from Sec_Scripts import dirs 
import os, sys, pidly
import matplotlib.pyplot as plt

idl = pidly.IDL('/uufs/chpc.utah.edu/sys/pkg/idl/8.4/idl84/bin/idl')

obs_dir = dirs('OBS')

if len(sys.argv) <= 1:
	print('You must pass an OBSID')
	sys.exit()

obsid = sys.argv[1]

cldir = obs_dir+'/'+obsid+'/event_cl'

try:
	print('Getting the coords of the Sat. pointing')
	RA_PNT, DEC_PNT, PA_PNT, EXP = fh(obsid,'A')
	print(RA_PNT, DEC_PNT, PA_PNT)
except:
	print('No A, so trying B')
	RA_PNT, DEC_PNT, PA_PNT, EXP = fh(obsid,'B')
	print(RA_PNT, DEC_PNT, PA_PNT)

print('Getting the coords. to the possible sources')
RA_src, DEC_src, Dist_src, flux_src = GET(RA_PNT, DEC_PNT, 3.5)

D_len = len(RA_src)
src_data = np.zeros((D_len,4))
for i in range(D_len):
	src_data[i][0] = RA_src[i]
	src_data[i][1] = DEC_src[i]
	src_data[i][2] = Dist_src[i]
	src_data[i][3] = flux_src[i]

print('Ordering:')
src_data = src_data*np.array([1,1,1,-1])
src_data = src_data[src_data[:,3].argsort()] # Sorts the list of possible sources by distance from the RA_PNT & DEC_PNT 
src_data = src_data*np.array([1,1,1,-1])

home = '/uufs/astro.utah.edu/common/home/u1019304/'

print('Plotting the shadow')
for i in range(D_len):
	idl('.r onesrc')
	idl('.r onesrc')
	idl('.r onesrc')
	rasrc = src_data[i][0]; decsrc = src_data[i][1]
	idl('onesrc,'+str(rasrc)+','+str(decsrc)+','+str(PA_PNT)+','+str(RA_PNT)+','+str(DEC_PNT)+',/plotit')
#	os.system('display '+home+'NuSTAR/Figures/blah.png')
	likes = raw_input('Was this a good source?  [Y/N]')
	if likes == 'Y':
		print('The info for that source is:')
		print(str(src_data[i]))
		for det in ['A','B']:
			os.system('mv '+home+'NuSTAR/Figures/blah'+det+'.fits '+home+'NuSTAR/Figures/'+obsid+'_'+det+str(i)+'.fits')
	else:
		os.system('rm '+home+'NuSTAR/Figure/blah*')
	blahzy = raw_input('Any key when ready for next image: '+str(i+1)+'/'+str(D_len)+' sources')


