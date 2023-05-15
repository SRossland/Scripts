import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib
from astropy.io import fits
from astropy.coordinates import SkyCoord as sc
from astropy import units as u

#obs_dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs/'
#file1 = 'Full_list_0_and_lt85_source.dat'
file1 = 'Zero_obs.txt'
file2 = 'src_obs.txt'
file3 = 'Largesrc_obs.txt'

#list_src = np.loadtxt(obs_dir+file1)
list_src = np.loadtxt(file1)
list_src2 = np.loadtxt(file2)
list_src3 = np.loadtxt(file3)

obs_full = []
for i in range(len(list_src)):
	obs_full.append(str(int(list_src[i])))
for i in range(len(list_src2)):
	obs_full.append(str(int(list_src2[i])))
for i in range(len(list_src3)):
	obs_full.append(str(int(list_src3[i])))

file_dir = '/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/xray/NuSTAR/OBS'

RA = np.zeros(len(obs_full))
Dec = np.zeros(len(obs_full))
exp = np.zeros(len(obs_full))

for i in range(len(obs_full)):
	obsid = obs_full[i][0:11]
	eventA = file_dir+'/'+obsid+'/event_cl/nu'+obsid+'A01_cl.evt'
	eventB = file_dir+'/'+obsid+'/event_cl/nu'+obsid+'B01_cl.evt'
	if os.path.isfile(eventA) == True:
		event = file_dir+'/'+obsid+'/event_cl/nu'+obsid+'A01_cl.evt'
	elif os.path.isfile(eventB) == True:
		event = file_dir+'/'+obsid+'/event_cl/nu'+obsid+'B01_cl.evt'	
	else:
		continue
		
	with fits.open(event) as hdul:
		if hdul[1].header['EXPOSURE'] <= 0.0:
			continue	
		exp[i] = hdul[1].header['EXPOSURE']			
		RA[i] = hdul[0].header['RA_PNT']
		Dec[i] = hdul[0].header['DEC_PNT']
#	print(exp[i])

idx = np.where(exp == 0.0)
print(idx)
#print(len(obs_full))
cmap = matplotlib.cm.get_cmap('viridis_r')
#cmap = matplotlib.cm.get_cmap('hot_r')
normalize = matplotlib.colors.LogNorm(vmin=7000., vmax=max(exp))

colors = [cmap(normalize(value)) for value in exp]

#print(RA[0:10])
#print(Dec[0:10])
c = sc(ra = RA*u.degree, dec = Dec*u.degree, frame='fk5')
c, c.galactic

#with open('OBS_pointing_and_exp_list.ps','a+') as O:
#        O.write("{:^11} | {:^9} | {:^9} | {:^11}".format('OBSID','RA','DEC','EXP'))

#for i in range(len(exp)):
#        with open('OBS_pointing_and_exp_list.ps','a+') as O:
#		ci = sc(ra = RA[i]*u.degree, dec = Dec[i]*u.degree, frame='fk5')
#                O.write("{:<11} | {:^9} | {:^9} | {:11}".format(obs_full[i][0:11],RA[i],Dec[i],exp[i]))
#		if abs(ci.galactic.b.deg) < 10.0:
#                        O.write('*Not in stack')
#                O.write("\n")

tick_labels = np.array([150,120,90,60,30,0,330,300,270,240,210])
tick_labels = np.remainder(tick_labels+360+0,360)

x = np.remainder(c.galactic.l.degree+360,360)
ind = x>180
x[ind] -=360
x = -x

alpha = 0.25

fig = plt.figure(figsize=(16,12))
fig.patch.set_facecolor('White')
ax = fig.add_subplot(111, projection="mollweide", axisbg = "White")
#ax.scatter(c.galactic.l.radian, c.galactic.b.radian, color=colors, alpha=alpha)
ax.scatter(np.radians(x),np.radians(c.galactic.b.degree),s = 20, linewidth=8, color='white' ,edgecolor=colors, alpha=alpha)
#ax.scatter(np.radians(x),np.radians(c.galactic.b.degree),color=colors,alpha=alpha)
ax.axhspan((-10.0/180.0)*np.pi,(10.0/180.0)*np.pi, alpha=0.5, color='red')
ax.set_xticklabels(tick_labels)
ax.set_xlabel("RA")
ax.set_ylabel("DEC")
ax.grid(True)
cax, _ = matplotlib.colorbar.make_axes(ax)
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)
cbar.set_label("Log Time of Exposure")
plt.show()


