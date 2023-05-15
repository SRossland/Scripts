#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os, sys
import matplotlib
import healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord as sc
from astropy import units as u
from pylab import cm


obslist = sys.argv[1]
ski = sys.argv[2]

if ski.lower() in ['yes','true']:
    skip = True
else:
    skip = False

obs_full = np.genfromtxt(obslist,dtype = '|U11')

#obs_dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs/'
#file1 = 'Full_list_0_and_lt85_source.dat'
#file1 = 'Zero_obs.txt'
#file2 = 'src_obs.txt'
#file3 = 'Largesrc_obs.txt'

#list_src = np.loadtxt(obs_dir+file1)
#list_src = np.loadtxt(file1)
#list_src2 = np.loadtxt(file2)
#list_src3 = np.loadtxt(file3)

#obs_full = []
#for i in range(len(list_src)):
#	obs_full.append(str(int(list_src[i])))
#for i in range(len(list_src2)):
#	obs_full.append(str(int(list_src2[i])))
#for i in range(len(list_src3)):
#	obs_full.append(str(int(list_src3[i])))

file_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS'

RA = np.zeros(len(obs_full))
Dec = np.zeros(len(obs_full))
exp = np.zeros(len(obs_full))

for i,ob in enumerate(obs_full):
    #obsid = obs_full[i][0:11]
    eventA = file_dir+'/'+ob+'/event_cl/nu'+ob+'A01_cl.evt'
    eventB = file_dir+'/'+ob+'/event_cl/nu'+ob+'B01_cl.evt'
    if skip:
        if os.path.isfile(file_dir+'/'+ob+'/event_defcl/excl.reg'):
            continue
    if os.path.isfile(eventA) == True:
        event = file_dir+'/'+ob+'/event_cl/nu'+ob+'A01_cl.evt'
    elif os.path.isfile(eventB) == True:
        event = file_dir+'/'+ob+'/event_cl/nu'+ob+'B01_cl.evt'	
    else:
        continue

    with fits.open(event) as hdul:
        if hdul[1].header['EXPOSURE'] <= 0.0:
            continue
        exp[i] = hdul[1].header['EXPOSURE']
        RA[i] = hdul[0].header['RA_PNT']
        Dec[i] = hdul[0].header['DEC_PNT']
#	print(exp[i])
##############################################################################################
# New shit

nside = 256
npix = hp.nside2npix(nside)
mesh = np.zeros(npix)
tmesh = np.zeros(npix)

c = sc(ra=RA*u.degree, dec=Dec*u.degree, frame='fk5')

#x = np.remainder(c.galactic.l.degree+360,360)
#ind = x>180
#x[ind] -=360
#x = -x
x = c.galactic.l.degree
y = c.galactic.b.degree

for i in range(len(Dec)):
    theta=np.radians(90.-y[i])
#  phi = y[i]*np.pi/180.
#  theta = x[i]
    phi = np.radians(x[i])
    pix = hp.ang2pix(nside, theta, phi)
    vector = hp.pix2vec(nside,pix)
    pix3 = hp.query_disc(nside,vector,3.*np.pi/180.,inclusive=True)
    pix1 = hp.query_disc(nside,vector,np.pi/180.,inclusive=True)
    tmesh[pix3] = exp[i]/1000.
    tmesh[pix1] = 0.

    mesh += tmesh
    tmesh *= 0.

covered = float(np.count_nonzero(mesh))/float(len(mesh))*100
print('Fraction of sky covered: '+str(covered))

phi = np.linspace(0., 2*np.pi,4000)
theta = np.ones(len(phi))*(90.-29.)*np.pi/180.
pixp = hp.ang2pix(nside,theta,phi)
theta = np.ones(len(phi))*(90.+29.)*np.pi/180.
pixn = hp.ang2pix(nside,theta,phi)
showmap = mesh
showmap[pixp] = max(mesh)
showmap[pixn] = max(mesh)

mask = np.ma.log(showmap)
mapmoll = mask.filled(showmap)

#cmap = cm.bwr
hp.mollview(mapmoll,norm='log', unit='Log(texp/1 ks)',title=' ', cbar=False, cmap='viridis', min=3.0, max=6.0)#, size=2000, cmap='jet') #min=8, max=170,
hp.graticule(dpar=30, dmer=30)



labs = ['150','120','90','60','30','210','240','270','300','330']
angs = [-210,-240,-270,-300,-330,210,240,270,300,330]
for i in range(len(angs)):
    hp.projtext(float(angs[i]),2,labs[i],lonlat=True)
for j in [-60.0, -30.0, 30.0, 60.0]:
    hp.projtext(10,j,str(j),lonlat=True)
#plt.show()

fig = plt.gcf()
ax = plt.gca()
image = ax.get_images()[0]
cmap = fig.colorbar(image, ax=ax, orientation='horizontal')
# Need to provide arguments to create an approriate colorbar
cmap.set_ticks([3.0,4.0,5.0,6.0])
cmap.ax.tick_params(labelsize=15)
for i in range(4):
    cmap.ax.get_xticklabels()[i].set_fontweight(1000)

plt.savefig('Mollweide_blanksky_30deg.png')




'''
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

'''
