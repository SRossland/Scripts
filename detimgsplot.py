#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntas: detimgsplot.py path/to/obslist/obslist low_energy high_energy path/to/write/directory

import os, sys, pidly, numpy as np, matplotlib.pyplot as plt
from astropy.io import fits


obsid = np.genfromtxt(sys.argv[1],'|U11')
elow = sys.argv[2]
ehigh = sys.argv[3]
curr_dir = sys.argv[4]
if len(sys.argv) > 5:
  idx = int(sys.argv[5])
else: 
  idx = 0
obs_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/'

idl = pidly.IDL('/uufs/chpc.utah.edu/sys/pkg/idl/8.4/idl84/bin/idl')

os.chdir(obs_dir)


j = 0

for i,ob in enumerate(obsid[idx:]):
  cldir = obs_dir+'OBS/'+ob+'/event_cl'
  eventA = cldir+'/nu'+ob+'A01_cl.evt'
  eventB = cldir+'/nu'+ob+'B01_cl.evt'
  idl("cd, current = dir")
  idl("dir = dir+'OBS/'")
  idl("obsid='"+ob+"'")
  idl("cldir = '"+cldir+"/'")
  if os.path.isfile(eventA):
    idl("mkimgs,cldir,obsid,'A',"+elow+","+ehigh)
    os.system('mv '+cldir+'/imA'+elow+'to'+ehigh+'keV.fits '+curr_dir)
    j+=1
  if os.path.isfile(eventB):
    idl("mkimgs,cldir,obsid,'B',"+elow+","+ehigh)
    os.system('mv '+cldir+'/imB'+elow+'to'+ehigh+'keV.fits '+curr_dir)
    j+=2

  if j == 3:
    os.system('ds9 '+curr_dir+'/*A*keV.fits -cmap HSV -scale limits 0 0.5 -zoom 0.75 -smooth '+curr_dir+'/*B*keV.fits -cmap HSV -scale limits 0 0.5 -smooth')
  if j == 1:
    os.system('ds9 '+curr_dir+'/*A*keV.fits -cmap HSV -scale limits 0 0.5 -zoom 0.75 -smooth')
  if j == 2:
    os.system('ds9 '+curr_dir+'/*B*keV.fits -cmap HSV -scale limits 0 0.5 -zoom 0.75 -smooth') 
 
  choi = raw_input('save?')
  if choi.lower() in ['y','yes','1','01','yeah','sure']:
    with open(curr_dir+'/blankdata_obsid_with_issues.txt','a+') as O:
      O.write(ob+'\n')
  os.system('rm '+curr_dir+'/im*'+elow+'to'+ehigh+'keV.fits')
  j = 0

