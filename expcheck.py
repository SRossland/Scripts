#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

import numpy as np
import os, sys, tempfile, shutil
from astropy.io import fits

obslist = sys.argv[1]
mode = sys.argv[2]

obs_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS'
expmapdir = obs_dir.replace('OBS','CXB/01/expmap')
logdir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs/'

obs = np.genfromtxt(logdir+obslist,'|U11')

# get exp time for full sun nosun and display it.

scratch = '/scratch/local/u1019304'
os.system('mkdir '+scratch)
os.chdir(scratch)

Path = tempfile.mkdtemp(dir = scratch)
os.chdir(Path)

print(os.getcwd())

if mode == 'NOOCC':
  se = '01'
else:
  se = '02'

for i, ob in enumerate(obs):
  print(ob)
  sep_dir_temp = os.getcwd()+'/'+ob+'/event_sep_cl/'+mode+'/'
  def_dir = obs_dir+'/'+ob+'/event_defcl/'
  sep_dir = obs_dir+'/'+ob+'/event_sep_cl/'
  os.system('cp -r '+obs_dir+'/'+ob+' '+Path)
  for det in ['A','B']:
    print(det)
    os.system('/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/Scripts/rerun_screen_pipe.sh '+Path+'/'+ob+' '+det+' '+mode)
    os.system('mv '+sep_dir+'nu'+ob+det+'**NOSUN.fits '+def_dir)
    os.system('mv '+sep_dir_temp+'NEW/nu'+ob+det+se+'cl.evt '+sep_dir+'nu'+ob+det+se+'_fullevts_NOSUN.fits')
    os.system('mv '+sep_dir_temp+'NEW/nu'+ob+det+se+'_gti.fits '+sep_dir+'nu'+ob+det+se+'01_gti_NOSUN.fits')
    #os.system('mv '+sep_dir_temp+'nu'+ob+det+'**NOSUN.fits '+sep_dir_temp+'/NEW')
    #os.system('mv '+sep_dir_temp+'NEW/nu'+ob+det+se+'_cl.evt '+sep_dir_temp+'TESTFILEMF.fits')
    #os.system('mv '+sep_dir_temp+'NEW/nu'+ob+det+se+'_gti.fits '+sep_dir_temp+'TESTFILEGTI.fits')  
#shutil.remtree(Path)

