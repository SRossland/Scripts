#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

import numpy as np
import os, sys, subprocess
from astropy.io import fits

obslist = sys.argv[1]
mode = sys.argv[2]
#scratch = sys.argv[3]
de = sys.argv[3]

obs_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS'
logdir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs/'

obs = np.genfromtxt(logdir+obslist,'|U11')
#os.chdir(obs_dir)

#if os.getcwd() != scratch:
#  print('Not the right Directory: '+os.getcwd())
#  print('Should be in: '+scratch)
#  sys.exit()

if mode == 'NOOCC':
  se = '01'
  OCC='SCIENCE'
else:
  se = '02'
  OCC='OCCULTATION'


for i, ob in enumerate(obs):
  print(ob)
  #replace_old(ob, de)
  sep_dir_temp = obs_dir+'/'+ob+'/event_sep_cl/'+mode+'/'
  def_dir = obs_dir+'/'+ob+'/event_defcl/'
  sep_dir = obs_dir+'/'+ob+'/event_sep_cl/'
  #os.system('cp -r '+obs_dir+'/'+ob+' '+scratch)
  #for det in ['B']:
  #for det in ['A','B']:
  gti_file = sep_dir_temp+'/nu'+ob+de+se+'01_gti_NOSUN.fits'
  evt_file = gti_file.replace('01_gti_','_fullevts_')
  if os.path.isfile(gti_file):
      os.system('cp -a '+gti_file+' '+def_dir)
      os.system('cp -a '+evt_file+' '+def_dir)
   #I should be copying the gti file and changing that one to read into cutting
      with fits.open(gti_file) as hdul:
	  start = hdul[1].data['START']
	  stop = hdul[1].data['STOP']
      
      start += 300
      stop -= 300
      idx = np.where((stop-start) < 0)
      start[idx] = None
      stop[idx] = None

      with fits.open(gti_file, mode='update') as hdul:
  	  #hdul[1].data['START'] += 300
          #hdul[1].data['STOP'] -= 300
	  hdul[1].data['START'] = start
	  hdul[1].data['STOP'] = stop
	  hdul.flush()
      
      # fix the neg values in the gti's  # for 02 gtiexpr='ELV<=5'
      os.system('nuscreen infile=./'+ob+'/event_sep_cl/'+mode+'/nu'+ob+de+se+'_fullevts_NOSUN.fits hkfile=./'+ob+'/event_cl/nu'+ob+de+'_fpm.hk mkffile=./'+ob+'/event_cl/nu'+ob+de+'.mkf outdir=./'+ob+'/event_sep_cl/NOOCC/TEST gtiscreen=yes usrgtifile=./'+ob+'/event_sep_cl/NOOCC/nu'+ob+de+se+'01_gti_NOSUN.fits evtscreen=no gtiexpr=SUNSHINE==0 outfile=TEST'+se+de+'.evt clobber=yes')
     
      #os.system('nuscreen obsmode='+OCC+' infile=./'+ob+'/event_sep_cl/'+mode+'/nu'+ob+de+se+'_fullevts_NOSUN.fits gtiscreen=yes evtscreen=yes gtiexpr=DEFAULT gradeexpr=DEFAULT statusexpr=DEFAULT usrgtifile=./'+ob+'/event_sep_cl/'+mode+'/nu'+ob+de+se+'01_gti_NOSUN.fits createattgti=yes createinstrgti=yes outdir=./'+ob+'/event_sep_cl/'+mode+'/TEST hkfile=./'+ob+'/event_cl/nu'+ob+de+'_fpm.hk mkffile=./'+ob+'/event_cl/nu'+ob+de+'.mkf outfile=TEST'+se+de+'.evt')

#      os.system('nuscreen infile=./'+ob+'/event_sep_cl/'+mode+'/nu'+ob+de+se+'_fullevts_NOSUN.fits gtiscreen=yes evtscreen=no gtiexpr=NONE gradeexpr="DEFAULT" statusexpr="DEFAULT" outdir=./'+ob+'/event_sep_cl/'+mode+'/TEST obsmode='+OCC+' usrgtifile=./'+ob+'/event_sep_cl/'+mode+'/nu'+ob+de+se+'01_gti_NOSUN.fits hkfile=./'+ob+'/event_cl/nu'+ob+de+'_fpm.hk mkffile=./'+ob+'/event_cl/nu'+ob+de+'.mkf outfile=TEST'+se+de+'.evt clobber=yes')
      os.system('cp -a '+sep_dir_temp+'TEST/TEST'+se+de+'.evt '+sep_dir_temp+'/nu'+ob+de+se+'_fullevts_NOSUN.fits')
      os.system('rm -rf '+sep_dir_temp+'TEST')
      #os.system('/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/Scripts/rerun_screen_pipe.sh '+obs_dir+'/'+ob+' '+de+' '+mode)
      #os.system('mv '+sep_dir+mode+'/nu'+ob+de+se+'_fullevts_NOSUN.fits '+def_dir)
      #os.system('mv '+sep_dir+mode+'/nu'+ob+de+se+'01_gti_NOSUN.fits '+def_dir)
      #os.system('cp '+sep_dir_temp+'NEW/nu'+ob+de+se+'_cl.evt '+sep_dir+mode+'/nu'+ob+de+se+'_fullevts_NOSUN.fits')
      #os.system('cp '+sep_dir_temp+'NEW/nu'+ob+de+se+'_gti.fits '+sep_dir+mode+'/nu'+ob+de+se+'01_gti_NOSUN.fits')
  else:
      continue
