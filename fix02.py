import os, sys, numpy as np
obs = np.genfromtxt('/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs/01_OBSLIST.txt',dtype='|U11')

obs_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/'

for i,ob in enumerate(obs):
  for det in ['A','B']:
    olde = obs_dir+ob+'/event_defcl/nu'+ob+det+'01_fullevts_NOSUN.fits'
    oldg = obs_dir+ob+'/event_defcl/nu'+ob+det+'0101_gti_NOSUN.fits'
    newe = obs_dir+ob+'/event_sep_cl/NOOCC/nu'+ob+det+'01_fullevts_NOSUN.fits'
    newg = obs_dir+ob+'/event_sep_cl/NOOCC/nu'+ob+det+'0101_gti_NOSUN.fits'
    defcl = obs_dir+ob+'/event_defcl'
    sepcl = obs_dir+ob+'/event_sep_cl/NOOCC'
    if os.path.isfile(newe):
	if os.path.isfile(olde):
	  archive = sepcl+'/archive'
	  if os.path.isdir(archive) == False: 
		os.system('mkdir '+archive)
	  if os.path.isfile(archive+'/'+newe):
		os.system('rm '+newe)
		os.system('mv '+olde+' '+sepcl)
	  else:
  	        os.system('mv '+newe+' '+archive)
	  	os.system('mv '+olde+' '+sepcl)
    if os.path.isfile(newg):
        if os.path.isfile(oldg):
          archive = sepcl+'/archive'
          if os.path.isdir(archive) == False:
                os.system('mkdir '+archive)
	  if os.path.isfile(archive+'/'+newg):
		os.system('rm '+newg)
		os.system('mv '+oldg+' '+sepcl)
	  else:
          	os.system('mv '+newg+' '+archive)
          	os.system('mv '+oldg+' '+sepcl)
          continue
    if os.path.isfile(olde):
	os.system('mv '+olde+' '+sepcl)
    if os.path.isfile(oldg):
	os.system('mv '+oldg+' '+sepcl)
	
#    else: 
#	with open('BAD_02s.txt','a+') as O:
#		O.write(ob+' '+det+'\n')
