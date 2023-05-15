#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntax: viewobs.py obslist writelist 

import os, sys, numpy as np

obsfile = sys.argv[1]
sf = sys.argv[2]
if len(sys.argv) == 4:
  start_index = int(sys.argv[3])
else: 
  start_index = 0

log_dir= '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs/'

obs_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/'

obs = np.genfromtxt(log_dir+obsfile,'|U11')

saved_obs = []

def saveList(obsid, savefile):
  with open(savefile, 'a+') as O:
    O.write(obsid+'\n')

for i,ob in enumerate(obs[start_index:]):
  print(i+start_index)
  eventA = obs_dir+ob+'/event_cl/nu'+ob+'A01_cl.evt'
  eventB = obs_dir+ob+'/event_cl/nu'+ob+'B01_cl.evt'
  excl = obs_dir+ob+'/event_defcl/excl.reg'

  os.system('ds9 '+eventA+' -cmap HSV -scale limits 0 3.0 -zoom 0.75 -smooth -regions '+excl+' '+eventB+' -cmap HSV -scale limits 0 4.0 -smooth -regions '+excl)
  quest = raw_input('Save?')
  if quest.lower() in ['yes','yeah', 'sure', '1', 'y']:
      saveList(ob,sf)



