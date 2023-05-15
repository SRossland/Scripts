#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntax: MakeDetErrorplot.py workingdir plot questions

# This program combines the making of a file to read all the det values and errors into another program and plot them.
# createDetnorm.py and plotDetnorms.py 

import os, sys, numpy as np, uuid

work_dir = sys.argv[1]
plsa = sys.argv[2]

if plsa.lower() not in ['save','plot']:
  print('Do you want to plot it or save it?')
  sys.exit()

twoplot = False
if len(sys.argv) > 3:
  if sys.argv[3].lower() == 'mcmc':
    mcmc = True
    if not os.path.isfile(work_dir+'/Adetnormerrors_mcmc.txt'):
      current_dir = os.getcwd()
      os.chdir(work_dir)
      os.system('mcmcDetnormprint.py')
      os.chdir(current_dir)
  else:
    twoplot = True
    if not os.path.isfile(work_dir+'/A_olddetnormvals.txt'):
      os.system('det_values_plot.py '+work_dir)


# I will change the directory to inputed directory above and run the programs from there. The final plot will be saved there. This way it allows me to use a generic name for all plots. 


os.chdir(work_dir)
tempfile = os.getcwd()

if not (os.path.isfile(tempfile+'/Adetnormerrors.txt')) and not (os.path.isfile(tempfile+'/Bdetnormerrors.txt')):
  if os.path.isfile(tempfile+'/Adetnormerrors.txt'):
    os.system('rm '+tempfile+'/Adetnormerrors.txt')
  if os.path.isfile(tempfile+'/Bdetnormerrors.txt'):
    os.system('rm '+tempfile+'/Bdetnormerrors.txt')
  os.system('createDetnormfile.py '+tempfile)

if mcmc:
  os.system('plotDetnorms.py '+tempfile+'/Adetnormerrors.txt '+tempfile+'/Bdetnormerrors.txt '+plsa+' mcmc '+tempfile+'/Adetnormerrors_mcmc.txt '+tempfile+'/Bdetnormerrors_mcmc.txt')
if twoplot:
  os.system('plotDetnorms.py '+tempfile+'/Adetnormerrors.txt '+tempfile+'/Bdetnormerrors.txt '+plsa+' '+tempfile+'/A_olddetnormvals.txt '+tempfile+'/B_olddetnormvals.txt')
else: 
  os.system('plotDetnorms.py '+tempfile+'/Adetnormerrors.txt '+tempfile+'/Bdetnormerrors.txt '+plsa) 


