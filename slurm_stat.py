#!/bin/python3
'''#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python'''

# Here's the rundown -- This script calls the batch submission process (batch.slurm). Ideally,
# that will then run image_stats.py, a python2 script that first creates a list of energies
# to pass, then passes them to a stacking program (image_stack.py) and a statistics script (count_stat_new.py)
# Once all those scripts have been completed, the slurm announces it is complete, and this script 
# closes out any temp directories created in the process (it is the image_stat.py program that moves all
# fits, images, and txt files to the final directory beforehand). 

# system: slurm_stat.py low_energy high_energy energy_step mode
# Definition: low_energy = the lowest energy to be passes through the process, NuSTAR is limited reliably
#			   to 3.0 keV
#	      high_energy = the highest energy to be passed. Should limit it to no higher than 40
# 	      energy_step = the factor that the low energy is multiplied by to created the intermidiate 
#			    steps, i.e. 3.0*energy_step = new_energy so range is -- 3.0 - new_energy, then 
#			    repeated until new_energy >= high_energy
#	      mode = 01 or 02 for science or obscurred data
#             full = skiping those obs with excl.reg files or not
# Note: The energy is multiplied by the energy_step and then converted to a energy bin to ensure each step
#       is a real energy value in NuSTAR (it has 0.04 keV steps between 1.6 - 165.4, or channel number from 
#       0-4095. the formulation is Energy = channel_number * 0.04 + 1.6)

import tempfile, os, sys
#import tempfile, shutil, os, sys

archive = sys.argv[1]
datadir = sys.argv[2]
mode = sys.argv[3]
le = sys.argv[4]
he = sys.argv[5]
energy_step = sys.argv[6]
skip = sys.argv[7]

first_dir = os.getcwd()
Path = tempfile.TemporaryDirectory(dir = first_dir)
os.chdir(Path.name)
os.system('image_stats.py '+archive+' '+datadir+' '+mode+' nosun '+le+' '+he+' BOTH '+os.getcwd()+' '+energy_step+' '+skip)
#os.system('batch.slurm '+os.getcwd()+' '+le+' '+he+' '+energy_step+' '+mode+' '+skip)
os.chdir(first_dir)
Path.cleanup()

# Could think about using a 'with' command ex:
# with tempfile.TemporaryDirectory(dir = './') as path:
#     os.chdir(path)
#     os.sytem('batch.slurm '+blah)
#     os.chdir('..')
'''
##### Python 2 code (tempfile worked differently in 2)
first_dir = os.getcwd()
dirpath = tempfile.mkdtemp()
os.chdir(dirpath)
os.system('batch.slurm /uufs/astro.utah.edu/common/home/u1019304/'+path_current+' '+le+' '+he+' '+energy_step)
os.chdir(first_dir)
shutil.rmtree(dirpath) '''

