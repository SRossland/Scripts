#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# Syntax: Flux_results.py file write_file return_length

# Purpose: This program takes the file given assuming it is a readout from ExcessFlux.py and 
#          returns a list of the highest flux observations for a specified amountto a given 
#          file. This will not discriminate between the adjusted and non-adjusted files. 

# Modules:
import numpy as np
import os, sys
import matplotlib.pyplot as plt

fi = sys.argv[1]
wfi = sys.argv[2]
ret_len = int(sys.argv[3])*-1

det = fi.split('_')[4]

obsid, flux = [],[]

with open(fi,'r+') as f:
  lines = f.readlines()

for line in lines:
  obsid.append(line.split()[0])
  flux.append(line.split()[1])

flux_float = [float(i) for i in flux]

sorted_list = sorted(zip(obsid,flux_float), key = lambda x: x[1])
obs,flx = zip(*sorted_list)

ret_flux = list(flx[ret_len:len(flx)])
ret_obs = list(obs[ret_len:len(obs)])

ret_flux.reverse()
ret_obs.reverse()

for i,blah in enumerate(ret_obs):
  with open(det+'_'+wfi,'a+') as O:
    O.write(blah+' '+str(ret_flux[i])+'\n')




