#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntax: mcmcDetnormprint.py 

import os, sys, numpy as np

with open('A_params_mcmc.txt','r+') as f:
  lines = f.readlines()

with open('A_params.txt','r+') as f:
  lines_A = f.readlines()
with open('B_params.txt','r+') as f:
  lines_B = f.readlines()


expA = lines_A[0].split()[-1]
expB = lines_B[0].split()[-1]

os.system('mv A_params.txt A_params_inital.txt')
os.system('mv B_params.txt B_params_inital.txt')

# write the new A/B_params.txt files

for line in lines:
  with open('A_params.txt','a+') as O:
    O.write(line.rstrip()+' '+expA+'\n')

with open('B_params_mcmc.txt','r+') as f:
  lines = f.readlines()

for line in lines:
  with open('B_params.txt','a+') as O:
    O.write(line.rstrip()+' '+expB+'\n')

# Do the dets now:

Adets = ['A_det0_params_mcmc.txt','A_det1_params_mcmc.txt','A_det2_params_mcmc.txt','A_det3_params_mcmc.txt']
Bdets = ['B_det0_params_mcmc.txt','B_det1_params_mcmc.txt','B_det2_params_mcmc.txt','B_det3_params_mcmc.txt']

with open('A_det0_params_mcmc.txt','r+') as f:
  linesA0 = f.readlines()
with open('A_det1_params_mcmc.txt','r+') as f:
  linesA1 = f.readlines()
with open('A_det2_params_mcmc.txt','r+') as f:
  linesA2 = f.readlines()
with open('A_det3_params_mcmc.txt','r+') as f:
  linesA3 = f.readlines()

with open('B_det0_params_mcmc.txt','r+') as f:
  linesB0 = f.readlines()
with open('B_det1_params_mcmc.txt','r+') as f:
  linesB1 = f.readlines()
with open('B_det2_params_mcmc.txt','r+') as f:
  linesB2 = f.readlines()
with open('B_det3_params_mcmc.txt','r+') as f:
  linesB3 = f.readlines()

fidetA, fidetB = 'Adetnormerrors_mcmc.txt','Bdetnormerrors_mcmc.txt'

for i,li in enumerate(linesA0):
  crap,crap2, a,b,c = linesA1[i].split()
  crap3,crap4, d,e,f = linesA2[i].split()
  crap5,crap6, g,h,j = linesA3[i].split()
  with open(fidetA,'a+') as O:
    O.write(li.rstrip()+' '+' '.join([a,b,c,d,e,f,g,h,j,expA])+'\n')


for i,li in enumerate(linesB0):
  crap,crap2, a,b,c = linesB1[i].split()
  crap3,crap4, d,e,f = linesB2[i].split()
  crap5,crap6, g,h,j = linesB3[i].split()
  with open(fidetB,'a+') as O:
    O.write(li.rstrip()+' '+' '.join([a,b,c,d,e,f,g,h,j,expA])+'\n')

print('fin')
