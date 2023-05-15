#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# syntax: thisprogram directory

# This program reads the directory for all longformat files from the CXB directory and pulls the norms for the detectors

import os, sys, numpy as np
import matplotlib.pyplot as plt

txtlist = sys.argv[1]

rawlist = os.listdir(txtlist)
longformatlist = []

with open(txtlist+'/A_params.txt') as f:
  first_line = f.readline().rstrip()
expA = first_line.split()[5]

with open(txtlist+'/B_params.txt') as f:
  first_line = f.readline().rstrip()
expB = first_line.split()[5]

for i,st in enumerate(rawlist):
  if 'longformat' in st:
    longformatlist.append(st)

Alist, Blist = [],[]

for i,st in enumerate(longformatlist):
  if 'A_' in st:
    Alist.append(st)
  if 'B_' in st:
    Blist.append(st)

det0A,det1A,det2A,det3A,det0B,det1B,det2B,det3B = [],[],[],[],[],[],[],[]


for i,fi in enumerate(Alist):
  f = open(txtlist+'/'+fi)
  for line in f.readlines():
    if 'det0 norm' in line:
      det0A.append(line.split()[3])
    if 'det1 norm' in line:
      det1A.append(line.split()[3])
    if 'det2 norm' in line:
      det2A.append(line.split()[3])
    if 'det3 norm' in line:
      det3A.append(line.split()[3])
  f.close()

for i,fi in enumerate(Blist):
  f = open(txtlist+'/'+fi)
  for line in f.readlines():
    if 'det0 norm' in line:
      det0B.append(line.split()[3])
    if 'det1 norm' in line:
      det1B.append(line.split()[3])
    if 'det2 norm' in line:
      det2B.append(line.split()[3])
    if 'det3 norm' in line:
      det3B.append(line.split()[3])
  f.close()


# Now we should have all the values of the dets for both A and B

# Now we need to create the X values, which will be the minimum value of the energy bin.

Ax, Bx = [],[]

for i,fi in enumerate(Alist):
  Ax.append(fi.split('_')[2])

for i,fi in enumerate(Blist):
  Bx.append(fi.split('_')[2])

for i in range(len(det0A)):

	with open(txtlist+'/A_olddetnormvals.txt','a+') as O:
		O.write(det0A[i]+' '+det1A[i]+' '+det2A[i]+' '+det3A[i]+' '+Ax[i]+' '+expA+'\n')
for i in range(len(det0B)):
	with open(txtlist+'/B_olddetnormvals.txt','a+') as O:
                O.write(det0B[i]+' '+det1B[i]+' '+det2B[i]+' '+det3B[i]+' '+Bx[i]+' '+expB+'\n')

# plotting
'''
plt.scatter(Ax,det0A,c='blue',label='det0')
plt.scatter(Ax,det1A,c='red',label='det1')
plt.scatter(Ax,det2A,c='green',label='det2')
plt.scatter(Ax,det3A,c='black',label='det3')
plt.title('A det norms for energies '+str(Ax[0])+' to '+str(Ax[-1]))
plt.xlabel('Energy keV')
plt.ylabel('Norm Value')
plt.legend()
plt.show()

plt.scatter(Bx,det0B,c='blue',label='det0')
plt.scatter(Bx,det1B,c='red',label='det1')
plt.scatter(Bx,det2B,c='green',label='det2')
plt.scatter(Bx,det3B,c='black',label='det3')
plt.title('B det norms for energies '+str(Bx[0])+' to '+str(Bx[-1]))
plt.xlabel('Energy keV')
plt.ylabel('Norm Value')
plt.legend()
plt.show()
'''
