#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntax: newaCXBvalues.py /path/to/directory  # This directory must have the long format txt files and the [A/B]_params.txt files.

import os, sys, numpy as np

direct = sys.argv[1]
if direct.lower() == 'current':
	direct = os.getcwd()

# Check for the needed files:

# Workflow: 	(1) check for files
#		(2) get energy limits from file names
#		(3) load new aCXB values
#		(4) check for gaps in energy values of the new aCXB values, if gaps occur, pull old values
#		(5) compile list in params file
#		(6) rename old file
#		(7) write new file

if (os.path.isfile('A_params.txt') and os.path.isfile('B_params.txt')):
  dets = ['A','B']
if (os.path.isfile('A_params.txt') and not os.path.isfile('B_params.txt')):
  dets = ['A']
if (os.path.isfile('B_params.txt') and not os.path.isfile('A_params.txt')):
  dets = ['B']
if not (os.path.isfile('A_params.txt') and os.path.isfile('B_params.txt')):
  print('No param files are found! EXITING!')
  sys.exit()

li = os.listdir(direct)
A_list, B_list, A_enerlow, B_enerlow, A_enerhigh, B_enerhigh = [],[],[],[],[],[]

for i,blah in enumerate(li):
  if 'longformat' in blah:
    if 'A_nosun' in blah:
      A_list.append(blah)
      A_enerlow.append(blah.split('_')[2])
      A_enerhigh.append(blah.split('_')[3])
    else:
      B_list.append(blah)
      B_enerlow.append(blah.split('_')[2])
      B_enerhigh.append(blah.split('_')[3])


if len(A_list) != len(B_list):
  print('Missing energy files: ')
  print('A energies:' )
  print(A_enerlow)
  print('B energies:' )
  print(B_enerlow)
  print('Missing energies are:')
  print((set(A_enerlow).difference(B_enerlow)))
  sys.exit()

lowe = A_enerlow
highe = A_enerhigh


# we need to sort the list of energies and files:

A_crap = sorted(zip(A_list, lowe, highe), key = lambda x: float(x[1]))
B_crap = sorted(zip(B_list, B_enerlow, B_enerhigh), key = lambda x: float(x[1]))

A_list, lowe, highe = zip(*A_crap)
B_list, B_lowcrap, B_highcrap = zip(*B_crap)

# Now we need to load the new energies
Afinal, Bfinal = [],[]
A_newvals, B_newvals = [],[]
for det in dets:
  with open(direct+'/'+det+'_params.txt') as f:
    if det == 'A': oldlines_A = f.readlines()
    else: oldlines_B = f.readlines()
  with open(direct+'/'+det+'_aCXBtemp_params.txt') as f:
    lines = f.readlines()
  lowtest = '0.0'
  for line in lines:

    if line == lines[-1]:
      if det == 'A':  A_newvals.append(line)
      else: B_newvals.append(line)
      continue

    if (line.split()[0] == lowtest.split()[0]) or (lowtest == '0.0'):
      lowtest = line
      continue
    else:
      if det == 'A':
        A_newvals.append(lowtest)
        lowtest = line
      else:
        B_newvals.append(lowtest)
        lowtest = line
  
# Something is going wrong in the print so I need to fix this part.

# find the gaps in info
gaps_A, gaps_B = [], []
for det in dets:
  if det == 'A':
    gaps_A = [i.split()[0] for i in A_newvals]
    gaps_A = list(sorted(set(lowe).difference(gaps_A)))
  if det == 'B':
    gaps_B = [i.split()[0] for i in B_newvals]
    gaps_B = list(sorted(set(lowe).difference(gaps_B)))


for det in dets:
  if det == 'A':
    j = 0
    for i,en in enumerate(lowe):
      if en in gaps_A:
        Afinal.append(oldlines_A[i])
        j += 1
        continue
      Afinal.append(A_newvals[i-j])
  else:
    j = 0
    for i,en in enumerate(lowe):
      if en in gaps_B:
        Bfinal.append(oldlines_B[i])
        j+=1
        continue
      Bfinal.append(B_newvals[i-j])

# Afinal and Bfinal should be my new values in order.

for det in dets:
  os.system('mv '+direct+'/'+det+'_params.txt '+direct+'/'+det+'_params_firstrun.txt')
  if os.path.isfile(direct+'/spec'+det+'.pha'):
    os.system('mv '+direct+'/spec'+det+'.pha '+direct+'/spec'+det+'_firstrun.pha')
  if det == 'A':
    for line in Afinal:
      with open(direct+'/'+det+'_params.txt','a+') as O:
        O.write(line.rstrip()+'\n')
  elif det == 'B':
    for line in Bfinal:
      with open(direct+'/'+det+'_params.txt','a+') as O:
        O.write(line.rstrip()+'\n')

print('Finished')
     







 
