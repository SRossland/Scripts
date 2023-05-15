#!/bin/python3

import numpy as np
import os, sys

tempfile = sys.argv[1]

def saveData(lo, hi, d, el, eh, e):
  return sorted(zip(lo,hi,d,el,eh,e), key = lambda x: x[0])

for i in range(4):
  for j,det in enumerate(['A','B']):
    fi = det+'_det'+str(i)+'_params.txt'
    with open(fi, 'r') as f:
      lines = f.readlines()
    lowe, highe, dat, errl, errh, exp = [],[],[],[],[],[]
    for line in lines:
      li = line.split()
      lowe.append(float(li[0]))
      highe.append(float(li[1]))
      dat.append(li[2])
      errl.append(li[3])
      errh.append(li[4])
      exp.append(li[5])
    if j == 0:
      if i == 0:  Adata0 = saveData(lowe, highe, dat, errl, errh, exp)
      if i == 1:  Adata1 = saveData(lowe, highe, dat, errl, errh, exp)
      if i == 2:  Adata2 = saveData(lowe, highe, dat, errl, errh, exp)
      if i == 3:  Adata3 = saveData(lowe, highe, dat, errl, errh, exp)
    if j == 1:
      if i == 0:  Bdata0 = saveData(lowe, highe, dat, errl, errh, exp)
      if i == 1:  Bdata1 = saveData(lowe, highe, dat, errl, errh, exp)
      if i == 2:  Bdata2 = saveData(lowe, highe, dat, errl, errh, exp)
      if i == 3:  Bdata3 = saveData(lowe, highe, dat, errl, errh, exp)



le1, le2 = [], []
for i,li in enumerate(Adata0):
  le1.append(li[0])
for i,li in enumerate(Adata3):
  le2.append(li[0])
print(set(le1).difference(le2))
writefi = 'detnormerrors.txt'
#print(Adat[4][69][4])

Adatafull, Bdatafull = [],[]
for i in range(len(Adata0)):
  li = str(Adata0[i][0])+' '+str(Adata0[i][1])
  for j in [Adata0,Adata1,Adata2,Adata3]:
    li = li+' '+j[i][2]+' '+j[i][3]+' '+j[i][4]
  li = li+' '+Adata0[i][5]
  Adatafull.append(li)

for i in range(len(Bdata0)):
  li = str(Bdata0[i][0])+' '+str(Bdata0[i][1])
  for j in [Bdata0,Bdata1,Bdata2,Bdata3]:
    li = li+' '+j[i][2]+' '+j[i][3]+' '+ j[i][4]
  li = li+' '+Bdata0[i][5]
  Bdatafull.append(li)


for i,li in enumerate(Adatafull):
#  with open(tempfile_A,'a+') as O:
  with open(tempfile+'/A'+writefi,'a') as O:
    O.write(li+'\n')

for i,li in enumerate(Bdatafull):
#  with open(tempfile_B,'a+') as O:
  with open(tempfile+'/B'+writefi,'a') as O:
    O.write(li+'\n')
    
