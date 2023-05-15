import numpy as np, os, sys
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord as sc

dirs = '/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/xray/NuSTAR/OBS'
logsdir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs'

def check(obsid, dirs):
    eventA = dirs+'/'+obsid+'/event_cl/nu'+obsid+'A01_cl.evt' ; eventB = dirs+'/'+obsid+'/event_cl/nu'+obsid+'B01_cl.evt'
    if os.path.isfile(eventA) == True:
      with fits.open(eventA) as hdul:
        exA = hdul[1].header['EXPOSURE']
        RA_A = hdul[0].header['RA_PNT']
        DEC_A = hdul[0].header['DEC_PNT']
    if os.path.isfile(eventB) == True:
      with fits.open(eventB) as hdul:
        exB = hdul[1].header['EXPOSURE']
        RA_B = hdul[0].header['RA_PNT']
        DEC_B = hdul[0].header['DEC_PNT']
    if ('RA_A' in locals()) == True:
      if 'RA_A' != 0.:
        RA = RA_A
        DEC = DEC_A
    elif ('RA_B' in locals()) == True:
      if 'RA_A' != 0.:
        RA = RA_B
        DEC = DEC_B
    else:
      gr = 0; exA = 0.; exB=0.; RA=0.; DEC = 0.
      return gr, exA, exB, RA, DEC
    c = sc(ra = RA*u.degree, dec = DEC*u.degree, frame='fk5')
    if abs(c.galactic.b.deg) >= 10.0:
      gr = 1
    else:
      gr = 0
    if ('exA' in locals()) == False:
      exA = 0.
    if ('exB' in locals()) == False:
      exB = 0.
    return gr, exA, exB, RA, DEC

f = open(logsdir+'/Zero_reg_May1.txt','r+')
g = open(logsdir+'/Src_reg_May1.txt','r+')
h = open(logsdir+'/Large_src_May1.txt','r+')
obsids = []
zero, srcobs, Largesrc = [],[],[]
for line in f.readlines():
  obsids.append(line.rstrip())
  zero.append(line.rstrip())
for line in g.readlines():
  obsids.append(line.rstrip())
  srcobs.append(line.rstrip())
for line in h.readlines():
  obsids.append(line.rstrip())
  Largesrc.append(line.rstrip())
f.close()
g.close()
h.close()
file = 'MIT_obslist.txt'
with open(file, 'w+') as O:
  O.write(('{:^13} | {:^11} | {:^11} | {:^16} | {:^16}').format('OBSID','RA','DEC','Exp. A','Exp. B'))
  O.write('\n')

expA = 0; expB = 0

#j = 0

#for i in range(len(obsids)):
#  if obsids[i][0] == '7':
#    continue
#  gr, exA, exB, ra, dec = check(obsids[i],dirs)
#  if gr == 1:
#    expA += exA; expB += exB; j += 1
#    with open(file, 'a+') as O:
#      O.write(('{:^13} | {:^11} | {:^11} | {:^16} | {:^16}').format(obsids[i],ra,dec,exA,exB))
#      O.write('\n')

breakdown = [zero, srcobs, Largesrc]
names = ['Zero_obs','src_obs','Largesrc_obs']
for i in range(3):
  obs = breakdown[i]
  k = 0; eA = 0; eB = 0
  for j in range(len(obs)):
    if obs[j][0] == '7':
       continue
    gr, exA, exB, blah, bla = check(obs[j],dirs)
    if gr == 1:
      with open(str(names[i])+'.txt','a+') as O:
        O.write(obs[j]+'\n')
      k += 1; eA += exA; eB += exB
  print(k, eA, eB)



#print(j)
#print(len(obsids))
#with open(file, 'a+') as O:
#  O.write('='*11+'\n')
#  O.write('Exposure for A: '+str(expA))
#  O.write('\n'+'Exposure for B: '+str(expB))
