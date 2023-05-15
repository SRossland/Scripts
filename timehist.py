import numpy as np
import os, sys
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time

dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs'
file = sys.argv[1]
file2 = sys.argv[2]
dirx = '/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/xray/NuSTAR/OBS'

expA = []
expB = []
dateA = []
dateB = []
obsid = []

f = open(dir+'/'+file,'r+')
for line in f.readlines():
	obsid.append(line[0:11])
f.close()
f = open(dir+'/'+file2,'r+')
for line in f.readlines():
	obsid.append(line[0:11])
f.close()
for i in range(len(obsid)):
	print(obsid[i])
	got_date = 'no'
	eventA = dirx+'/'+obsid[i]+'/event_cl/nu'+obsid[i]+'A01_cl.evt'
	eventB = dirx+'/'+obsid[i]+'/event_cl/nu'+obsid[i]+'B01_cl.evt'

	if os.path.isfile(eventA) == True:
		with fits.open(eventA) as hdulA:
			expA.append(hdulA[1].header['EXPOSURE'])
			da = hdulA[1].header['DATE-OBS']
			date = Time(da, format='isot', scale='utc')
			dateA.append(date.decimalyear)
	if os.path.isfile(eventB) == True:
		with fits.open(eventB) as hdulB:
			expB.append(hdulB[1].header['EXPOSURE'])
			da = hdulB[1].header['DATE-OBS']
			date = Time(da, format='isot', scale='utc')
			dateB.append(date.decimalyear)


expA = np.asarray(expA,dtype=float)
expB = np.asarray(expB,dtype=float)
dateA = np.asarray(dateA,dtype=float)
dateB = np.asarray(dateB,dtype=float)

thic=180
print(dateA[:10])
print(max(expA)-min(expA))

binsA = (max(dateA)-min(dateA))/thic
binsB = (max(dateB)-min(dateB))/thic
print(binsA,binsB)

A_d = zip(dateA,expA)
B_d = zip(dateB,expB)

Acounts = np.zeros(int(binsA)+1)
Bcounts = np.zeros(int(binsB)+1)

A_d.sort()
B_d.sort()

date_min = min(dateA)
date_max = date_min+thic
i = 0
for j in range(len(A_d)):
	if A_d[j][0] > date_max:
		date_min = date_max
		date_max = date_min+thic
		i += 1
	Acounts[i] += A_d[j][1]

date_min = min(dateB)
date_max = date_min+thic
i = 0
for j in range(len(B_d)):
        if B_d[j][0] > date_max:
                date_min = date_max
                date_max = date_min+thic
                i += 1
        Bcounts[i] += B_d[j][1] 


#date_tick = Time(dateA, format='mjd')
#histA = np.histogram(Acounts,bins=int(binsA)+1,density=False)
#histB = np.histogram(Bcounts,bins=int(binsB)+1,density=False)

plt.hist(dateA,weights=expA,alpha=0.5,label='FPMA')
plt.hist(dateB,weights=expB,alpha=0.5,label='FPMB')
#plt.xticks(date_tick.decimalyear)
plt.title('Expsosure Distribution by Date')
plt.xlabel('Date')
plt.ylabel('Exposure (sec)')
plt.show()
