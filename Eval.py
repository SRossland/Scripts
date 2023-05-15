#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./Eval.py list

# This will compare the list of download targets with the current obsids
# in the OBS file currently downloaded (checks if we got them all).
# This will be done by getting the obsid numbers and comparing them to 
# available directories under OBS

import sys, os

listin = sys.argv[1]

dir = os.getcwd()
listout = dir+'/LogFile/'+listin
file = open(listout,'r')
txt = file.readlines()
obs = []

for i in range(len(txt)):
	if len(txt[i]) > 80:
		a = txt[i].split('|')
		obs.append(a)

obsid = []
dater = []

for i in range(len(obs)):
	ob = obs[i]
	obsid.append(ob[5])
	dater.append(ob[11])
	ob = None

date = []
for i in dater:
	a = str(i).replace(' ','')
	date.append(a)

print(date)

# So, Obsid should have a title of OBSID and all id's from the download list
# Lets compare them

k = len(obsid)
dirc = dir+'/OBS'
j = 0

for i in range(1,k):
	if os.path.isdir(dirc+'/'+obsid[i]) == False:
		if len(date[i]) != 0:
			print(obsid[i])
			j += 1

print('End',j)
