#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# syntax: ./nodupes.py

# Purpose:  quickly checks for duplicates in the obsid list for both good and bad obs

# Created by: Steven P. Rossland 2018

import numpy as np, os, sys

dir = '/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/xray/NuSTAR/LogFile'

if os.path.isfile(dir+'/goodobs_nodupes.dat') == True:
	wlist = raw_input('good dupe list exists, give a file name to rename list:  ')
	mlist, crap = wlist.split(".")
	newfile = mlist+'.dat' 
	if os.path.isfile(dir+'/'+newfile) == True:
		print("That file also exists...wtf are you doing?")
		sys.exit()
	os.system('mv '+dir+'/goodobs_nodupes.dat '+dir+'/'+newfile)

if os.path.isfile(dir+'/badobs_nodupes.dat') == True:
        wlist = raw_input('bad dupe list exists, give a file name to rename list:  ')
        mlist, crap = wlist.split(".")
        newfile = mlist+'.dat'
        if os.path.isfile(dir+'/'+newfile) == True:
                print("That file also exists...wtf are you doing?")
                sys.exit()
        os.system('mv '+dir+'/badobs_nodupes.dat '+dir+'/'+newfile)


oldobs = []
glah = dir+'/goodobs.dat'
f = open(glah,'r+')
for line in f.readlines():
	if len(line) == 0:
		break
	ob = line[0:11]
	oldobs.append(ob)
f.close()

bad_oldobs = []
blah = dir+'/badobs.dat'
g = open(blah, 'r+')
for line in g.readlines():
	if len(line) == 0:
		break
	ob = line[0:11]
	bad_oldobs.append(ob)
g.close()

if len(oldobs) == len(set(oldobs)):
	os.system('mv '+glah+' '+dir+'/goodobs_nodupes.dat')
else:	
	newgoodobs = list(set(oldobs))
	for i in range(len(newgoodobs)):
		with open(dir+'/goodobs_nodupes.dat','a+') as o:
			o.write(newgoodobs[i]+'\n')
	os.system('rm -f '+glah)

if len(bad_oldobs) == len(set(bad_oldobs)):
	os.system('mv '+blah+' '+dir+'/badobs_nodupes.dat')
else:
	newbadobs = list(set(bad_oldobs))
	for i in range(len(newbadobs)):
		with open(dir+'/badobs_nodupes.dat','a+') as b:
			b.write(newbadobs[i]+'\n')
	os.system('rm -f '+blah)

print("duplicates erradicated")
