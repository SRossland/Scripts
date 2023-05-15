#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./sourcecheck.py

# This is a script to check if the existence an excl file is present, if so
# then it will ignore the obs, if it doesn't have one, then those will be 
# passed to a list to be used in the stack for summer of 2018
# A second part of the code is going to tell us how much exp time is available
# This exp time will be split into 2 sections, OCCULTATED and NONOCCULTATED

import os, sys
from astropy.io import fits

dir = '/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/xray/NuSTAR/'
exclless = []
excl = []

########## Get list of good obs ##################


file = dir+'LogFile/goodobs_nodupes.dat'

gobslist = []

f = open(file, 'r+')
for line in f.readlines():
	gobslist.append(line.rstrip())
f.close()

exptotA01 = 0 
exptotB01 = 0
exptotA02 = 0
exptotB02 = 0
######### Check for source excl file ############

with open('../logs/OBSDATA1.txt','w+') as o:
	o.write('{:^11} | {:^15} | {:^15} | {:^15} | {:^15} | {:^5} | {:^8}'.format('OBSID','SRC Exp A','SRC Exp B','OCC Exp A','OCC Exp B','SRC','RAD'))
	o.write('\n')

try :
	for i in range(len(gobslist)):
		print(i)	
		cldir = dir+'OBS/'+gobslist[i]+'/event_cl'
		### get exp times ###
		### On source time ###
		srcfileA = cldir+'/nu'+gobslist[i]+'A01_cl.evt'
		srcfileB = cldir+'/nu'+gobslist[i]+'B01_cl.evt'
		if os.path.isfile(srcfileA) == True:
			imgarA = fits.open(srcfileA)
			hdrA = imgarA[0].header
			expA = hdrA['Exposure']
			exptotA01 = exptotA01 + float(expA)
			imgarA.close()
		else:
			expA = '0'
		
		if os.path.isfile(srcfileB) == True:
			imgarB = fits.open(srcfileB)
			hdrB = imgarB[0].header
			expB = hdrB['Exposure']
			exptotB01 = exptotB01 + float(expB)
			imgarB.close()
		else:
			expB = '0'

		path = dir+'OBS/'+gobslist[i]+'/event_defcl'
		pathA = path+'/nu'+gobslist[i]+'A02_cl.evt'
		pathB = path+'/nu'+gobslist[i]+'B02_cl.evt'
		if os.path.isfile(pathA) == True:
			imgarA = fits.open(pathA)
			hdrA = imgarA[0].header
			expA02 = hdrA['Exposure']
			exptotA02 = exptotA02 + float(expA02)
			imgarA.close()
		else:
			expA02 = '0'

		if os.path.isfile(pathB) == True:
			imgarB = fits.open(pathB)
			hdrB = imgarB[0].header
			expB02 = hdrB['Exposure']
			exptotB02 = exptotB02 + float(expB02)		
		else:
			expB02 = '0'

		if os.path.isfile(path+'/excl.reg') == True:
			excl_file = path+'/excl.reg'
			fe = open(excl_file,'r+')
			crap = []
			for line in fe.readlines():
				l = line.rstrip()
				crap.append(l)
			a,b,rad = crap[3].split(",")
			rad = rad.rstrip(")")
			excl.append(gobslist[i])
			with open('../logs/OBSDATA1.txt','a+') as o:
				o.write('{:11} | {:15} | {:15} | {:15} | {:15} | {:^5} | {:8}'.format(gobslist[i],expA,expB,expA02,expB02,'yes',rad))
				o.write('\n')
		else:
			exclless.append(gobslist[i])
			with open('../logs/OBSDATA1.txt','a+') as o:
                                o.write('{:11} | {:15} | {:15} | {:15} | {:15} | {:^5} | {:8}'.format(gobslist[i],expA,expB,expA02,expB02,'no','0'))
				o.write('\n')

	with open('../logs/OBSDATA1.txt','a+') as o:
		o.write('Total Stats:'+'\n')
		o.write('    Total OBS: '+str(len(gobslist))+'\n')
		o.write('    OBS with SRC:  '+str(len(excl))+'\n')
		o.write('    OBS w/o SRC:   '+str(len(exclless))+'\n')
		o.write('    Total SRC time A:  '+str(exptotA01)+'\n')
		o.write('    Total OCC time A:  '+str(exptotA02)+'\n')
		o.write('    Total SRC time B:  '+str(exptotB01)+'\n')
		o.write('    Total OCC time B:  '+str(exptotB02)) 
	
	for i in range(len(excl)):
		with open('../logs/ExclOBS.dat','a+') as o:
			o.write(excl[i]+'\n')
	
	for i in range(len(exclless)):
		with open('../logs/NOExclOBS.dat','a+') as o:
			o.write(exclless[i]+'\n')

except Exception as e:
	print(e)
	print("Error, check file")
	sys.exit()
		
