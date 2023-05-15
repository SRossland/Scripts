#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: maskdetimages.py 

# file is the fits file you with to mask off and is assumed to be a det image

import sys, os, numpy as np
from astropy.io import fits

fifi = raw_input("Fits file to mask:")
name = raw_input("new file name [include .fits]:  ")

dirs = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs/bgdstackimages02/'
	
fi = dirs+fifi

if os.path.isfile(fi) == True:
	filename = dirs+fifi
else: 
	print('File does not exist: Exiting')
	sys.exit()
	
hdul = fits.open(filename)  # This is not working, but I do know fits is loading, and the file can be read in interactive mode
data = hdul[0].data
print('A')

print('B')
data[:23,:] = 0
data[339:,:] = 0
data[179:186,:] = 0
data[:,:23] = 0
data[:,342:] = 0
data[:,179:186] = 0


# This space is open to create a mask of high points (use a statistical approach).

print('C')
hdulw = fits.PrimaryHDU(data)
hdulw.writeto(dirs+name)

hdul.close()
