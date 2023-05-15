#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# syntax: count_stat.py det

# This purpose of this script is to use the counts in an image from 3-5keV (currently) in a least squares approach to find the values of:
#  Mi,j = a*Ai,j + b1 + b2 + b3 + b4, where M is the model value, A is the physical location of the pixel, a is a variable to find, and b's are the base values of the individual dets.


import os, sys, numpy as np
from astropy.io import fits
import scipy.optimize as opt
#from scipy.optimize import minimize
import matplotlib.pyplot as plt


detp = sys.argv[1]

# Things I need to do for this program:
  # pass the file to be processed in the command line
  # pass if it's 01 or 02 files


##################
# files to use for counts:  pixmapA/B in nuskybgd auxil
# 			    edgemaskA/B in auxil
#			    

# local files
nusky_dir = '/uufs/astro.utah.edu/common/home/u1019304/my_idl/nuskybgd/auxil/'
edgeA = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/fullmaskA.fits'
edgeB = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/fullmaskB.fits'
local_dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/'


#######################

# This part is the image count.  The image is found in the auxil file

# NEED TO CHANGE THIS TO PASSING THE FILE IN THE COMMAND LINE####***********IMPORTANT

#dataA = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/A_3_5SUN_masked.fits'
#dataB = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/B_3_5SUN_masked.fits'
#dataA = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/spec/dataANOSUN_2.fits'
#dataB = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/spec/dataBNOSUN_2.fits'
#dataA = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/spec/DataA_NOSUN_02_3_5keV.fits'
#dataB = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/spec/DataB_NOSUN_02_3_5keV.fits'
dataA = '/uufs/astro.utah.edu/common/home/u1019304/DataA_01_full_3_8keV.fits'
dataB = '/uufs/astro.utah.edu/common/home/u1019304/DataB_01_full_3_8keV.fits'

## Test data
#dataA = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/dataA_flat.fits'
#dataB = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/dataB_flat.fits'
#dataA = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/dataA_flat_detvar.fits'
#dataB = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/dataB_flat_detvar.fits'
#dataA = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/dataA_apgrad.fits'
#dataB = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/dataB_apgrad.fits'
#dataA = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/dataA_apgrad_detvar.fits'
#dataB = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/dataB_apgrad_detvar.fits'

if detp == 'A':
    edge = edgeA
    data = dataA
else:
    edge = edgeB
    data = dataB

# Load data

with fits.open(nusky_dir+'detnum'+detp+'.fits') as hA:
  dA = hA[1].data

with fits.open(edge) as hA:
  eA = hA[0].data

with fits.open(data) as hA:
  datA = hA[0].data

with fits.open(nusky_dir+'det'+detp+'_det1.img') as gA:
  gradA = gA[0].data

#################################

# fix the Aij, Bij arrays 

Aij = gradA[len(gradA)/2 - 180:len(gradA)/2 + 180, len(gradA)/2 - 180:len(gradA)/2 + 180]*eA

# So, for each value where the mask is not 0, we will run the fuction and minimize to that value.

dA += 1

vA = dA*eA


##################
# Counts of Aij for A 

a0 = np.zeros((360,360))
a1 = np.zeros((360,360))
a2 = np.zeros((360,360))
a3 = np.zeros((360,360))


a0[vA == 1]=1
a1[vA == 2]=1
a2[vA == 3]=1
a3[vA == 4]=1

#Bij = np.sum(eB)

##################

PA = np.copy(datA)

#################
# Rebin the data arrays.  have to ensure that each array
# has at least a few counts in each bin

binA = np.zeros((360,360))

gA = 0 # this is the sum of the counts in the Aij, Bij arrays
gA_mod = 0
kA = 1 # this value will be used as a bin value
min_value = 300.0

for i in range(360):
    for j in range(360):
        gA += PA[i,j]
        if (gA >= min_value): #and (gA_mod >= min_value):
            binA[i,j] = kA
            gA = 0
            kA += 1
        else:
            binA[i,j] = kA 
        if (i == 359) and (j == 359) and (gA < min_value):
            As = np.where(binA == kA)
            binA[As] = kA - 1
            kA -= 1
            

#print(kA)

#plt.matshow(binA)
#plt.savefig(local_dir+'bin'+detp+'.png')
#hduAbin = fits.PrimaryHDU(binA)
#hduAbin.writeto(local_dir+'bin'+detp+'.fits')


bins_of_A = np.zeros(kA)
bins_of_modA = np.zeros(kA)

for i in range(kA):
    vals_of_PA = np.where(binA == i+1)
    bins_of_A[i] = np.sum(PA[vals_of_PA])
    if bins_of_A[i] == 0.0:
        print('failed to create bins')
        print(str(i) + ' out of '+str(kA)+' bins')
        sys.exit()

'''for i in range(kB):
    vals_of_PB = np.where(binB == i+1)
    bins_of_B[i] = np.sum(PB[vals_of_PB])'''
    

#################
def detfn(param,detnum):
    global PA
    moAdet = param*detnum
    idx = np.where((PA*detnum) != 0)
    #idx = np.where(detnum == 1)
    #bins_of_det = np.zeros(int(np.sum(detnum)))
    #bins_of_moddet = np.zeros(int(np.sum(detnum)))
    #for i in range(len(bins_of_det)):
    #    bins_of_det[i] = PA[vals_of_det[0][i],vals_of_det[1][i]]
    #    bins_of_moddet[i] = moAdet[vals_of_det[0][i],vals_of_det[1][i]]
    para = moAdet[idx] - PA[idx] + PA[idx]*(np.log1p(PA[idx])-np.log1p(moAdet[idx]))
    #para = moAdet[idx] - PA[idx] + PA[idx]*(np.log(PA[idx]) - np.log(moAdet[idx],where=moAdet[idx]>0))
    return (2*np.sum(para))

def fnA(params):
    x,y,z,w,t = params
    global PA,Aij,a0,a1,a2,a3
#    numerA = (PA-(x*Aij+y*a0+z*a1+w*a2+t*a3))**2
    moA = x*Aij+y*a0+z*a1+w*a2+t*a3
    for i in range(kA):
        moValsA = np.where(binA == i+1)
        bins_of_modA[i] = np.sum(moA[moValsA])
    numerA = bins_of_modA - bins_of_A + bins_of_A*(np.log(bins_of_A) - np.log(bins_of_modA))
    return (2 * np.sum(numerA))
# np.sum(np.divide(numerA,PA,out=np.zeros_like(numerA),where=PA!=0))


x0_d = np.zeros(4)
detsn = [a0,a1,a2,a3]

for k in range(4):
    resA = opt.minimize(detfn, x0_d[k],args = detsn[k],method='Nelder-Mead',tol=1e-06,options={'maxiter':50000,'maxfev':5000,'disp':True})
    x0_d[k] = resA.x
    with open(local_dir+'count_stat_01_nosun.txt', 'a+') as O:
        O.write('Telescope-'+detp+'; detector-'+str(k)+':'+'\n')
        O.write(str(resA)+'\n')
    
    #resdet = opt.minimize(detfn, x0_d[k],args = detsn[k],method='Nelder-Mead',options={'maxiter':5000,'disp':True})

# The next comments are for apgrad
#    x0_d[k] = resdet.x

#print(str(x0_d))
x0 = np.zeros(5)
x0[0] = 0.1
x0[1],x0[2],x0[3],x0[4] = x0_d

resA = opt.minimize(fnA,x0,method = 'Nelder-Mead', tol=1e-8, options = {'maxiter' : 5000, 'maxfev':5000,'disp' : True})

# The double ## indicate the active ones for sun/full version of counts stat

##with open(local_dir+'count_stat_printout.txt', 'w+') as O:
##    O.write('Telescope '+detp+':'+'\n')
##    O.write(str(resA)+'\n')

x,y,z,w,t = resA.x
modA = x*Aij+y*a0+z*a1+w*a2+t*a3

#x,y,z,t = x0_d
#modA = x*a0+y*a1+z*a2+t*a3

residA = PA - modA

hduAmod = fits.PrimaryHDU(modA)
hduAres = fits.PrimaryHDU(residA)


hduAmod.writeto(local_dir+detp+'_Model_01_nosun.fits')
hduAres.writeto(local_dir+detp+'_resid_01_nosun.fits')


#plt.matshow(PA)
#plt.matshow(modA)
#plt.matshow(residA)
#plt.show()

plt.plot(residA)
plt.savefig(local_dir+detp+'_resid_01_nosun.png')

#plt.show()

