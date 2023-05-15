#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# syntax: count_stat_new.py det file lowenergy highenergy sun/nosun/full

# Finds the statistics of an image for Mi,j = a*Ai,j + b1 + b2 + b3 + b4, where M is the model value, A is the physical location of the pixel, a is a variable to find, and b's are the base values of the individual dets.

# Tried both scipy optimize and curve_fit, however, without the ability to set boundries, it is 
# difficult to get reliable values out, and they seem to be very sensitive to the the initial values

# Now I am going to try lmfit
#  UPDATE:  
#          I need to redesign the code to upload the exposure map, divide the exp_map by the exposure and create a mask out of it.  For
#	   non excl.reg files it won't do anything really, but for those with, it will give the norm to the data file.
#	   Also, redo the binning image to create and actual image.


import os, sys, numpy as np
from astropy.io import fits
import scipy.optimize as opt
from scipy.optimize import curve_fit
import matplotlib

matplotlib.use('Agg')

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
from scipy.stats import norm


detp = sys.argv[1]
data = sys.argv[2]
lowe = sys.argv[3]
highe = sys.argv[4]
sep = sys.argv[5]
homedir = sys.argv[6]
offx = int(float(sys.argv[7]))
offy = int(float(sys.argv[8]))

def checkfile(fi):
  if os.path.isfile(fi) == False:
	return "BAD"
  else:
	return "GOOD"

if checkfile(data) == "BAD":
  print('Data file not found, please check file path')
  sys.exit()


##################
# files to use for counts:  pixmapA/B in nuskybgd auxil
#                           edgemaskA/B in auxil
#                           

# local files
nusky_dir = '/uufs/astro.utah.edu/common/home/u1019304/my_idl/nuskybgd/auxil/'
edge = '/uufs/astro.utah.edu/common/home/u1019304/temp/fullmask'+detp+'_final2.fits'
#edge = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/fullmask'+detp+'.fits'
#edgeB = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/fullmaskB.fits'
local_dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil/'

# Load data
with fits.open(data) as exp:
  Exposure = exp[0].header['EXPOSURE']

with fits.open('/uufs/astro.utah.edu/common/home/u1019304/temp/detmap'+detp+'.fits') as hA:
  dA = hA[0].data

with fits.open(edge) as hA:
  eA = hA[0].data  # To make a variable mask, make this mask a fraction of 1, this counts for the variable exposure time in excl.reg obs.

with fits.open(data) as hA:
  datA = hA[0].data

with fits.open(nusky_dir+'det'+detp+'_det1.img') as gA:
  gradA = gA[0].data

with fits.open('/uufs/astro.utah.edu/common/home/u1019304/temp/pixmap'+detp+'.fits') as pm:
  pixmap = pm[0].data  #check for new pix map image possibly make :)

#Load in the expmap. This should be the exp version of the data file.  Used on one line, ref: [jkl;] to find
exp_map = data.replace('Data','Exp')
with fits.open(exp_map) as hdul:
  exp_mask = hdul[0].data

eA_M = np.copy(exp_mask)
#eA_M /= Exposure
# may think about eA_M /= np.max(eA_M)
eA_M /= np.max(eA_M)
################################


# fix the Aij, Bij arrays 
if detp == 'A':
  Aij = gradA[(len(gradA)/2 - 180) - (5 + offy):(len(gradA)/2 + 180) - (5 + offy), (len(gradA)/2 - 180) + (8 + offx):(len(gradA)/2 + 180) + (8 + offx)]*eA
else:
  Aij = gradA[(len(gradA)/2 - 180) - (16 + offy):(len(gradA)/2 + 180) - (16 + offy), (len(gradA)/2 - 180) + (10 + offx):(len(gradA)/2 + 180) + (10 + offx)]*eA

# Here I am going to mess with the Aij to make it smaller in value

Aij *= 0.0001
#Aij *= 0.0001 # standard
# This is for a test in telescope A to see if I have teh same issues I am seeing

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
a3[vA == 4]=1 # To test if A3 is causing the issue of low points in the fitting, setting this to 0


#A_0 = Aij * a0
#A_1 = Aij * a1
#A_2 = Aij * a2
#A_3 = Aij * a3

##################

PA = np.copy(datA)*eA

#idx = np.where(PA != 0)
#print(np.min(PA[idx]))
#print(np.mean(PA[idx]))


#################
# This is where this script diverges from count_stat.py, we will bin in a 1-D array 
# First attempt: create a tuple of data values that I can append and if it meets the
# criteria, then it will be converted into a numpy array.

# Create a long 1-D array that will represent the gradiant bins
def get_bins():
  A_bins, count_bins = [], []
  count, Aij_count = 0, 0.
  kA = 0
  min_value = 10.0 # min value in each bin, if it's below this, it will combine bins
  max_bin_length = 100000 # maximum number of bins, if it's longer, we will up the min count value

  bin_map = np.zeros((360,360))  # This array will be used to create the image
  bin_idx = 0  # This value is used as a counter to map from the bins to the image

  det0_bins,det1_bins,det2_bins,det3_bins = [],[],[],[]
  d0count, d1count, d2count, d3count = 0,0,0,0

# Think about using ravel to flatten out the 2-d arrays.
# Part of sanity checks, I am switching the i and j indexing to see if it does something to the resulting images.
  for i in range(360):
    for j in range(360):
      bin_map[j,i] = bin_idx   # assigning bin values to pixels
      count += PA[j,i]     #
      Aij_count += Aij[j,i]*eA_M[j,i]  # NEED THE FRACTIONAL EXPOSURE ON THESE COUNT VALUES FOR MODEL ONLY!!!!!!!!
    ### keeping count of the counts in active dets pixels
      d0count += a0[j,i]*eA_M[j,i]
      d1count += a1[j,i]*eA_M[j,i]
      d2count += a2[j,i]*eA_M[j,i]
      d3count += a3[j,i]*eA_M[j,i]
      if count >= min_value:
	bin_idx += 1                 # For bin map index increase
        count_bins.append(count)
        A_bins.append(Aij_count)
        count,Aij_count = 0, 0.
        kA += 1
        det0_bins.append(d0count)
        det1_bins.append(d1count)
        det2_bins.append(d2count)
        det3_bins.append(d3count)
        d0count,d1count,d2count,d3count = 0,0,0,0
      if (i == 359) and (j == 359) and (count < min_value):  # Fix aij issue, not a value to put into the value
      #Aij_count += A_bins[-1:]  # If the last value in the array below the min count, it will
        A_bins[-1:] += Aij_count   # be added to the last bin. This may want to be redone
        det3_bins[-1:] += d3count  # it is also a det3 bin at this point
      if len(A_bins) > max_bin_length:
        i, j, count, Aij_count, d0count, d1count, d2count, d3count = 0,0,0, 0., 0,0,0,0
        count_bins, A_bins, det0_bins, det1_bins, det2_bins, det3_bins = [],[],[],[],[],[]
        min_value += 10      
        kA = 0
  return kA, det0_bins, det1_bins, det2_bins, det3_bins, A_bins, count_bins, min_value, bin_map
################  Removed code from loop
'''               ##################  WARNING  #####################
      if a0[j,i] == 1:   # All of this is to make the det 1-D arrays Need to count how many times the det array is on
        det0_bins.append()
      else: det0_bins.append(0)
      if a1[j,i] == 1:
        det1_bins.append(1)
      else: det1_bins.append(0)
      if a2[j,i] == 1:
        det2_bins.append(1)
      else: det2_bins.append(0)
      if a3[j,i] == 1:
        det3_bins.append(1)
      else: det3_bins.append(0)'''

################
# Need to find a way to seperate the bins into det binned arrays
# Think I found a way, it looks at the masks for each det and if 
# it exists in thier array, the count goes into there.

kA, det0_bins, det1_bins, det2_bins, det3_bins, A_bins, count_bins, min_value, bin_map = get_bins()



#print(min_value, len(A_bins), kA)
#print(min(A_bins))
# Check to see if the length of bins is consistant

if len(A_bins) != kA:
  print('inconsistant bin length')

if len(A_bins) != len(count_bins):
  print('A and count bins are of unequal length')

# Check for zeros in the bins

if A_bins.count(0) > 0:
  print("There are zeros in the bins")
  print(A_bins.count(0))
  print(A_bins.count(2.0))
  sys.exit()

# make it into a numpy array

A_bins = np.asarray(A_bins)
count_bins = np.asarray(count_bins)
det0_bins = np.asarray(det0_bins)
det1_bins = np.asarray(det1_bins)
det2_bins = np.asarray(det2_bins)
det3_bins = np.asarray(det3_bins)

######################
# Function to be passed to the optimization function
######################
def fn(params):
#  x,y,z,w,t = params
  x = params['aCXB']
  y = params['det0']
  z = params['det1']
  w = params['det2']
  t = params['det3']
#  global A_bins,det0_bins,det1_bins,det2_bins,det3_bins 
  mod = (x*A_bins) + (y*det0_bins) + (z*det1_bins) + (w*det2_bins) + (t*det3_bins)
  numerA = mod - count_bins + count_bins*(np.log(count_bins)-np.log(mod))
  return (2 * np.sum(numerA))

# This is to find the error values around the points, this is done
# by holding the aCXB norm constant (x_cerror) while allowing the other
# norms to be passed through the minimization function
def fn_cerror(params):
  y = params['det0']
  z = params['det1']
  w = params['det2']
  t = params['det3']
  mod = (x_c*A_bins) + (y*det0_bins) + (z*det1_bins) + (w*det2_bins) + (t*det3_bins)
  numerA = mod - count_bins + count_bins*(np.log(count_bins)-np.log(mod))
  return (2 * np.sum(numerA))
####################################

def detfn(param, detbin):
  mod = param*detbin
  idx = np.where((count_bins*detbin) != 0)
  numer = mod[idx] - count_bins[idx] + count_bins[idx]*(np.log(count_bins[idx])-np.log(mod[idx]))
  return(2 * np.sum(numer))

def f(X,x,y,z,w,t):
  x0 = X[:,0]
  x1 = X[:,1]
  x2 = X[:,2]
  x3 = X[:,3]
  x4 = X[:,4] 
  return (x*x0) + (y*x1) + (z*x2) + (w*x3) + (t*x4)

def residual(p):
  return (2 * np.sum(p-count_bins+count_bins*(np.log(count_bins)-np.log(p))))

def lmmin(pars,data=None):
  mod = pars['aCXB']*A_bins + pars['det0']*det0_bins + pars['det1']*det1_bins + pars['det2']*det2_bins + pars['det3']*det3_bins
  if data is None:
    return mod
  return mod - data
#  return (2*np.sum(mod-count_bins+count_bins*(np.log(count_bins)-np.log(mod))))
  
######################

# I will try the lmfit.minimize method
#if os.path.isfile('/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/CXB/01/nosun/3.0_40.16_71steps/'+detp+'_nosun_'+lowe+'_'+highe+'_keVparams_longformat.txt'):
#  print('TRUE')
#else: 
#  print('FALSE')

#f = open('/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/CXB/01/nosun/3.0_40.16_71steps/old_errors/'+detp+'_nosun_'+lowe+'_'+highe+'_keVparams_longformat.txt', 'r+')
  

#for line in f.readlines():
#  if 'aCXB norm' in line:
#    aCXB_initval = line.split()[3]
#  if 'det0 norm' in line:
#    det0_initval = line.split()[3]
#  if 'det1 norm' in line:
#    det1_initval = line.split()[3]
#  if 'det2 norm' in line:
#    det2_initval = line.split()[3]
#  if 'det3 norm' in line:
#    det3_initval = line.split()[3]
   
#f.close()
with open('/uufs/astro.utah.edu/common/home/u1019304/temp/normmeanvals.txt', 'r+') as f:
  lines = f.readlines()

aCXB_initval, det0_initval, det1_initval, det2_initval, det3_initval = lines[0].split()


aCXB_initval = float(aCXB_initval)
det0_initval = float(det0_initval)
det1_initval = float(det1_initval)
det2_initval = float(det2_initval)
det3_initval = float(det3_initval)

aCXB_initval *= Exposure; aCXB_initval *= (float(highe)-float(lowe))
det0_initval *= Exposure; det0_initval *= (float(highe)-float(lowe))
det1_initval *= Exposure; det1_initval *= (float(highe)-float(lowe))
det2_initval *= Exposure; det2_initval *= (float(highe)-float(lowe))
det3_initval *= Exposure; det3_initval *= (float(highe)-float(lowe))

fit_params = Parameters()
fit_params.add('aCXB', value=aCXB_initval, min=0.0)
fit_params.add('det0', value=det0_initval, min=0.0)
fit_params.add('det1', value=det1_initval, min=0.0)
fit_params.add('det2', value=det2_initval, min=0.0)
fit_params.add('det3', value=det3_initval, min=0.0)

#fit_params = Parameters()
#fit_params.add('aCXB', value=0.1, min=0.0)
#fit_params.add('det0', value=1.0, min=0.0)
#fit_params.add('det1', value=1.0, min=0.0)
#fit_params.add('det2', value=1.0, min=0.0)
#fit_params.add('det3', value=1.0, min=0.0)

#print(fit_params['aCXB'].value,fit_params['det0'].value,fit_params['det1'].value,fit_params['det2'].value,fit_params['det3'].value)


#Y = np.column_stack([A_bins, det0_bins, det1_bins, det2_bins, det3_bins])

#print(len(fit_params), len(Y))
'''
out = minimize(lmmin, fit_params, kws={'data': count_bins}) 
fit = lmmin(out.params)

report_fit(out, show_correl=True)
'''
'''
# So here is where I will be trying the scipy least squares:
bnds = ([0.0,1000.],[0.0,1000.],[0.0,1000.],[0.0,1000.],[0.0,1000.])
initial_guess = [0.1,5.0,5.0,5.0,5.0]
Xval = np.column_stack([A_bins,det0_bins,det1_bins,det2_bins,det3_bins])
popt, pcov = opt.curve_fit(f, Xval, count_bins, p0=initial_guess, bounds = bnds)

print(popt)
'''

'''
x,y,z,w,t = popt
#fullvall = np.column_stack([Aij,a0,a1,a2,a3])
mod = x*Aij+y*a0+z*a1+w*a2+t*a3
plt.matshow(mod)
plt.show()

resid = PA - mod
plt.matshow(resid)
plt.show()
# Initial values for the guesses of parameters
x0 = np.zeros(5)
detsn = [det0_bins,det1_bins,det2_bins,det3_bins]




for k in range(4):
  resdet = opt.minimize(detfn, 10.0, args=detsn[k], method='Nelder-Mead',tol=1e-06,options={'maxiter':50000,'maxfev':5000,'disp':True})
  x0[k+1] = resdet.x
  print(resdet.x)
'''  


#x0 = np.zeros(5)
#x0[0] = 0.1
#x0[1] = 1.0
#x0[2] = 1.0
#x0[3] = 1.0
#x0[4] = 1.0

#bnds = ([0.0,1000.],[0.0,1000.],[0.0,1000.],[0.0,1000.],[0.0,1000.])
 #!!!!!!!!!!!!!!!!!#################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
resA = minimize(fn, fit_params, method = 'nelder', tol=1e-15)
 #!!!!!!!!!!!!!!!!!################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#resA = opt.minimize(fn, x0, method = 'Nelder-Mead',tol=1e-8,options= {'maxiter' : 5000, 'maxfev':5000,'disp' : True})
#resA = opt.minimize(fn,x0,method = 'L-BFGS-B', bounds = bnds)
"""
print(resA.x)

x,y,z,w,t = resA.x
"""

#model = x*Aij+y*a0+z*a1+w*a2+t*a3
#model = fit_params['aCXB'].value*Aij+fit_params['det0'].value*a0+fit_params['det1'].value*a1+fit_params['det2'].value*a2+fit_params['det3'].value*a3  # THIS NEEDS TO BE MULTI BY THE EXP IMAGE

cstat = fn(fit_params)
#print(cstat)
#resid = model - PA



yaxis = np.arange(0,len(count_bins))
resid_x = lmmin(fit_params, count_bins)
#print(np.std(resid_x))
#print(np.std(resid_x[:-1]))
#report_fit(resA, show_correl=False)

#resid[resid >= np.std(resid)*4] = 0
# Here i need to find the raw pixel value, and exclude all of those real pixel values that are in that raw value
# I'll have to find which values are 
# The file I need to load is pixmap in the auxil file of nuskybgd
# I need to remove pixels from data (PA)
# Plan:  find the det it is in and the raw pixel value of them
#        find the raw pixel value of those pixels and exclude them from that det in PA
#####################################
mult = 11
# old value was 7. Tried 9, but it was still a little much. now on 11
# I need to change the eA mask type to a bool mask type, this allows the comparison to work since it will only
# look at the "active" pixels in the PA matrix (those that are not excluded by the mask)

######
##mask = eA.astype('bool')
##idx_std_row, idx_std_col = np.where(np.abs(PA) >= np.sqrt(model)*mult) # why still np.sqrt? because we are looking at cstat mean
#idx_std_row, idx_std_col = np.where(np.abs(PA) >= np.sqrt(np.mean(PA[mask]))*mult) # instead of mean, use model value there!
#idx_std_row, idx_std_col = np.where(np.abs(resid) >= np.std(resid)*mult)
##Amatrix = [a0,a1,a2,a3]  # small count regime we need to capture the statistics of fluctuations
##for i in range(len(idx_std_row)):
##  id_A = int(dA[idx_std_row[i]][idx_std_col[i]])-1
##  pixmap_fil = pixmap*Amatrix[id_A]
##  idx_pix_row, idx_pix_col = np.where(pixmap_fil == pixmap[idx_std_row[i]][idx_std_col[i]])
##  if len(idx_pix_row) > 0:
##    newbinning = True
##    for j in range(len(idx_pix_row)):
##      PA[idx_pix_row[j]][idx_pix_col[j]] = 0
##      eA[idx_pix_row[j]][idx_pix_col[j]] = 0
##  else: newbinning = False
# The values are not being properly subtracted from data matrix, it's leaving LARGE areas of 0 that are then
# incorporated into the model and not just considered zero...also the subtraction is ditching too many points
######################################
#resid = model - PA

##########reprocess###############
# redo masks and new binning:
'''  I need to write out each mask '''

##if newbinning:
##  Aij *= eA
##  a0 *= eA
##  a1 *= eA
##  a2 *= eA
##  a3 *= eA
 # don't do a rebin, just remove the bins that have 0 counts i.e. np.where(bins == 0);
##  kA, det0_bins, det1_bins, det2_bins, det3_bins, A_bins, count_bins, min_value, bin_map = get_bins()
#print(min_value, len(A_bins), kA)
##  A_bins = np.asarray(A_bins)
##  count_bins = np.asarray(count_bins)
##  det0_bins = np.asarray(det0_bins)
##  det1_bins = np.asarray(det1_bins)
##  det2_bins = np.asarray(det2_bins)
##  det3_bins = np.asarray(det3_bins)

#bin_arrays = np.stack((count_bins,A_bins,det0_bins,det1_bins,det2_bins,det3_bins))

bin_file = homedir+'/BINS'+detp+'_'+'01_nosun'+'_'+str(lowe)+'_'+str(highe)+'keV.fits'
fits.writeto(bin_file, bin_map)
with fits.open(bin_file, mode = 'update') as hdul:
  hdul[0].header['COMMENT'] = 'bins image for energy of '+str(lowe)+' to '+str(highe)
  hdul[0].header['COMMENT'] = 'A total number of bins is: '+str(kA) 
  hdul.flush()
# new params
#resA = minimize(fn, fit_params, method = 'nelder', tol=1e-15)
#model = fit_params['aCXB']*Aij+fit_params['det0']*a0+fit_params['det1']*a1+fit_params['det2']*a2+fit_params['det3']*a3 # HERE TOO
#cstat = fn(fit_params)
#print('Cstat: '+str(cstat))
#yaxis = np.arange(0,len(count_bins))
#resid_x = lmmin(fit_params, count_bins)
#print(np.std(resid_x))
#print(np.std(resid_x[:-1]))

#report_fit(resA, show_correl=True)
##################################
#resid = (model - PA)
#################################
### Make a fits file of the model#
#################################

#model *= (eA * eA_M) # This line is to multiply the model by the updated mask and the exposure mask to create a model jkl; What the hell is this?


#fii = homedir+'/Model'+detp+'_'+'01'+'_'+sep+'_'+str(lowe)+'_'+str(highe)+'keV.fits'
#fits.writeto(fii, model)
#with fits.open(fii, mode='update') as hdul:
#  hdul[0].header['ACXBNORM'] = str(fit_params['aCXB'].value)
#  hdul[0].header['DET0NORM'] = str(fit_params['det0'].value)
#  hdul[0].header['DET1NORM'] = str(fit_params['det1'].value)
#  hdul[0].header['DET2NORM'] = str(fit_params['det2'].value)
#  hdul[0].header['DET3NORM'] = str(fit_params['det3'].value)
#  hdul[0].header['COMMENT'] = 'Model for '+detp+' Energy: '+str(lowe)+' '+str(highe)
#  hdul.flush()

#################################
#plt.matshow(eA)
#plt.show() 

def gaussian(x, a, b, c):
    return a*np.exp(-np.power(x - b, 2)/(2*np.power(c, 2)))

idx_for_hist = np.where(eA > 0)
#(mu, sigma) = norm.fit(resid[idx_for_hist])
#lnspc = np.linspace(np.min(resid[idx_for_hist]), np.max(resid[idx_for_hist]), len(np.ravel(resid[idx_for_hist])))

#lnspc = np.linspace(np.min(resid[idx_for_hist]), np.max(resid[idx_for_hist]), len(np.ravel(resid[idx_for_hist])))

#pars, cov = curve_fit(gaussian, lnspc, np.ravel(resid[idx_for_hist]))
#print(pars)
#pdf_g = norm.pdf(lnspc, mu, sigma)
#plt.plot(lnspc, pdf_g)
#plt.plot(lnspc,gaussian(lnspc,*pars))
#n, bins, patches = plt.hist(np.ravel(np.around(resid[idx_for_hist],decimals=2)), bins='auto', normed=True, facecolor='green')
#plt.title("Binned distribution of residual counts in "+lowe+" to "+highe+" keV")
#plt.xlabel("Counts per bin")
#plt.ylabel("# of bins")
#plt.savefig(homedir+"/hist_"+detp+"_"+lowe+"_"+highe+".png")
#plt.close()
#xvals = np.arange(0,360*360)
#plt.scatter(xvals,np.sort(np.ravel(resid)))
#plt.show()

#plt.matshow(model)
#plt.title("Background model for "+lowe+" to "+highe+" keV")
#plt.savefig(homedir+"/model_image_"+detp+"_"+sep+"_"+lowe+"_"+highe+"keV.png")
#plt.close()

#plt.matshow(resid)
#plt.title("Residual map for "+lowe+" to "+highe+" keV")
#plt.savefig(homedir+"/resid_image_"+detp+"_"+sep+"_"+lowe+"_"+highe+"keV.png")
#plt.close()
#np.savetxt('resid'+detp+'.txt',resid)


############### From here, it is a C-stat fit program ########### should move to own program


def fitcstat(x,y,z,w,t):
#  global A_bins,det0_bins,det1_bins,det2_bins,det3_bins 
  mod = (x*A_bins) + (y*det0_bins) + (z*det1_bins) + (w*det2_bins) + (t*det3_bins)
  numerA = mod - count_bins + count_bins*(np.log(count_bins)-np.log(mod))
  return (2 * np.sum(numerA))

stand = Parameters()
stand.add('aCXB', value = fit_params['aCXB'])
stand.add('det0', value = fit_params['det0'])
stand.add('det1', value = fit_params['det1'])
stand.add('det2', value = fit_params['det2'])
stand.add('det3', value = fit_params['det3'])

fit_params_errors = Parameters()
fit_params_errors.add('det0', value=fit_params['det0'].value, min=0.0)
fit_params_errors.add('det1', value=fit_params['det1'].value, min=0.0)
fit_params_errors.add('det2', value=fit_params['det2'].value, min=0.0)
fit_params_errors.add('det3', value=fit_params['det3'].value, min=0.0)

def resetParams(val1, val2, val3, val4):
  fit_params_errors['det0'].set(value = val1, min=0.0)
  fit_params_errors['det1'].set(value = val2, min=0.0)
  fit_params_errors['det2'].set(value = val3, min=0.0)
  fit_params_errors['det3'].set(value = val4, min=0.0)

midval = fit_params['aCXB']
xvals,x1,x2 = [],[],[]
yvals,y1,y2 = [],[],[]
x_c = fit_params['aCXB']
#y_c = fit_params['det0']
#z_c = fit_params['det1']
#w_c = fit_params['det2']
#t_c = fit_params['det3'] 
xvals.append(x_c.value)
yvals.append(cstat)
newcstat = np.copy(cstat)
m_cstat = np.copy(cstat)
mult_val = 10.0 
onesigup = 0
onesigdown = 0
newdetvals = []
########### While loop #################
while newcstat < mult_val+m_cstat:
#  x,y,z,w,t = params
  x_c += 0.0001
#  fit_params_errors = Parameters()
#  fit_params_errors.add('det0', value=fit_params['det0'].value, min=0.0)
#  fit_params_errors.add('det1', value=fit_params['det1'].value, min=0.0)
#  fit_params_errors.add('det2', value=fit_params['det2'].value, min=0.0)
#  fit_params_errors.add('det3', value=fit_params['det3'].value, min=0.0)
  resA = minimize(fn_cerror, fit_params_errors, method = 'nelder', tol=1e-15)
  newcstat = fitcstat(x_c,fit_params_errors['det0'].value,fit_params_errors['det1'].value,fit_params_errors['det2'].value,fit_params_errors['det3'].value) 
  # Need to build in an exception to catch lower cstat values if they exist -- reset cstat value to new min
  # Do i need to reset the parameter values before each run? ##############################
  if newcstat < m_cstat: 
	print('exception exisits '+data+' high')
	m_cstat = None
	m_cstat = newcstat
        if len(newdetvals) > 0: newdetvals = []
	newdetvals.append(fit_params_errors['det0'].value)
	newdetvals.append(fit_params_errors['det1'].value)
	newdetvals.append(fit_params_errors['det2'].value)
	newdetvals.append(fit_params_errors['det3'].value)
	newaCXB = x_c
        stand['aCXB'].set(value=newaCXB)
        stand['det0'].set(value=newdetvals[0])
        stand['det1'].set(value=newdetvals[1])
        stand['det2'].set(value=newdetvals[2])
        stand['det3'].set(value=newdetvals[3])
	xvals,yvals,x1,y1 =  [],[],[],[]
	xvals.append(x_c.value)
	yvals.append(m_cstat)
	blah = open(homedir+'/'+detp+"_high_params_errors.txt",'w+')
 	blah.close()
	continue
  xvals.append(float(x_c))
  yvals.append(float(newcstat))
  x1.append(float(x_c))
  y1.append(float(newcstat))
  with open(homedir+'/'+detp+"_high_params_errors.txt",'a+') as O:
    O.write(str(x_c)+' '+str(fit_params_errors['det0'].value)+' '+str(fit_params_errors['det1'].value)+' '+str(fit_params_errors['det2'].value)+' '+str(fit_params_errors['det3'].value)+'\n')
#  fit_params_errors=None

if len(newdetvals) > 0:
	resetParams(newdetvals[0],newdetvals[1],newdetvals[2],newdetvals[3])
	x_c = None
	x_c = newaCXB
	newcstat = np.copy(m_cstat)
else:
	fit_params_errors = None
	fit_params_errors = Parameters()
	fit_params_errors.add('det0', value=fit_params['det0'].value, min=0.0)
	fit_params_errors.add('det1', value=fit_params['det1'].value, min=0.0)
	fit_params_errors.add('det2', value=fit_params['det2'].value, min=0.0)
	fit_params_errors.add('det3', value=fit_params['det3'].value, min=0.0)
	newcstat=np.copy(cstat)
	x_c = fit_params['aCXB']
xvals_2 = xvals[:]
yvals_2 = yvals[:]

while newcstat < mult_val+m_cstat:
  x_c -= 0.0001
#  fit_params_errors = Parameters()
#  fit_params_errors.add('det0', value=fit_params['det0'].value, min=0.0)
#  fit_params_errors.add('det1', value=fit_params['det1'].value, min=0.0)
#  fit_params_errors.add('det2', value=fit_params['det2'].value, min=0.0)
#  fit_params_errors.add('det3', value=fit_params['det3'].value, min=0.0)
  resA = minimize(fn_cerror, fit_params_errors, method = 'nelder', tol=1e-15)
  newcstat = fitcstat(x_c,fit_params_errors['det0'].value,fit_params_errors['det1'].value,fit_params_errors['det2'].value,fit_params_errors['det3'].value)
  # Need to build in an exception to catch lower cstat values if they exist
  if newcstat < m_cstat: 
	print('exception exisits '+data+' low')
	m_cstat = None
        m_cstat = newcstat
        if len(newdetvals) > 0: newdetvals = []
        newdetvals.append(fit_params_errors['det0'].value)
        newdetvals.append(fit_params_errors['det1'].value)
        newdetvals.append(fit_params_errors['det2'].value)
        newdetvals.append(fit_params_errors['det3'].value)
        newaCXB = x_c
        stand['aCXB'].set(value=newaCXB)
        stand['det0'].set(value=newdetvals[0])
        stand['det1'].set(value=newdetvals[1])
        stand['det2'].set(value=newdetvals[2])
        stand['det3'].set(value=newdetvals[3])
        x2,y2 =  [],[]
        xvals = xvals_2[:]
        yvals = yvals_2[:]
        blah = open(homedir+'/'+detp+"_low_params_errors.txt",'w+')
        blah.close()
	continue
  xvals.append(float(x_c))
  yvals.append(float(newcstat))
  y2.append(float(newcstat))
  x2.append(float(x_c))
  with open(homedir+'/'+detp+"_low_params_errors.txt",'a+') as O:
    O.write(str(x_c)+' '+str(fit_params_errors['det0'].value)+' '+str(fit_params_errors['det1'].value)+' '+str(fit_params_errors['det2'].value)+' '+str(fit_params_errors['det3'].value)+'\n')

fi_high = homedir+'/'+detp+"_high_params_errors.txt"
fi_low = homedir+'/'+detp+"_low_params_errors.txt"
fwrite = homedir+'/'+detp+"_params_errors.txt"
with open(fi_high) as fh:
	lines = fh.readlines()
	with open(fwrite,'w') as fw:
		fw.write(str(lines))
with open(fi_low) as fl:
	lines = fl.readlines()
	with open(fwrite,'a') as fw:
		fw.write(str(lines))

os.system('rm '+fi_high+' '+fi_low) 

# if there are new lower cstat values, I will probably have to redo some images and plots. I will have to redo some files also. so check.


#  fit_params_errors = None

################################################################

######## Functions for errors #######################################
def parabola(x,a,b,c):
  return a*x**2+b*x+c

def solu(cs,a,b,cv):
  c = cv-cs-1
  pos = (-b+np.sqrt(b**2-4*a*c))/(2*a)
  neg = (-b-np.sqrt(b**2-4*a*c))/(2*a)
  return max(pos,neg)

def sold(cs,a,b,cv):
  c = cv-cs-1
  pos = (-b+np.sqrt(b**2-4*a*c))/(2*a)
  neg = (-b-np.sqrt(b**2-4*a*c))/(2*a)
  return min(pos,neg)
#####################################################################
#for i in range(len(xvals)):
#  with open('parabola.txt','a+') as O:
#    O.write(str(xvals[i])+' '+str(yvals[i])+'\n')

####### Curve fit #################################
p,c = curve_fit(parabola,xvals,yvals)

#p = [0.005,.01,10000]
pars1, cov1 = opt.curve_fit(parabola,x1,y1,p)
pars2, cov2 = opt.curve_fit(parabola,x2,y2,p)

yfit1,yfit2 = [],[]

for i in range(len(x1)): 
    fitvalu = parabola(x1[i],*pars1)
    yfit1.append(fitvalu)

for i in range(len(x2)): 
    fitvalu = parabola(x2[i],*pars2)
    yfit2.append(fitvalu)
#################################################

############ 1-sig error values ###############
sigup = solu(m_cstat,*pars1)
sigdown = sold(m_cstat,*pars2)    
###############################################

plt.scatter(xvals,yvals,s=1)
plt.plot(x1,yfit1,label = 'fit1',linewidth=2)
plt.plot(x2,yfit2,label = 'fit2',linewidth=2)
plt.axvline(x=fit_params['aCXB'])
plt.title("Cstat error values for 10 sigma in "+lowe+" to "+highe+" keV")
plt.xlabel("aCXB")
plt.ylabel(r'$\Delta$'+"Cstat")
plt.legend()
#plt.show()
plt.savefig(homedir+"/Cstat_errors_"+detp+"_"+sep+"_"+lowe+"_"+highe+".png")
plt.close()
# save the plots and fit a parabola to both sides of the graph to recall exact values. 

############ Save Values to file ###################

# The format of the file will be the aCXB and det values with notes on the gradiant values used
# next will be the parameters for the parabolic fit with the range of xvals

with open(homedir+'/'+detp+"_params.txt",'a+') as O:
  O.write(lowe+' '+highe+' '+str(fit_params['aCXB'].value)+' '+str(sigup)+' '+str(sigdown)+' '+str(Exposure)+'\n')

with open(homedir+'/'+detp+'_'+sep+'_'+lowe+'_'+highe+'_'+'keVparams_longformat.txt','a+') as O: 
  O.write(" Gradiant file: "+nusky_dir+'det'+detp+'_det1.img divided by 0.0001'+'\n')
  O.write("Bins: "+str(len(count_bins))+" Total count:"+str(np.sum(count_bins))+'\n')
  O.write("aCXB norm value: "+str(fit_params['aCXB'].value)+'\n')
  O.write("det0 norm value: "+str(fit_params['det0'].value)+'\n')
  O.write("det1 norm value: "+str(fit_params['det1'].value)+'\n')
  O.write("det2 norm value: "+str(fit_params['det2'].value)+'\n')
  O.write("det3 norm value: "+str(fit_params['det3'].value)+'\n')
  O.write("#"*25+'\n')
  O.write("The following fit parameters are for the error values based on a parabola ax^2+bx+c"+'\n')
  O.write("X-range: "+str(min(xvals))+" - "+str(max(xvals))+'\n')
  O.write("Min Cstat value: "+str(cstat)+'\n')
  O.write("The full range: "+str(p)+'\n')
  O.write("Params for <= aCXB value: "+str(pars1)+'\n')
  O.write("params for > aCXB value: "+str(pars2)+'\n')

model_data = (stand['aCXB'].value*Aij+stand['det0'].value*a0+stand['det1'].value*a1+stand['det2'].value*a2+stand['det3'].value*a3)*eA_M  # THIS NEEDS TO BE MULTI BY THE EXP IMAGE
model = (stand['aCXB'].value*Aij+stand['det0'].value*a0+stand['det1'].value*a1+stand['det2'].value*a2+stand['det3'].value*a3)

fiii = homedir+'/Model_norm_'+detp+'_'+'01'+'_'+sep+'_'+str(lowe)+'_'+str(highe)+'keV.fits'
fii = homedir+'/Model'+detp+'_'+'01'+'_'+sep+'_'+str(lowe)+'_'+str(highe)+'keV.fits'
fits.writeto(fii, model_data)
with fits.open(fii, mode='update') as hdul:
  hdul[0].header['ACXBNORM'] = str(stand['aCXB'].value)
  hdul[0].header['DET0NORM'] = str(stand['det0'].value)
  hdul[0].header['DET1NORM'] = str(stand['det1'].value)
  hdul[0].header['DET2NORM'] = str(stand['det2'].value)
  hdul[0].header['DET3NORM'] = str(stand['det3'].value)
  hdul[0].header['COMMENT'] = 'Model for '+detp+' Energy: '+str(lowe)+' '+str(highe)
  hdul.flush()

fits.writeto(fiii, model)
with fits.open(fiii, mode='update') as hdul:
  hdul[0].header['ACXBNORM'] = str(stand['aCXB'].value)
  hdul[0].header['DET0NORM'] = str(stand['det0'].value)
  hdul[0].header['DET1NORM'] = str(stand['det1'].value)
  hdul[0].header['DET2NORM'] = str(stand['det2'].value)
  hdul[0].header['DET3NORM'] = str(stand['det3'].value)
  hdul[0].header['COMMENT'] = 'Model for '+detp+' Energy: '+str(lowe)+' '+str(highe)
  hdul.flush()

resid = model_data - PA

plt.matshow(model)
plt.title("Background model for "+lowe+" to "+highe+" keV")
plt.savefig(homedir+"/model_image_"+detp+"_"+sep+"_"+lowe+"_"+highe+"keV.png")
plt.close()

plt.matshow(resid)
plt.title("Residual map for "+lowe+" to "+highe+" keV")
plt.savefig(homedir+"/resid_image_"+detp+"_"+sep+"_"+lowe+"_"+highe+"keV.png")
plt.close()

n, bins, patches = plt.hist(np.ravel(np.around(resid[idx_for_hist],decimals=2)), bins='auto', normed=True, facecolor='green')
plt.title("Binned distribution of residual counts in "+lowe+" to "+highe+" keV")
plt.xlabel("Counts per bin")
plt.ylabel("# of bins")
plt.savefig(homedir+"/hist_"+detp+"_"+lowe+"_"+highe+".png")
plt.close()
