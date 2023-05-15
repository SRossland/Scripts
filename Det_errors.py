#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# syntax: Det_errors.py  data_file bins_file detmap_file expmap_file aCXBmap_file norm_file

# Purpose:  This program is to find the error values associated with the det norms found by the count_stat_new.py program
# by reintroducing the fits files for the Data and the Bins. 

import numpy as np, os, sys, matplotlib.pyplot as plt
from astropy.io import fits
import scipy.optimize as opt
from scipy.optimize import curve_fit, fsolve
from lmfit import Parameters, minimize

data_file = sys.argv[1]
bins_file = sys.argv[2]
det_map_file = sys.argv[3]
exp_map_file = sys.argv[4]
Aij_file = sys.argv[5]
norm_file = sys.argv[6]
offx = 0  #int(sys.argv[7])
offy = 0  #int(sys.argv[8])

data_split = data_file.split('/')[-1]
lowe = data_split.split('_')[3]
hig = data_split.split('_')[4]
highe, crap = hig.split('keV.fits')

detp = data_split.split('_')[0][4]

############# Definitions ###################

# C-stat fit

def fit_Cstat(x,y,z,w,t):
  mod = ((x*A_bins) + (y*det0_bins) + (z*det1_bins) + (w*det2_bins) + (t*det3_bins))
  numerA = mod - count_bins + count_bins*(np.log(count_bins)-np.log(mod))
  return (2 * np.sum(numerA))


def resetParams(val0, val1, val2, val3, val4):
  fit_params['aCXB'].set(value = val0)
  fit_params['det0'].set(value = val1, min=0.0)
  fit_params['det1'].set(value = val2, min=0.0)
  fit_params['det2'].set(value = val3, min=0.0)
  fit_params['det3'].set(value = val4, min=0.0)

def fn_cerror(para):
  y = para['det0']
  z = para['det1']
  w = para['det2']
  t = para['det3']
  mod = ((x_c*A_bins) + (y*det0_bins) + (z*det1_bins) + (w*det2_bins) + (t*det3_bins))
  numerA = mod - count_bins + count_bins*(np.log(count_bins)-np.log(mod))
  return (2 * np.sum(numerA))

def fn_cerror_det0(para):
  x = para['aCXB']
  z = para['det1']
  w = para['det2']
  t = para['det3']
  mod = ((x*A_bins) + (y_c*det0_bins) + (z*det1_bins) + (w*det2_bins) + (t*det3_bins)) 
  numerA = mod - count_bins + count_bins*(np.log(count_bins)-np.log(mod))
  return (2 * np.sum(numerA))

def fn_cerror_det1(para):
  x = para['aCXB']
  y = para['det0']
  w = para['det2']
  t = para['det3']
  mod = ((x*A_bins) + (y*det0_bins) + (z_c*det1_bins) + (w*det2_bins) + (t*det3_bins)) 
  numerA = mod - count_bins + count_bins*(np.log(count_bins)-np.log(mod))
  return (2 * np.sum(numerA))

def fn_cerror_det2(para):
  x = para['aCXB']
  z = para['det1']
  y = para['det0']
  t = para['det3']
  mod = ((x*A_bins) + (y*det0_bins) + (z*det1_bins) + (w_c*det2_bins) + (t*det3_bins)) 
  numerA = mod - count_bins + count_bins*(np.log(count_bins)-np.log(mod))
  return (2 * np.sum(numerA))

def fn_cerror_det3(para):
  x = para['aCXB']
  z = para['det1']
  w = para['det2']
  y = para['det0']
  mod = ((x*A_bins) + (y*det0_bins) + (z*det1_bins) + (w*det2_bins) + (t_c*det3_bins)) 
  numerA = mod - count_bins + count_bins*(np.log(count_bins)-np.log(mod))
  return (2 * np.sum(numerA))

def setParameters(name1, var1, name2, var2, name3, var3, name4, var4):
  global fit_params
  fit_params = None
  fit_params = Parameters()
  fit_params.add(name1, value=var1)
  fit_params.add(name2, value=var2, min = 0.0)
  fit_params.add(name3, value=var3, min = 0.0)
  fit_params.add(name4, value=var4, min = 0.0) 

def resetaCXBparams(val1, val2, val3, val4):
  fit_params_acxb['det0'].set(value=val1, min=0.0)
  fit_params_acxb['det1'].set(value=val2, min=0.0)
  fit_params_acxb['det2'].set(value=val3, min=0.0)
  fit_params_acxb['det3'].set(value=val4, min=0.0)
# I need a function that can run to find the error


def getData(fi):
  with fits.open(fi) as hdul:
    var = hdul[0].data
  return var

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

def runNewaCXB():
  #print('NEW aCXB NORM:')
  global newcstat, m_cstat, x_c, fit_params_acxb
  newcstat = np.copy(m_cstat)
  newdetvals_a = []
  fit_params_acxb = None
  fit_params_acxb = Parameters()
  fit_params_acxb.add('det0',value=params['det0'].value, min=0.0)
  fit_params_acxb.add('det1',value=params['det1'].value, min=0.0)
  fit_params_acxb.add('det2',value=params['det2'].value, min=0.0)
  fit_params_acxb.add('det3',value=params['det3'].value, min=0.0)
  xvals_a,x1_a,x2_a = [],[],[]
  yvals_a,y1_a,y2_a = [],[],[]
  x_c = np.copy(fit_params['aCXB'].value)
  newaCXB = np.copy(x_c)
  onesigup = 0
  onesigdown = 0
  while newcstat < mult_val+m_cstat:
    x_c += 0.0001
    #print(fit_params_acxb['det0'].value)
    resA = minimize(fn_cerror, fit_params_acxb, method = 'nelder', tol=1e-15)
    #print(fit_params_acxb['det0'].value)
    newcstat = fit_Cstat(x_c,fit_params_acxb['det0'].value,fit_params_acxb['det1'].value,fit_params_acxb['det2'].value,fit_params_acxb['det3'].value)
    #print(newcstat,m_cstat)
    if newcstat < m_cstat:
        #print('newvalue')
        m_cstat = None
        m_cstat = np.copy(newcstat)
        if len(newdetvals_a) > 0: newdetvals_a = []
        newdetvals_a.append(fit_params_acxb['det0'].value)
        newdetvals_a.append(fit_params_acxb['det1'].value)
        newdetvals_a.append(fit_params_acxb['det2'].value)
        newdetvals_a.append(fit_params_acxb['det3'].value)
        newaCXB = np.copy(x_c)
        xvals_a,yvals_a,x1_a,y1_a =  [],[],[],[]
        xvals_a.append(x_c)
        yvals_a.append(m_cstat)
        blah = open("aCXB_high_params_errors.txt",'w+')
        blah.close()
        continue
    xvals_a.append(float(x_c))
    yvals_a.append(float(newcstat))
    x1_a.append(float(x_c))
    y1_a.append(float(newcstat))
    with open("aCXB_high_params_errors.txt",'a+') as O:
      O.write(str(x_c)+' '+str(fit_params_acxb['det0'].value)+' '+str(fit_params_acxb['det1'].value)+' '+str(fit_params_acxb['det2'].value)+' '+str(fit_params_acxb['det3'].value)+'\n')

  if len(newdetvals_a) > 0:
        resetaCXBparams(newdetvals_a[0],newdetvals_a[1], newdetvals_a[2], newdetvals_a[3])
        x_c = None
        x_c = np.copy(newaCXB)
        newcstat = np.copy(m_cstat)
  else:
        fit_params_acxb = None
        fit_params_acxb = Parameters()
        fit_params_acxb.add('det0',value=params['det0'].value, min=0.0)
        fit_params_acxb.add('det1',value=params['det1'].value, min=0.0)
        fit_params_acxb.add('det2',value=params['det2'].value, min=0.0)
        fit_params_acxb.add('det3',value=params['det3'].value, min=0.0)
        x_c = None
        x_c = fit_params['aCXB'].value
        newcstat = np.copy(m_cstat)
  xvals_2_a = xvals_a[:]
  yvals_2_a = yvals_a[:]
  newdetvals_a = []
  while newcstat < mult_val+m_cstat:
    x_c -= 0.0001
    resA = minimize(fn_cerror, fit_params_acxb, method = 'nelder', tol=1e-15)
    newcstat = fit_Cstat(x_c,fit_params_acxb['det0'].value,fit_params_acxb['det1'].value,fit_params_acxb['det2'].value,fit_params_acxb['det3'].value)
    if newcstat < m_cstat:
        m_cstat = None
        m_cstat = np.copy(newcstat)
        if len(newdetvals_a) > 0: newdetvals_a = []
        newdetvals_a.append(fit_params_acxb['det0'].value)
        newdetvals_a.append(fit_params_acxb['det1'].value)
        newdetvals_a.append(fit_params_acxb['det2'].value)
        newdetvals_a.append(fit_params_acxb['det3'].value)
        newaCXB = np.copy(x_c)
        x2_a,y2_a =  [],[]
        xvals_a = xvals_2_a[:]
        yvals_a = yvals_2_a[:]
        blah = open("aCXB_low_params_errors.txt",'w+')
        blah.close()
        continue
    xvals_a.append(float(x_c))
    yvals_a.append(float(newcstat))
    y2_a.append(float(newcstat))
    x2_a.append(float(x_c))
    #print(str(newcstat))
    with open("aCXB_low_params_errors.txt",'a+') as O:
      O.write(lowe+' '+highe+' '+str(x_c)+' '+str(fit_params_acxb['det0'].value)+' '+str(fit_params_acxb['det1'].value)+' '+str(fit_params_acxb['det2'].value)+' '+str(fit_params_acxb['det3'].value)+'\n')

  if len(newdetvals_a) > 0:
        params['aCXB'].set(value = newaCXB)
        params['det0'].set(value = newdetvals_a[0], min=0.0)
        params['det1'].set(value = newdetvals_a[1], min=0.0)
        params['det2'].set(value = newdetvals_a[2], min=0.0)
        params['det3'].set(value = newdetvals_a[3], min=0.0)

  fi_high_a = "aCXB_high_params_errors.txt"
  fi_low_a = "aCXB_low_params_errors.txt"
  fwrite_a = detp+"_newaCXB_params_errors.txt"

  with open(fi_high_a) as fh:
        lines = fh.readlines()
        with open(fwrite_a,'a') as fw:
                fw.write(str(lines))
  with open(fi_low_a) as fl:
        lines = fl.readlines()
        with open(fwrite_a,'a') as fw:
                fw.write(str(lines))

  os.system('rm '+fi_high_a+' '+fi_low_a)

  p,c = curve_fit(parabola,xvals_a,yvals_a)
  pars1, cov1 = opt.curve_fit(parabola,x1_a,y1_a,p)
  pars2, cov2 = opt.curve_fit(parabola,x2_a,y2_a,p)
  sigup = solu(m_cstat,*pars1)
  sigdown = sold(m_cstat,*pars2)

  with open(detp+"_aCXBtemp_params.txt",'a+') as O:
    O.write(lowe+' '+highe+' '+str(newaCXB)+' '+str(sigup)+' '+str(sigdown)+' '+str(exp)+'\n')

  #print('done')
  return m_cstat

def getDeciPlace(num):
  s = '{:.16f}'.format(num).split('.')[1]
  return len(s) - len(s.lstrip('0')) + 1


def findStep(x):
  if x == 0:
    root = fsolve(StepCstat0,y_c)
    ste = (np.abs(root) - y_c)/20
    #print(y_c,ste)
  if x == 1:
    root = fsolve(StepCstat1,z_c)
    ste = (np.abs(root) - z_c)/20
  if x == 2:
    root = fsolve(StepCstat2,w_c)
    ste = (np.abs(root) - w_c)/20
  if x == 3:
    root = fsolve(StepCstat3,t_c)
    ste = (np.abs(root) - t_c)/20
  if float(ste) == 0: ste = 0.00001
  st = 10**-getDeciPlace(float(ste))
  return st

def StepCstat0(step_val):
  global m_cstat, mult_val, x_c
  mod = ((x_c*A_bins) + (step_val*det0_bins) + (z_c*det1_bins) + (w_c*det2_bins) + (t_c*det3_bins))
  numerA = mod - count_bins + count_bins*(np.log(count_bins)-np.log(mod))
  return ((2 * np.sum(numerA))-(m_cstat+mult_val))

def StepCstat1(step_val):
  global m_cstat, mult_val, x_c
  mod = ((x_c*A_bins) + (y_c*det0_bins) + (step_val*det1_bins) + (w_c*det2_bins) + (t_c*det3_bins))
  numerA = mod - count_bins + count_bins*(np.log(count_bins)-np.log(mod))
  return ((2 * np.sum(numerA))-(m_cstat+mult_val))

def StepCstat2(step_val):
  global m_cstat, mult_val, x_c
  mod = ((x_c*A_bins) + (y_c*det0_bins) + (z_c*det1_bins) + (step_val*det2_bins) + (t_c*det3_bins))
  numerA = mod - count_bins + count_bins*(np.log(count_bins)-np.log(mod))
  return ((2 * np.sum(numerA))-(m_cstat+mult_val))

def StepCstat3(step_val):
  global m_cstat, mult_val, x_c
  mod = ((x_c*A_bins) + (y_c*det0_bins) + (z_c*det1_bins) + (w_c*det2_bins) + (step_val*det3_bins))
  numerA = mod - count_bins + count_bins*(np.log(count_bins)-np.log(mod))
  return ((2 * np.sum(numerA))-(m_cstat+mult_val))

##############################################
# Need to try to load in the model image and see if I can pull those values that are zero
#mask_file = 'Model'+detp+'_01_nosun_'+lowe+'_'+highe+'keV.fits'
#if os.path.isfile(mask_file):
#  mod = getData(mask_file)
#elif os.path.isfile('fits_file/'+mask_file):
#  mod = getData('fits_file/'+mask_file)
#else: 
#  print('Model file not found, ensure it exists!')
#  sys.exit()

#mask = np.ones((360,360))
#mask[mod == 0] = 0

# Instead of doing what I am doing above, it would be easier to just take the mask used originally

mask_file = '/uufs/astro.utah.edu/common/home/u1019304/temp/fullmask'+detp+'_final.fits' # Hard coded file
mask = getData(mask_file)


# load data
# Data to load: count_bins, A_bins, det0-3 bins, e_m binned 
# I'll have to load det map, data, expmap, and bin map.

data = getData(data_file)
bins = getData(bins_file)
det_map = getData(det_map_file)
exp_map = getData(exp_map_file)
gradA = getData(Aij_file)

det_map = det_map*mask
exp_map = exp_map*mask
data = data*mask

if detp == 'A':
  Aij = gradA[(len(gradA)/2 - 180) - (1 + offy):(len(gradA)/2 + 180) - (1 - offy), (len(gradA)/2 - 180) + (4 + offx):(len(gradA)/2 + 180) + (4 + offx)]*mask
else:
  Aij = gradA[(len(gradA)/2 - 180) - (12 + offy):(len(gradA)/2 + 180) - (12 - offy), (len(gradA)/2 - 180) + (6 + offx):(len(gradA)/2 + 180) + (6 + offx)]*mask

#Aij[exp_map == 0] = 0  # need to double check this
Aij *= 0.0001

a0 = np.zeros((360,360)); a1 = np.zeros((360,360)); a2 = np.zeros((360,360)); a3 = np.zeros((360,360))
a0[det_map == 0] = 1
a1[det_map == 1] = 1
a2[det_map == 2] = 1
a3[det_map == 3] = 1

a0 *= mask; a1 *= mask; a2 *= mask; a3 *= mask


# need to fix the exp_map to make sure it is normalized

with fits.open(data_file) as hdul:
  exp = hdul[0].header['EXPOSURE']

expmask = np.copy(exp_map)
expmask /= exp

# FORGOT TO MULTIPLY THE ARRAYS BY THE EXPOSURE MAP:
############## This one line has been fucking with my life!!!!!
Aij *= expmask; a0 *= expmask; a1 *= expmask; a2 *= expmask; a3 *= expmask
##############
# rebin up the data, Acxb, dets, and exp according to the bins
# bins is an array from 0 to how many bins there are (-1) 

A_bins, det0_bins, det1_bins, det2_bins, det3_bins, e_m, count_bins = [],[],[],[],[],[],[]
for i in range(int(np.max(bins))+1):
  idx = np.where(bins == i)
  A_bins.append(np.sum(Aij[idx]))
  e_m.append(np.sum(expmask[idx]))
  count_bins.append(np.sum(data[idx]))
  det0_bins.append(np.sum(a0[idx]))
  det1_bins.append(np.sum(a1[idx]))
  det2_bins.append(np.sum(a2[idx]))
  det3_bins.append(np.sum(a3[idx]))
  

# if the last bin has a total count of 0 in the count bins, this removes that bin and adds it to the bin prior
if count_bins[-1] == 0:
  for matrix in [A_bins, det0_bins, det1_bins, det2_bins, det3_bins, e_m, count_bins]:
    var1 = matrix.pop()
    var2 = matrix.pop()
    matrix.append(var1+var2)

A_bins = np.asarray(A_bins)
det0_bins = np.asarray(det0_bins)
det1_bins = np.asarray(det1_bins)
det2_bins = np.asarray(det2_bins)
det3_bins = np.asarray(det3_bins)
e_m = np.asarray(e_m)
count_bins = np.asarray(count_bins)

idx_count = np.mean(count_bins)
idx_A = np.mean(A_bins)
idx_em = np.mean(e_m)


f = open(norm_file,'r+')
for line in f.readlines():
  if 'aCXB norm' in line:
    aCXB_val = line.split()[3]
  if 'det0 norm' in line:
    det0_val = line.split()[3]
  if 'det1 norm' in line:
    det1_val = line.split()[3]
  if 'det2 norm' in line:
    det2_val = line.split()[3]
  if 'det3 norm' in line:
    det3_val = line.split()[3]
f.close()

#print(aCXB_val, det0_val, det1_val, det2_val, det3_val)
# find the errors for det0 - det3

################ Globals ####################
params = Parameters()
params.add('aCXB', value=float(aCXB_val))
params.add('det0', value=float(det0_val), min=0.0)
params.add('det1', value=float(det1_val), min=0.0)
params.add('det2', value=float(det2_val), min=0.0)
params.add('det3', value=float(det3_val), min=0.0)
#print(det1_val)
#fit_params = Parameters()
#fit_params.add('aCXB', value=float(aCXB_val))
#fit_params.add('det0', value=float(det0_val), min=0.0)
#fit_params.add('det1', value=float(det1_val), min=0.0)
#fit_params.add('det2', value=float(det2_val), min=0.0)
#fit_params.add('det3', value=float(det3_val), min=0.0)


x_c = params['aCXB'].value; y_c = params['det0'].value; z_c = params['det1'].value; w_c = params['det2'].value; t_c = params['det3'].value
cstat = fit_Cstat(params['aCXB'].value,params['det0'].value,params['det1'].value,params['det2'].value,params['det3'].value)
newcstat = float(np.copy(cstat))
m_cstat = float(np.copy(cstat))
mult_val = 10.0
detarray = ['det0','det1','det2','det3']
#############################################

#### main error function ################
def doCstatRun():
  # define globals
  global y_c, z_c, w_c, t_c, newcstat, m_cstat, newdetvals, fit_params, x_c
  xvals, x1, x2 = [],[],[]
  yvals, y1, y2 = [],[],[]
  newdetvals = []
  fit_params = None
  for i in [0,1,2,3]:
    if i == 0: 
      setParameters('aCXB', params['aCXB'].value, 'det1', params['det1'].value*1.06, 'det2', params['det2'].value, 'det3', params['det3'].value)
    if i == 1: 
      setParameters('aCXB', params['aCXB'].value, 'det0', params['det0'].value, 'det2', params['det2'].value, 'det3', params['det3'].value)
    if i == 2: 
      setParameters('aCXB', params['aCXB'].value, 'det0', params['det0'].value, 'det1', params['det1'].value, 'det3', params['det3'].value)
    if i == 3: 
      setParameters('aCXB', params['aCXB'].value, 'det0', params['det0'].value, 'det1', params['det1'].value, 'det2', params['det2'].value)

    for direction in ['high','low']:
      if i == 0:
        xvals.append(float(y_c))
        if direction == 'high': x1.append(float(y_c))
        if direction == 'low': x2.append(float(y_c))
#	fit_params = Parameters()
#	fit_params.add('aCXB', value=float(aCXB_val))
	#fit_params.add('det0', value=float(det0_val), min=0.0)
#	fit_params.add('det1', value=float(det1_val), min=0.0)
#	fit_params.add('det2', value=float(det2_val), min=0.0)
#	fit_params.add('det3', value=float(det3_val), min=0.0)

      if i == 1:
        xvals.append(float(z_c))
        if direction == 'high': x1.append(float(y_c))
        if direction == 'low': x2.append(float(y_c))
	fit_params = Parameters()
        fit_params.add('aCXB', value=float(aCXB_val))
        fit_params.add('det0', value=float(det0_val), min=0.0)
        #fit_params.add('det1', value=float(det1_val), min=0.0)
        fit_params.add('det2', value=float(det2_val), min=0.0)
        fit_params.add('det3', value=float(det3_val), min=0.0)

      if i == 2:
        xvals.append(float(w_c))
        if direction == 'high': x1.append(float(w_c))
        if direction == 'low': x2.append(float(w_c))
	fit_params = Parameters()
        fit_params.add('aCXB', value=float(aCXB_val))
        fit_params.add('det0', value=float(det0_val), min=0.0)
        fit_params.add('det1', value=float(det1_val), min=0.0)
        #fit_params.add('det2', value=float(det2_val), min=0.0)
        fit_params.add('det3', value=float(det3_val), min=0.0)

      if i == 3:
        xvals.append(float(z_c))
        if direction == 'high': x1.append(float(t_c))
        if direction == 'low': x2.append(float(t_c))
	fit_params = Parameters()
        fit_params.add('aCXB', value=float(aCXB_val))
        fit_params.add('det0', value=float(det0_val), min=0.0)
        fit_params.add('det1', value=float(det1_val), min=0.0)
        fit_params.add('det2', value=float(det2_val), min=0.0)
        #fit_params.add('det3', value=float(det3_val), min=0.0)

      yvals.append(m_cstat)
      # Need to find step values that will span from low cstat, to cstat+some_value in X many steps:
      # Option 1: solve for cstat+some_value for the parameter i'm adjusting and split that gap between inital value and final value by X steps
      # Option 2: Euler's method or 2nd order Runge-kutta method (both require derivatives and in this instance is not a good choice)
      print(m_cstat)
      if direction == 'high': step = np.abs(findStep(i)); fi = 'high_param_errors.txt'; y1.append(m_cstat)
      elif direction == 'low': step = np.abs(findStep(i))*-1; fi = 'low_param_errors.txt'; y2.append(m_cstat)
      while newcstat < mult_val+m_cstat:
        if i == 0: # index for which det we are in and find cstat values with those parameters
          print(step,'step')
          y_c += float(step)
          #print(y_c, fit_params['aCXB'].value, fit_params['det1'].value, fit_params['det2'].value, fit_params['det3'].value)
          resA = minimize(fn_cerror_det0,fit_params, method='nelder', tol=1e-15)
          #print(y_c, fit_params['aCXB'].value, fit_params['det1'].value, fit_params['det2'].value, fit_params['det3'].value)
          newcstat = fit_Cstat(fit_params['aCXB'].value,y_c,fit_params['det1'].value,fit_params['det2'].value,fit_params['det3'].value) 
          print(newcstat,'newcstat')
        elif i == 1:
          z_c += step
          resA = minimize(fn_cerror_det1,fit_params, method='nelder', tol=1e-15)
          newcstat = fit_Cstat(fit_params['aCXB'].value,fit_params['det0'].value,z_c,fit_params['det2'].value,fit_params['det3'].value)
        elif i == 2:
          w_c += step
          resA = minimize(fn_cerror_det2,fit_params, method='nelder', tol=1e-15)
          newcstat = fit_Cstat(fit_params['aCXB'].value,fit_params['det0'].value,fit_params['det1'].value,w_c,fit_params['det3'].value)
        elif i == 3:
          t_c += step
          resA = minimize(fn_cerror_det3,fit_params, method='nelder', tol=1e-15)
          newcstat = fit_Cstat(fit_params['aCXB'].value,fit_params['det0'].value,fit_params['det1'].value,fit_params['det2'].value,t_c)
    # see if new cstat value is lower than previous one
        if newcstat < m_cstat:
          m_cstat = None
          m_cstat = np.copy(newcstat)
          if len(newdetvals) > 0: newdetvals = []
          # run newacxb val
          if i == 0: 
            newdetvals.append(fit_params['aCXB'].value); newdet0 = float(np.copy(y_c)); newdetvals.append(newdet0); newdetvals.append(fit_params['det1'].value); newdetvals.append(fit_params['det2'].value); newdetvals.append(fit_params['det3'].value)
            xvals,yvals,x1,y1,x2,y2 =  [],[],[],[],[],[]
          if i == 1: 
            newdetvals.append(fit_params['aCXB'].value); newdetvals.append(fit_params['det0'].value); newdet1 = float(np.copy(z_c)); newdetvals.append(newdet1); newdetvals.append(fit_params['det2'].value); newdetvals.append(fit_params['det3'].value)
            xvals,yvals,x1,y1,x2,y2 =  [],[],[],[],[],[]
          if i == 2: 
            newdetvals.append(fit_params['aCXB'].value); newdetvals.append(fit_params['det0'].value); newdetvals.append(fit_params['det1'].value); newdet2 = float(np.copy(w_c)); newdetvals.append(newdet2); newdetvals.append(fit_params['det3'].value)
            xvals,yvals,x1,y1,x2,y2 =  [],[],[],[],[],[]
          if i == 3: 
            newdetvals.append(fit_params['aCXB'].value); newdetvals.append(fit_params['det0'].value); newdetvals.append(fit_params['det1'].value); newdetvals.append(fit_params['det2'].value); newdet3 = float(np.copy(t_c)); newdetvals.append(newdet3)
            xvals,yvals,x1,y1,x2,y2 =  [],[],[],[],[],[]
          blah = open('high_param_errors.txt','w+')
          blah.close()
          blah = open('low_param_errors.txt','w+')
          blah.close()
# place reset values here:  
          params['aCXB'].set(value = newdetvals[0])
          params['det0'].set(value = newdetvals[1], min=0.0)
          params['det1'].set(value = newdetvals[2], min=0.0)
          params['det2'].set(value = newdetvals[3], min=0.0)
          params['det3'].set(value = newdetvals[4], min=0.0)	
	  prev_mstat = np.copy(m_cstat)
          m_cstat = runNewaCXB()
          newcstat = np.copy(m_cstat)
 	  #print(prev_mstat,m_cstat,'stats')
# Need to reset params in case there is chance newaCXB function doesn't find a new low value than the one found here:
          return False  # this was continue
        if i == 0: 
          xvals.append(float(y_c)) 
          if direction == 'high': x1.append(float(y_c))
          if direction == 'low': x2.append(float(y_c))
        if i == 1: 
          xvals.append(float(z_c)) 
          if direction == 'high': x1.append(float(z_c))
          if direction == 'low': x2.append(float(z_c))
        if i == 2: 
          xvals.append(float(w_c)) 
          if direction == 'high': x1.append(float(w_c))
          if direction == 'low': x2.append(float(w_c))
        if i == 3: 
          xvals.append(float(t_c))
          if direction == 'high': x1.append(float(t_c))
          if direction == 'low': x2.append(float(t_c))
        yvals.append(float(newcstat))
        if direction == 'high': y1.append(float(newcstat))
        if direction == 'low': y2.append(float(newcstat))
        # write parameters to temp file to be combined later
        with open(fi,'a+') as O:
          if i == 0:
            O.write(str(fit_params['aCXB'].value)+' '+str(y_c)+' '+str(fit_params['det1'].value)+' '+str(fit_params['det2'].value)+' '+str(fit_params['det3'].value)+'\n')
          if i == 1:
            O.write(str(fit_params['aCXB'].value)+' '+str(fit_params['det0'].value)+' '+str(z_c)+' '+str(fit_params['det2'].value)+' '+str(fit_params['det3'].value)+'\n')
          if i == 2:
            O.write(str(fit_params['aCXB'].value)+' '+str(fit_params['det0'].value)+' '+str(fit_params['det1'].value)+' '+str(w_c)+' '+str(fit_params['det3'].value)+'\n')
          if i == 3:
            O.write(str(fit_params['aCXB'].value)+' '+str(fit_params['det0'].value)+' '+str(fit_params['det1'].value)+' '+str(fit_params['det2'].value)+' '+str(t_c)+'\n')
      
      
      #reset parameters for next run:
      if i == 0:    
        setParameters('aCXB', params['aCXB'].value, 'det1', params['det1'].value, 'det2', params['det2'].value, 'det3', params['det3'].value)
        y_c = np.copy(params['det0'].value)
      if i == 1:
        setParameters('aCXB', params['aCXB'].value, 'det0', params['det0'].value, 'det2', params['det2'].value, 'det3', params['det3'].value)
        z_c = np.copy(params['det1'].value)
      if i == 2:
        setParameters('aCXB', params['aCXB'].value, 'det0', params['det0'].value, 'det1', params['det1'].value, 'det3', params['det3'].value)
        w_c = np.copy(params['det2'].value)
      if i == 3:
        setParameters('aCXB', params['aCXB'].value, 'det0', params['det0'].value, 'det1', params['det1'].value, 'det2', params['det2'].value)
        t_c = np.copy(params['det3'].value)
      newcstat = np.copy(m_cstat)

    # write out the norm values found:
    fi_high = 'high_param_errors.txt'
    fi_low = 'low_param_errors.txt'
    fwrite = detp+'_'+detarray[i]+'_params_errors.txt'

    if not os.path.isfile(fi_high): open(fi_high,'a').close()
    with open(fi_high) as fh:
        lines = fh.readlines()
        for line in lines:
          with open(fwrite,'w') as fw:
                fw.write(str(line)+'\n')
    if not os.path.isfile(fi_low): open(fi_low,'a').close()
    with open(fi_low) as fl:
        lines = fl.readlines()
        for line in lines:
          with open(fwrite,'a') as fw:
                fw.write(str(line)+'\n')

    os.system('rm '+fi_high+' '+fi_low)
#need to write out values to parameter file
    # find the one sig errors 
    p,c = curve_fit(parabola,xvals,yvals)
    up = True; down = True     #
    if len(x1) > 2:        #
      pars1, cov1 = opt.curve_fit(parabola,x1,y1,p)
      sigup = solu(m_cstat,*pars1)
      up = False      #
    if len(x2) > 2:          #
      pars2, cov2 = opt.curve_fit(parabola,x2,y2,p)
      sigdown = sold(m_cstat,*pars2)
      down = False   #

    if up:     #
	sigup = solu(m_cstat, *p)    #
    if down:   #
        sigdown = sold(m_cstat, *p)  #
    
#    sigup = solu(m_cstat,*pars1)
#    sigdown = sold(m_cstat,*pars2)

    # write the energy and errors
    with open(detp+'_'+detarray[i]+"_params.txt",'a+') as O:
      if i == 0:
        O.write(lowe+' '+highe+' '+str(params['det0'].value)+' '+str(sigup)+' '+str(sigdown)+' '+str(exp)+'\n')
      if i == 1:
        O.write(lowe+' '+highe+' '+str(params['det1'].value)+' '+str(sigup)+' '+str(sigdown)+' '+str(exp)+'\n')
      if i == 2:
        O.write(lowe+' '+highe+' '+str(params['det2'].value)+' '+str(sigup)+' '+str(sigdown)+' '+str(exp)+'\n')
      if i == 3:
        O.write(lowe+' '+highe+' '+str(params['det3'].value)+' '+str(sigup)+' '+str(sigdown)+' '+str(exp)+'\n')
     
  #reset variables:
    xvals, x1, x2 = [],[],[]
    yvals, y1, y2 = [],[],[]
    newcstat = np.copy(m_cstat)
    y_c = params['det0'].value; z_c = params['det1'].value; w_c = params['det2'].value; t_c = params['det3'].value; x_c = params['aCXB'].value
  return True
#'''
####  MAIN   ##########

def clearFile(fi):
  with open(fi, "r") as f:
    lines = f.readlines()
  with open(fi, "w") as f:
    for line in lines:
      if lowe in line:
        if highe in line:
          continue
      f.write(line)


while True:
  if doCstatRun():
    break
  # if I get a return of False, I need to clear any existing files that may have been written to:
  else:
  #if not doCstatRun(): # == False:
    for de in detarray:
      param_fi = detp+'_'+de+'_params.txt'
      param_ferr = detp+'_'+de+'_params_errors.txt'
      if os.path.isfile(param_fi):
      	clearFile(detp+'_'+de+'_params.txt')
      if os.path.isfile(param_ferr):
      	clearFile(detp+'_'+de+'_params_errors.txt')

print('finished')
















