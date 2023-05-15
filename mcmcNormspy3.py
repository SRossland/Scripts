#! /uufs/astro.utah.edu/common/home/u1019304/VENV3.7.9/bin/python3.7

import emcee, os, sys
import numpy as np, matplotlib.pyplot as plt
import corner, time
from astropy.io import fits

start_time = time.time()
############# Load data

data_file = sys.argv[1]
bins_file = sys.argv[2]
det_map_file = sys.argv[3]
exp_map_file = sys.argv[4]
Aij_file = sys.argv[5]
norm_file = sys.argv[6]

data_split = data_file.split('/')[-1]
lowe = data_split.split('_')[3]
hig = data_split.split('_')[4]
highe, crap = hig.split('keV.fits')

detp = data_split.split('_')[0][4]

############


########### Definitions

def lnprior(p):
	# This will be the parameters for the fitting, we can set limits 
	# any value less than 0.0 and each norm value will have it's own max
    x_c, y_c, z_c, w_c, t_c = p
    if (x_c <= 0.0 or x_c > 1.0 or y_c <= 0.0 or z_c <= 0.0 or w_c <= 0.0 or t_c <= 0.0):
        return -np.inf
    return 0

def lnlike(p,A,det0,det1,det2,det3,count):
	# This is the log liklihood from the previous method
    x_c, y_c, z_c, w_c, t_c = p
    model = ((x_c*A)+(y_c*det0)+(z_c*det1)+(w_c*det2)+(t_c*det3))
    numerA = model - count + count*(np.log(count) - np.log(model))
    lp = -2*np.sum(numerA)
    return lp

def lnprob(p, A, det0, det1, det2, det3, count):
    lp = lnprior(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(p, A, det0, det1, det2, det3, count)

# A change I implemented was to switch the lnlike to lp = (-) and the lnporb to lp (+) 
# This change may help find the value I am looking for 

def getData(fi):
    with fits.open(fi) as hdul:
        var = hdul[0].data
    return var
###########


mask_file = '/uufs/astro.utah.edu/common/home/u1019304/temp/fullmask'+detp+'_final.fits'
mask = getData(mask_file)

data = getData(data_file)
bins = getData(bins_file)
det_map = getData(det_map_file)
exp_map = getData(exp_map_file)
gradA = getData(Aij_file)

det_map = det_map*mask
exp_map = exp_map*mask
data = data*mask

Aij = gradA[int(len(gradA)/2) - 180:int(len(gradA)/2) + 180, int(len(gradA)/2) - 180:int(len(gradA)/2) + 180]*mask

Aij *= 0.0001

a0 = np.zeros((360,360)); a1 = np.zeros((360,360)); a2 = np.zeros((360,360)); a3 = np.zeros((360,360))
a0[det_map == 0] = 1
a1[det_map == 1] = 1
a2[det_map == 2] = 1
a3[det_map == 3] = 1

with fits.open(data_file) as hdul:
    exp = hdul[0].header['EXPOSURE']

expmask = np.copy(exp_map)
expmask /= exp

Aij *= expmask; a0 *= expmask; a1 *= expmask; a2*= expmask; a3*= expmask

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

f = open(norm_file, 'r+')
for line in f.readlines():
    if 'aCXB norm' in line:
        aCXB_val = float(line.split()[3])
    if 'det0 norm' in line:
        det0_val = float(line.split()[3])
    if 'det1 norm' in line:
        det1_val = float(line.split()[3])
    if 'det2 norm' in line:
        det2_val = float(line.split()[3])
    if 'det3 norm' in line:
        det3_val  = float(line.split()[3])
f.close()


inital = np.array([aCXB_val, det0_val, det1_val, det2_val, det3_val])

########### Main

# Need to use the values from the previous measurement, so load those

Nwalker, Ndim = 100,5

p0 = [np.array(inital)+1.e-6*np.random.randn(Ndim) for i in range(Nwalker)]

sampler = emcee.EnsembleSampler(Nwalker,Ndim,lnprob,args=(A_bins, det0_bins, det1_bins, det2_bins, det3_bins, count_bins))

#pos,prob,state = sampler.run_mcmc(p0, 500)
########################## TESTING
# can look at the walker trace here if you want

#trace
#fig, axs = plt.subplots(5)
#axs[0].plot(sampler.chain[:,:,0].T, '-', color='k', alpha=0.3)
#axs[0].axhline(aCXB_val, color='blue')
#axs[1].plot(sampler.chain[:,:,1].T, '-', color='k', alpha=0.3)
#axs[1].axhline(det0_val, color='blue')
#axs[2].plot(sampler.chain[:,:,2].T, '-', color='k', alpha=0.3)
#axs[2].axhline(det1_val, color='blue')
#axs[3].plot(sampler.chain[:,:,3].T, '-', color='k', alpha=0.3)
#axs[3].axhline(det2_val, color='blue')
#axs[4].plot(sampler.chain[:,:,4].T, '-', color='k', alpha=0.3)
#axs[4].axhline(det3_val, color='blue')
#plt.show()

# Values
#a,d0,d1,d2,d3 = np.median(sampler.flatchain,axis=0)
#print(a,d0,d1,d2,d3)
# Errors
#flat_samples = sampler.get_chain(discard = 100, thin=15, flat=True)
#print(flat_samples.shape)
#for i in range(Ndim):
#    mcmc = np.percentile(flat_samples[:,i], [10,50,90])
#    q = np.diff(mcmc)
#    print(mcmc[1],mcmc[1]+q[1],mcmc[1]-q[0])
#sys.exit()
##################
sampler.reset()


pos, prob, state = sampler.run_mcmc(p0, 10000)

flat_samples = sampler.get_chain(discard = 80, thin=10, flat=True)
printfiles = [detp+'_params_mcmc.txt', detp+'_det0_params_mcmc.txt',detp+'_det1_params_mcmc.txt',detp+'_det2_params_mcmc.txt',detp+'_det3_params_mcmc.txt']
for i in range(Ndim):
    mcmc = np.percentile(flat_samples[:,i], [16,50,84])
    q = np.diff(mcmc)
    #print(mcmc[1],mcmc[1]+q[1],mcmc[1]-q[0])
    with open(printfiles[i],'a+') as O:
        O.write(lowe+' '+highe+' '+str(mcmc[1])+' '+str(mcmc[1]+q[1])+' '+str(mcmc[1]-q[0])+'\n')

#print(np.median(sampler.flatchain, axis=0))
#print(pos, prob, state)
################################   

fig, axs = plt.subplots(5)
axs[0].plot(sampler.chain[:,:,0].T, '-', color='k', alpha=0.3)
axs[0].axhline(aCXB_val, color='blue')
axs[1].plot(sampler.chain[:,:,1].T, '-', color='k', alpha=0.3)
axs[1].axhline(det0_val, color='blue')
axs[2].plot(sampler.chain[:,:,2].T, '-', color='k', alpha=0.3)
axs[2].axhline(det1_val, color='blue')
axs[3].plot(sampler.chain[:,:,3].T, '-', color='k', alpha=0.3)
axs[3].axhline(det2_val, color='blue')
axs[4].plot(sampler.chain[:,:,4].T, '-', color='k', alpha=0.3)
axs[4].axhline(det3_val, color='blue')
fig.savefig('chains_'+detp+'_'+lowe+'to'+highe+'.png')

#print(sampler.flatchain)
fig = corner.corner(flat_samples, show_titles=True, labels=['aCXB','DET0','DET1','DET2','DET3'],truths = [aCXB_val, det0_val, det1_val, det2_val, det3_val], quantiles=[0.16,0.5,0.84])
fig.savefig('corner_norms_'+detp+'_'+lowe+'to'+highe+'.png')

plt.close('all')
#print('---------%s seconds---------' % (time.time() - start_time))









# Fin
