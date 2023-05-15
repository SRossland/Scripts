#!/bin/python

import numpy as np, matplotlib.pyplot as plt
import os, sys

Adata = sys.argv[1]
Bdata = sys.argv[2]
saveplot = sys.argv[3]

if saveplot.lower() not in ['save','plot']: 
	print('You need to specify if you are saving a figure or plotting!')
	sys.exit()

twoplot = False
mcmc = False
if len(sys.argv) > 4:
  if sys.argv[4] == 'mcmc':
    Adata_mc = sys.argv[5]
    Bdata_mc = sys.argv[6]
    Adatmc = np.genfromtxt(Adata_mc)
    Bdatmc = np.genfromtxt(Bdata_mc)
    mcmc = True
  else:
    Adata_opt = sys.argv[4]
    Bdata_opt = sys.argv[5]
    twoplot = True  

Adat = np.genfromtxt(Adata)
Bdat = np.genfromtxt(Bdata)

binwidthA = (Adat[:,1]-0.04)-Adat[:,0]
binwidthB = (Bdat[:,1]-0.04)-Bdat[:,0]

if len(set(Adat[:,14])) == 1:
  expA = np.copy(Adat[0,14])/1e6
else:
  expA = np.copy(Adat[:,14])/1e6

if len(set(Bdat[:,14])) == 1:
  expB = np.copy(Bdat[0,14])/1e6
else:
  expB = np.copy(Bdat[:,14])/1e6

#Av0 = np.copy(Adat[:,2]);Av1 = np.copy(Adat[:,5]); Av2 = np.copy(Adat[:,8]); Av3 = np.copy(Adat[:,11])
#Bv0 = np.copy(Bdat[:,2]);Bv1 = np.copy(Bdat[:,5]); Bv2 = np.copy(Bdat[:,8]); Bv3 = np.copy(Bdat[:,11])

#Ae0u = np.copy(Adat[:,3]); Ae1u = np.copy(Adat[:,6]); Ae2u = np.copy(Adat[:,9]); Ae3u = np.copy(Adat[:,12])
#Be0u = np.copy(Bdat[:,3]); Be1u = np.copy(Bdat[:,6]); Be2u = np.copy(Bdat[:,9]); Be3u = np.copy(Bdat[:,12])
#Ae0l = np.copy(Adat[:,4]); Ae1l = np.copy(Adat[:,7]); Ae2l = np.copy(Adat[:,10]); Ae3l = np.copy(Adat[:,13])
#Be0l = np.copy(Bdat[:,4]); Be1l = np.copy(Bdat[:,7]); Be2l = np.copy(Bdat[:,10]); Be3l = np.copy(Bdat[:,13])

#Av0 /= binwidthA; Av1 /= binwidthA; Av2 /= binwidthA; Av3 /= binwidthA
#Bv0 /= binwidthB; Bv1 /= binwidthB; Bv2 /= binwidthB; Bv3 /= binwidthB
#Av0 /= expA; Av1 /= expA; Av2 /= expA; Av3 /= expA
#Bv0 /= expB; Bv1 /= expB; Bv2 /= expB; Bv3 /= expB

#Ae0u /= binwidthA; Ae1u /= binwidthA; Ae2u /= binwidthA; Ae3u /= binwidthA
#Ae0u /= expA; Ae1u /= expA; Ae2u /= expA; Ae3u /= expA
#Ae0l /= binwidthA; Ae1l /= binwidthA; Ae2l /= binwidthA; Ae3l /= binwidthA
#Ae0l /= expA; Ae1l /= expA; Ae2l /= expA; Ae3l /= expA

#Be0u /= binwidthB; Be1u /= binwidthB; Be2u /= binwidthB; Be3u /= binwidthB
#Be0u /= expB; Be1u /= expB; Be2u /= expB; Be3u /= expB
#Be0l /= binwidthB; Be1l /= binwidthB; Be2l /= binwidthB; Be3l /= binwidthB
#Be0l /= expB; Be1l /= expB; Be2l /= expB; Be3l /= expB

# change the limits to a difference instead of a value:

#Ae0u = Ae0u - Av0; Ae1u = Ae1u - Av1; Ae2u = Ae2u - Av2; Ae3u = Ae3u - Av3
#Be0u = Be0u - Bv0; Be1u = Be1u - Bv1; Be2u = Be2u - Bv2; Be3u = Be3u - Bv3
#Ae0l = (Ae0l - Av0) * -1; Ae1l = (Ae1l - Av1) * -1; Ae2l = (Ae2l - Av2) * -1; Ae3l = (Ae3l - Av3) * -1
#Be0l = (Be0l - Bv0) * -1; Be1l = (Be1l - Bv1) * -1; Be2l = (Be2l - Bv2) * -1; Be3l = (Be3l - Bv3) * -1

# Xvalues

xval = Adat[:,0]+((Adat[:,1]-Adat[:,0])/2)

xerrl = xval-Adat[:,0]
xerru = Adat[:,1]-xval 

#################################
# Definitions to generally get the data needed so I can do other things:

def getDataValues(dat, exp, binwidth):
	v0 = np.copy(dat[:,2]);v1 = np.copy(dat[:,5]); v2 = np.copy(dat[:,8]); v3 = np.copy(dat[:,11])
	e0u = np.copy(dat[:,3]); e1u = np.copy(dat[:,6]); e2u = np.copy(dat[:,9]); e3u = np.copy(dat[:,12])
	e0l = np.copy(dat[:,4]); e1l = np.copy(dat[:,7]); e2l = np.copy(dat[:,10]); e3l = np.copy(dat[:,13])
	v0 /= binwidth; v1 /= binwidth; v2 /= binwidth; v3 /= binwidth
	v0 /= exp; v1 /= exp; v2 /= exp; v3 /= exp
	e0u /= binwidth; e1u /= binwidth; e2u /= binwidth; e3u /= binwidth
	e0u /= exp; e1u /= exp; e2u /= exp; e3u /= exp
	e0l /= binwidth; e1l /= binwidth; e2l /= binwidth; e3l /= binwidth
	e0l /= exp; e1l /= exp; e2l /= exp; e3l /= exp
	e0u = e0u - v0; e1u = e1u - v1; e2u = e2u - v2; e3u = e3u - v3
	e0l = (e0l - v0) * -1; e1l = (e1l - v1) * -1; e2l = (e2l - v2) * -1; e3l = (e3l - v3) * -1
	return v0, e0u, e0l, v1, e1u, e1l, v2, e2u, e2l, v3, e3u, e3l 
################################

Av0, Ae0u, Ae0l, Av1, Ae1u, Ae1l, Av2, Ae2u, Ae2l, Av3, Ae3u, Ae3l = getDataValues(Adat, expA, binwidthA)
Bv0, Be0u, Be0l, Bv1, Be1u, Be1l, Bv2, Be2u, Be2l, Bv3, Be3u, Be3l = getDataValues(Bdat, expB, binwidthB)
if mcmc:
	Av0c, Ae0uc, Ae0lc, Av1c, Ae1uc, Ae1lc, Av2c, Ae2uc, Ae2lc, Av3c, Ae3uc, Ae3lc = getDataValues(Adatmc, expA, binwidthA)
	Bv0c, Be0uc, Be0lc, Bv1c, Be1uc, Be1lc, Bv2c, Be2uc, Be2lc, Bv3c, Be3uc, Be3lc = getDataValues(Bdatmc, expB, binwidthB)




# get optional data



def getValues(fi):
  f = open(fi,'r+')
  A,B,C,D,E = [],[],[],[],[]
  for line in f.readlines():
    a,b,c,d,e,g = line.split()
    A.append(float(a))
    B.append(float(b))
    C.append(float(c))
    D.append(float(d))
    E.append(float(e))
  G = float(g)
  f.close()
  return A,B,C,D,E,G

def binWidth(vals):
  binwid = []
  wid = 0
  for i,val in enumerate(vals):
    if i == len(vals)-1:
      binwid.append(wid)
    else:
      wid = round((vals[i+1]-0.04)-val,2)
      binwid.append(wid)
  return binwid

if twoplot:
	det0A_01,det1A_01,det2A_01,det3A_01,Ax_01,expA01 = getValues(Adata_opt)
	det0B_01,det1B_01,det2B_01,det3B_01,Bx_01,expB01 = getValues(Bdata_opt)

	expA01 /= 1e6
	expB01 /= 1e6
	A01vals = sorted(zip(Ax_01,det0A_01,det1A_01,det2A_01,det3A_01), key = lambda x: float(x[0]))
	B01vals = sorted(zip(Bx_01,det0B_01,det1B_01,det2B_01,det3B_01), key = lambda x: float(x[0]))
	#A01vals.sort()
	#B01vals.sort()
	Ax01,det0A01,det1A01,det2A01,det3A01 = zip(*A01vals)
	Bx01,det0B01,det1B01,det2B01,det3B01 = zip(*B01vals)
	bin_width_A_01 = binWidth(Ax01)
	bin_width_B_01 = binWidth(Bx01)
	bwA01 = np.asarray(bin_width_A_01)
	bwB01 = np.asarray(bin_width_B_01)
	d0A01 = np.asarray(det0A01)
	d1A01 = np.asarray(det1A01)
	d2A01 = np.asarray(det2A01)
	d3A01 = np.asarray(det3A01)
	Ax01 = np.asarray(Ax01)
	d0B01 = np.asarray(det0B01)
	d1B01 = np.asarray(det1B01)
	d2B01 = np.asarray(det2B01)
	d3B01 = np.asarray(det3B01)
	Bx01 = np.asarray(Bx01)
	d0A01 = d0A01/bwA01
	d1A01 = d1A01/bwA01
	d2A01 = d2A01/bwA01
	d3A01 = d3A01/bwA01
	d0B01 = d0B01/bwB01
	d1B01 = d1B01/bwB01
	d2B01 = d2B01/bwB01
	d3B01 = d3B01/bwB01
	d0A01 /= expA01
	d1A01 /= expA01
	d2A01 /= expA01
	d3A01 /= expA01
	d0B01 /= expB01
	d1B01 /= expB01
	d2B01 /= expB01
	d3B01 /= expB01

# plot

fig, axs = plt.subplots(2,2)
axs[0,0].errorbar(xval, Av0, xerr=(xerrl,xerru), yerr=(Ae0l,Ae0u),fmt='or',markersize=3,label='After')
if twoplot: axs[0,0].scatter(Ax01,d0A01,c='blue',label='Before', marker='x')
if mcmc: axs[0,0].errorbar(xval, Av0c, xerr=(xerrl,xerru), yerr=(Ae0lc,Ae0uc),fmt='xb',markersize=3,label='mcmc')
axs[0,0].set_title('Detector 0')
axs[0,0].legend(loc=1)
axs[0,1].errorbar(xval, Av1, xerr=(xerrl,xerru), yerr=(Ae1l,Ae1u),fmt='or',markersize=3,label='After')
if twoplot: axs[0,1].scatter(Ax01,d1A01,c='blue',label='Before', marker='x')
if mcmc: axs[0,1].errorbar(xval, Av1c, xerr=(xerrl,xerru), yerr=(Ae1lc,Ae1uc),fmt='xb',markersize=3,label='mcmc')
axs[0,1].set_title('Detector 1')
axs[0,1].legend(loc=1)
axs[1,0].errorbar(xval, Av2, xerr=(xerrl,xerru), yerr=(Ae2l,Ae2u),fmt='or',markersize=3,label='After')
if twoplot: axs[1,0].scatter(Ax01,d2A01,c='blue',label='Before', marker='x')
if mcmc: axs[1,0].errorbar(xval, Av2c, xerr=(xerrl,xerru), yerr=(Ae2lc,Ae2uc),fmt='xb',markersize=3,label='mcmc')
axs[1,0].set_title('Detector 2')
axs[1,0].legend(loc=1)
axs[1,1].errorbar(xval, Av3, xerr=(xerrl,xerru), yerr=(Ae3l,Ae3u),fmt='or',markersize=3,label='After')
if twoplot: axs[1,1].scatter(Ax01,d3A01,c='blue',label='Before', marker='x')
if mcmc: axs[1,1].errorbar(xval, Av3c, xerr=(xerrl,xerru), yerr=(Ae3lc,Ae3uc),fmt='xb',markersize=3,label='mcmc')
axs[1,1].set_title('Detector 3')
axs[1,1].legend(loc=1)
if twoplot: fig.suptitle('Telescope A det norms before and after error finding')
if mcmc: fig.suptitle('Telescope A det norms: Least squares vs MCMC')
else: fig.suptitle('Telescope A det norms after error finding')
for ax in axs.flat:
  ax.set(xlabel='Energy keV', ylabel='Norm val (cts/Ms/binwidth)')
for ax in axs.flat:
  ax.label_outer()

if saveplot.lower() == 'save':
  plt.savefig('A_deterror_normvalues.png')
if saveplot.lower() == 'plot': 
  plt.show()

fig, axs = plt.subplots(2,2)
axs[0,0].errorbar(xval, Bv0, xerr=(xerrl,xerru), yerr=(Be0l,Be0u),fmt='or',markersize=3,label='After')
if twoplot: axs[0,0].scatter(Bx01,d0B01,c='blue',label='Before', marker='x')
if mcmc: axs[0,0].errorbar(xval, Bv0c, xerr=(xerrl,xerru), yerr=(Be0lc,Be0uc),fmt='xb',markersize=3,label='mcmc')
axs[0,0].set_title('Detector 0')
axs[0,0].legend(loc=1)
axs[0,1].errorbar(xval, Bv1, xerr=(xerrl,xerru), yerr=(Be1l,Be1u),fmt='or',markersize=3,label='After')
if twoplot: axs[0,1].scatter(Bx01,d1B01,c='blue',label='Before', marker='x')
if mcmc: axs[0,1].errorbar(xval, Bv1c, xerr=(xerrl,xerru), yerr=(Be1lc,Be1uc),fmt='xb',markersize=3,label='mcmc')
axs[0,1].set_title('Detector 1')
axs[0,1].legend(loc=1)
axs[1,0].errorbar(xval, Bv2, xerr=(xerrl,xerru), yerr=(Be2l,Be2u),fmt='or',markersize=3,label='After')
if twoplot: axs[1,0].scatter(Bx01,d2B01,c='blue',label='Before', marker='x')
if mcmc: axs[1,0].errorbar(xval, Bv2c, xerr=(xerrl,xerru), yerr=(Be2lc,Be2uc),fmt='xb',markersize=3,label='mcmc')
axs[1,0].set_title('Detector 2')
axs[1,0].legend(loc=1)
axs[1,1].errorbar(xval, Bv3, xerr=(xerrl,xerru), yerr=(Be3l,Be3u),fmt='or',markersize=3,label='After')
if twoplot: axs[1,1].scatter(Bx01,d3B01,c='blue',label='Before', marker='x')
if mcmc: axs[1,1].errorbar(xval, Bv3c, xerr=(xerrl,xerru), yerr=(Be3lc,Be3uc),fmt='xb',markersize=3,label='mcmc')
axs[1,1].set_title('Detector 3')
axs[1,1].legend(loc=1)
if twoplot: fig.suptitle('Telescope B det norms before and after error finding')
if mcmc: fig.suptitle('Telescope B det norms: Least squares vs MCMC')
else: fig.suptitle('Telescope B det norms after error finding')
for ax in axs.flat:
  ax.set(xlabel='Energy keV', ylabel='Norm val (cts/Ms/binwidth)')
for ax in axs.flat:
  ax.label_outer()

if saveplot.lower() == 'save':
  plt.savefig('B_deterror_normvalues.png')
if saveplot.lower() == 'plot':
  plt.show()

