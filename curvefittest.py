#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./curvefittest.py A
""" xcurve is the x data
    ycurve is the y data
    A is a param array assumed to be of length 4"""

import sys
import os
import string
import numpy as np
from scipy.optimize import curve_fit

A = sys.argv[1]
A = float(A)
fname = 'tempidldata.txt'

obs = 'obs.txt'
obsi = np.loadtxt(obs)
obsid = str(int(obsi))

datamatrix = np.loadtxt(fname)
xcurve = datamatrix[:,0]
ycurve = datamatrix[:,1]
err = datamatrix[:,2]

##############################################
"""This is the function y(x) = a*sin(b*t+c)+d that will be used to fit"""
def funcsine(t,a,b,c,d):
    return a*np.sin(b*t+c)+d
##############################################
"""f(x) = Ax + B"""
def funcslope(x,A,B):
    return A*x + B
##############################################
"""Gaussian (This one will never work, but is used as a sanity check)"""
def gaussian(x,a,b,c):
    return a * (np.exp(b-x)**2)/(c**2)
##############################################
"""Future models, cause hey! who doesn't love a good model?!?"""

##############################################
"""Chi-Squared calculation"""
def redchi(obs,mod,err,dof):
    chi = np.sum(((obs-mod)**2)/err**2)
    rchi = chi/dof
    return rchi
##############################################

def main():
    
    stddev = np.std(ycurve)
    mean = np.median(ycurve)

    Af = np.array([0.,0.,0.,0.])
    Af[0] = stddev
    Af[1] = A
    Af[2] = 1
    Af[3] = mean
    
    """Curve fit corner"""
    poptsine, pcovsine = curve_fit(funcsine,xcurve,ycurve,p0=Af,maxfev=10000)
    poptslope, pcovslope = curve_fit(funcslope,xcurve,ycurve,maxfev=5000)
    #poptgaus, pcovgaus = curve_fit(gaussian,xcurve,ycurve,maxfev=10000)
    
    """future modles just need to be added above in calculation and computed/listed in ys"""
    yfitsine = funcsine(xcurve,*poptsine)
    yfitslope = funcslope(xcurve,*poptslope)
    #yfitgaus = gaussian(xcurve,*poptgaus)
    yfitline = np.array([mean]*len(ycurve))
    ys = [yfitsine,yfitslope,yfitline]
    
    """add a zero for each model and rc with dof"""
    rc = np.array([0.,0.,0.])
    rc[0] = redchi(ycurve,yfitsine,err,(len(xcurve)-len(poptsine)))
    rc[1] = redchi(ycurve,yfitslope,err,(len(xcurve)-len(poptslope)))
    #rc[2] = redchi(ycurve,yfitgaus,err,len(poptgaus))
    rc[2] = redchi(ycurve,yfitline,err,len(xcurve)-1)
    rc2 = np.abs(1-rc)
    val = np.argmin(rc2)
    
    dir = os.getcwd()
    with open(dir+'/OBS'+'/'+obsid+'/Proc_log.txt','a+') as OTF:
        OTF.write("\n"+str(rc)+"\n    *sin: "+str(poptsine)+"\n    slope: "+str(poptslope) \
                  +"\n    median: "+str(mean))
              
    yfit = np.copy(ys[val])
    os.system('rm -f tempidldata.txt')
    np.savetxt('tempcft.txt',yfit)
main()
