#!/usr/bin/python
# syntax: ./randnumgen.py
""" x is a value given for the maximum t in a sin function. The return is values of t and a semi-random
    distribution of y around sin(x)"""

import sys
import os
import string
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


fname = 'tempidldata.txt'
datamatrix = np.loadtxt(fname)
xcurve = datamatrix[:,0]
ycurve = datamatrix[:,1]
print(datamatrix)
print(xcurve)
print(ycurve)
A = np.pi*2
x = xcurve
y = ycurve
plt.plot(xcurve,ycurve,'bo')

##############################################
"""This is the funct is driving me nuts"""
def funcy(x,A):
    
    t = x
    
    y = 0.2*np.sin(A*t+1.8)+0.5*np.random.normal(-1,1,len(t))
   
    return y,t
##############################################
def func(t,a,b,c,d):
    return a*np.sin(b*t+c)+d
##############################################



##############################################
def main():


    y, t = funcy(x,A)
    
    stddev = np.std(y)
    mean = np.median(y)
    print(stddev)
    print(mean)
    
    Af = np.array([0,0,0,0])
    Af[0] = stddev
    Af[1] = A
    Af[2] = 3
    Af[3] = mean

    popt, pcov = curve_fit(func,t,y,p0=Af)
    print(popt)
    
    yfit = 0.2*np.sin(A*t+1.5)
    yfit2 = func(t,*popt)
    
    #plt.plot(t,yfit,'g-')
    #plt.plot(t,yfit2,'r-')
    
    plt.show()

main()
