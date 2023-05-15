import os, sys, numpy as np, matplotlib.pyplot as plt
from astropy.io import fits
from lmfit import Parameters

def testFunc():
    global m_cstat, mode
    i = 0
    while i < 15:
      print(i)
      i += 1
      if i == 5: 
        mode = 'quit' 
        break
      m_cstat = param['one'].value
      param['one'].value += 0.25
      print('still running') 
param = Parameters()
param.add('one',value = 0.0)
param.add('two',value = 0.0)
param.add('three',value = 0.0)

testlist = []
m_cstat = 4.0
br = 'no'
mode = 'pass'

while m_cstat < 10:
  testFunc()
  if mode == 'quit': 
    print('it quit') 
    break
  if mode == 'pass': print('it finished?')


