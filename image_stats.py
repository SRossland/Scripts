#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# Syntax: image_stat.py path/to/obs obslist mode time le he det homedir

# This is a wrapper of sorts to create images in an energy bandwidth, stack those images together, 
# get the statistics of those images and a model assoicated with, and repeat the process.

# A new import will be pexpect, it is an interactive module that should allow me to interact with
# previous scripts without the need to rewrite them.

import pexpect, os, sys, numpy as np

# The two other programs that I will be using in this are:
#  (1) image_stack.py 
#  (2) count_stat_new.py

#  to use them they have a few arguments unique to each:
#  (1) image_stack.py obsdir obsid mode el_low el_high det(optional)
#  (2) count_stat_new.py det imagefile(created from (1))

obsdir = sys.argv[1]
obslist = sys.argv[2]
mode = sys.argv[3]
sep = sys.argv[4]
lowe = sys.argv[5]
highe = sys.argv[6]
det = sys.argv[7]
homedir = sys.argv[8]
energygap = sys.argv[9]
skip = sys.argv[10]
file_name_tag = sys.argv[11]
offx = sys.argv[12]
offy = sys.argv[13]

if det == "BOTH": 
  det1 = 'BOTH'
  det2 = ['A','B']
else:
  det1 = det
  det2 = det

def Energy(x):
  return x*0.04+1.6

def PI(x):
  return (x-1.6)/0.04

# Create the energy bands:
lowenergy, highenergy = [],[]
i = float(lowe)
j = 1
while i < float(highe):
  j += 1
  lowenergy.append(i)
  i *= float(energygap)
  i = round(Energy(int(PI(i))),2)
  highenergy.append(i)

highe = str(highenergy[-1])

# The file system is XRAY/NuSTAR/CXB/01/full/3_30_9steps/fits_file

end_file_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/CXB/'+mode
#full_end = end_file_dir+'/full/'+lowe+'_'+highe+'_'+str(j)+'steps'	#These lines are commented out to force only nosun files
#sun_end = end_file_dir+'/sun/'+lowe+'_'+highe+'_'+str(j)+'steps'	#
#nosun_end = end_file_dir+'/nosun/'+lowe+'_'+highe+'_'+str(j)+'steps'	#
#fi = [full_end,sun_end,nosun_end]					#
#fi = end_file_dir + '/nosun/'+lowe+'_'+highe+'_'+str(j)+'steps_'+file_name_tag

######## loop over the energies ##########
for m in range(len(lowenergy)):
  os.system("image_stack.py "+obsdir+" "+obslist+" "+mode+" "+str(lowenergy[m])+" "+str(highenergy[m])+" "+det1+" "+skip+" "+homedir)
##  for de in det2:
##    os.system("count_stat_new.py "+de+" "+homedir+"/Data"+de+"_"+mode+"_"+sep+"_"+str(lowenergy[m])+"_"+str(highenergy[m])+"keV.fits "+str(lowenergy[m])+" "+str(highenergy[m])+" "+sep+" "+homedir+" "+offx+" "+offy)




##for fil in ['fits_file','images']:
  #for k in range(len(fi)): # DO NOT UNCOMMENT THIS UNLESS YOU SET FI TO A LIST, IT IS BAD!!!!!!!!!!!
##    os.system('mkdir -p '+fi+'/'+fil)

#mo = ['full','sun','nosun']
mo = 'nosun'
##for i in range(len(fi)):
##os.system("mv "+homedir+"/*_"+mo+"*keV.fits "+fi+'/fits_file')
##os.system("mv "+homedir+"/Model*_"+mo+"*.fits "+fi+'/fits_file')
##os.system("mv "+homedir+"/*.png "+fi+'/images')
##os.system("mv "+homedir+"/*params* "+fi)

# The above lines use to have the first 2 os.system lines indented for the loop and mo and fi had index i. The last two lines had index:
# fi[mo.index(sep)], fancy way to send to appropriate directories. 
