#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# syntax: ./image_stack.py obsdir obsids mode el_low el_high DET 

import numpy as np, os, sys, datetime, math
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord as sc
#################################
# This program works as intended.  It creates 6 images, though I am still working on the norm one.
# As a way to check my work, I am applying a line that will let you skip the excl.reg 
# obs so you can see just blank sky. I can then compare the two images used and see what is going on

#skip = True # This will skip any obs that has a excl.reg file.
#skip = False
#################################
obsdir = sys.argv[1]
obslist = sys.argv[2]
elow = float(sys.argv[3])
ehigh = float(sys.argv[4])
detstodo = sys.argv[5]
homedir = sys.argv[6]


if detstodo == 'BOTH':
	det = ['A','B']
else: 
	det = [sys.argv[5]]
##################################

mask_dir = '/uufs/astro.utah.edu/common/home/u1019304/temp'
with fits.open(mask_dir+'/fullmaskA_final.fits') as mAhdu:
    edgemaskA = mAhdu[0].data

with fits.open(mask_dir+'/fullmaskB_final.fits') as mBhdu:
    edgemaskB = mBhdu[0].data

#edgemaskA, edgemaskB = np.ones((360,360)), np.ones((360,360))
el_choi = (elow-1.6)/0.04
hl_choi = (ehigh-1.6)/0.04

# These next lines should be checked by seeing if the files already exist, if they do, then read those in.
Afull = np.zeros((360,360)); Asun = np.zeros((360,360)); Anosun = np.zeros((360,360))
Bfull = np.zeros((360,360)); Bsun = np.zeros((360,360)); Bnosun = np.zeros((360,360))

expAfull = np.zeros((360,360)); expAsun = np.zeros((360,360)); expAnosun = np.zeros((360,360))
expBfull = np.zeros((360,360)); expBsun = np.zeros((360,360)); expBnosun = np.zeros((360,360))

normAfull = np.zeros((360,360)); normAsun = np.zeros((360,360)); normAnosun = np.zeros((360,360))
normBfull = np.zeros((360,360)); normBsun = np.zeros((360,360)); normBnosun = np.zeros((360,360))
##################################
def create_circ_mask(h,w,centerx,centery,radius):
	# Here is am again using X and Y, but they are local and should be treated as such by python
	Y,X = np.ogrid[:h,:w]
	dist_from_center = np.sqrt((X-centerx)**2 + (Y-centery)**2)
        mask = dist_from_center <= radius
        return mask
#################################

''' This is for creating the translation between det1 and skycoords to map the excl.reg'''
''' Need to pull the _det.fits, .attorb, usr.gti files.  The det1 is a file that maps'''
''' pixel 350,350 from det1 to sky.  This should let me create a full pixel map from det1'''
''' to skycoords that is the expomap, I can do it either by building up the time, or I can'''
''' create a map of the full time and exclude the dt from the excl.reg.'''
###########   IMPORTANT: MOST OF THE FOLLOWING CODE IS OBSOLETE AND HAS BEEN MOVED TO ANOTHER CODE (EXP_STACK.PY OR SOMETHING)
# Pull the files
# create the expomap (remember it needs to be flipped over the y-axis and rotated)
# map to skycoords
# exclude excl.reg
# map back to det1
# sum and mask by edgemask 
def get_exp(fi, val = 1):
	with fits.open(fi) as hdul:
		ret = hdul[val].header['EXPOSURE']
	return ret

def get_data(fi, key, data='data', val=1): #This is a general fetch def to get data from any fits file
  with fits.open(fi) as hdul:
    ret = hdul[val].data[key]
  return ret

def get_gti(fi):  #Gets the start/stop time from the gti files in the obs dir
  with fits.open(fi) as hdul:
    gti_start = hdul[1].data['START']
    gti_stop = hdul[1].data['STOP']
  return gti_start,gti_stop

def get_detpos(fi): #Gets the reference pixel coords and time for the DET1 files in obs/event_cl
  with fits.open(fi) as hdul:
    time = hdul[1].data['TIME']
    X = hdul[1].data['X_DET1']
    Y = hdul[1].data['Y_DET1']
  return time,X,Y

def get_roll(fi): # Gets the Roll angle and roll time (for reference) from the attorb file in obs/event_cl
  with fits.open(fi) as hdul:
    roll = hdul[1].data['ROLL']
    roll_time = hdul[1].data['TIME']
  return roll, roll_time

def rotatePoint(centerPoint,point, angle): # rotates the center pixel of the excl.reg back to eliminate the roll 
  angle = math.radians(angle)
  temp_point = point[0]-centerPoint[0] , point[1]-centerPoint[1]
  temp_point = ( temp_point[0]*math.cos(angle)-temp_point[1]*math.sin(angle) , temp_point[0]*math.sin(angle)+temp_point[1]*math.cos(angle))
  temp_point = temp_point[0]+centerPoint[0] , temp_point[1]+centerPoint[1]

  return temp_point

def getExcl(fi): # Retrieves the X Y and radius of the excl.reg file (assuming physical!!!!)
  f = open(fi,'r+')
  blah = []
  for line in f.readlines():
    blah.append(line)
  x = []
  y = []
  rad = []
  for i in range(3, len(blah)):
      pos = blah[i][7:-2].rstrip().split(',')
      x.append(float(pos[0]))
      y.append(float(pos[1]))
      rad.append(float(pos[2]))
  return x,y,rad

expomap = np.zeros((1000,1000))
# Here's the idea, get the point of the 350,350 pixel in sky coords, get the point of the center of the excl.reg
# then rotate that point around the 350,350 back and then translate it back to det1 coords. Using that point
# create the excl.reg on the det1 img for that time period and then repeat. I know this is inefficient, but 
# I need it to work now.

def circ_mask(h, w, centerx, centery,radius):
        Y, X = np.ogrid[:h,:w]
        dist_from_center = np.sqrt((X-centerx)**2 + (Y-centery)**2)
        mask = dist_from_center <= radius
        return mask

def update(fil, ary, obsid, ex):
  with fits.open(fil,mode='update') as hdul:
    hdul[0].data = ary
    hdul[0].header['EXPOSURE'] += ex
    hdul[0].header['COMMENT'] = obsid
    hdul.flush()


#################################
mode = '02'
BAD_PNT = 0.0
k = 1

#times = ['full','sun','nosun']  # All times have been changed, to go to original code, they need the index [j] added to them
times = 'sun'
soc = 'OCC'

imgA = [Afull,Asun,Anosun]
imgB = [Bfull,Bsun,Bnosun]

expA = [expAfull,expAsun,expAnosun]
expB = [expBfull,expBsun,expBnosun]

normA = [normAfull,normAsun,normAnosun]
normB = [normBfull,normBsun,normBnosun]

eAf,eAs,eAn = 0,0,0
eBf,eBs,eBn = 0,0,0

exptimeA = [eAf,eAs,eAn]
exptimeB = [eBf,eBs,eBn]

obsAf, obsAs, obsAn = [],[],[]
obsBf, obsBs, obsBn = [],[],[]

obsA = [obsAf, obsAs, obsAn]
obsB = [obsBf, obsBs, obsBn]

obs_count = 0

failed_array = np.ones((3,3))*-42.0
###
# Here I will create a blank image to update during each iteration of the loop if an image does not already exist.

#curr = os.getcwd()
# replaced curr with homedir in the fi statement below
tpe = ['Data','Exp','Norm']

for i in range(3):
  for j,typ in enumerate(tpe):
    for dets in det:
      fi = homedir+'/'+typ+dets+'_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
      if os.path.isfile(fi): continue
      else:  
	fits.writeto(fi,np.zeros((360,360)))
	with fits.open(fi,mode='update') as hdul:
		hdul[0].header['COMMENT'] = 'OBSIDS'
		hdul[0].header['EXPOSURE'] = 0
		hdul.flush()
'''  
  BlankData.writeto('./DataA_'+mode+'_'+times[i]+'_'+str(int(elow))+'_'+str(int(ehigh))+'keV.fits')
  BlankData.writeto('./DataB_'+mode+'_'+times[i]+'_'+str(int(elow))+'_'+str(int(ehigh))+'keV.fits')
  BlankData.writeto('./ExpA_'+mode+'_'+times[i]+'_'+str(int(elow))+'_'+str(int(ehigh))+'keV.fits')
  BlankData.writeto('./ExpB_'+mode+'_'+times[i]+'_'+str(int(elow))+'_'+str(int(ehigh))+'keV.fits')
  BlankData.writeto('./NormA_'+mode+'_'+times[i]+'_'+str(int(elow))+'_'+str(int(ehigh))+'keV.fits')
  BlankData.writeto('./NormB_'+mode+'_'+times[i]+'_'+str(int(elow))+'_'+str(int(ehigh))+'keV.fits')
'''
#Here we want to read in those files to the arrays to the loop doesn't have to change that much
# I know this is ugly, i'll try to clean it up later

for i in range(3):
  fiA = homedir+'/DataA'+'_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
  fiB = homedir+'/DataB'+'_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
  ExA = homedir+'/ExpA_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
  ExB = homedir+'/ExpB_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
  NoA = homedir+'/NormA_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
  NoB = homedir+'/NormB_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'

  if (detstodo == 'A') or (detstodo == 'BOTH'):
    with fits.open(fiA,mode='update') as hdul:
      hdul[0].header['DATE'] = str(datetime.datetime.now())
      imgA[i] = hdul[0].data
      obsA[i] = hdul[0].header['COMMENT']
      hdul.flush()
    with fits.open(ExA,mode='update') as hdul:
      hdul[0].header['DATE'] = str(datetime.datetime.now())
      expA[i] = hdul[0].data
      hdul.flush()
    with fits.open(NoA,mode='update') as hdul:
      hdul[0].header['DATE'] = str(datetime.datetime.now())
      normA[i] = hdul[0].data
      hdul.flush()

  if (detstodo == 'B') or (detstodo == 'BOTH'):
    with fits.open(fiB,mode='update') as hdul:
      hdul[0].header['DATE'] = str(datetime.datetime.now())
      imgB[i] = hdul[0].data
      obsB[i] = hdul[0].header['COMMENT']
      hdul.flush()
    with fits.open(ExB,mode='update') as hdul:
      hdul[0].header['DATE'] = str(datetime.datetime.now())
      expB[i] = hdul[0].data
      hdul.flush()
    with fits.open(NoB,mode='update') as hdul:
      hdul[0].header['DATE'] = str(datetime.datetime.now())
      normB[i] = hdul[0].data
      hdul.flush()


try:

	obs_orig = np.genfromtxt(obslist,dtype="|U11")
	for dets in det:
	  #print(len(obs_orig))
	  #for j in range(3):  #I think this should be just j = 2
            j = 2
	    if dets == 'A':
	      A_obs = np.copy(obsA[j]); A_obs = list(A_obs)
              if 'OBSIDS' in A_obs: A_obs.remove('OBSIDS')
              obs = list(set(obs_orig)-set(A_obs))
	      #print(len(obs),len(A_obs))
            else:
	      B_obs = np.copy(obsB[j]); B_obs = list(B_obs)
              if 'OBSIDS' in B_obs: B_obs.remove('OBSIDS')
              obs = list(set(obs_orig)-set(B_obs))
	    for i in range(len(obs)):
	###******* Here is a shortcut to ensure the files are only nosun, commented out is the old code that did all three	
	#	    if j == 0: 
	#	    event = obsdir+'/'+obs[i]+'/event_cl/nu'+obs[i]+dets+mode+'_cl.evt'

	#	    else:      #The below line was indented   <--------********
		    event = obsdir+'/'+obs[i]+'/event_sep_cl/'+soc+'/nu'+obs[i]+dets+mode+'_fullevts_'+times.upper()+'.fits'
		    # Check for source excl.reg file if mode is 01. The idea of this
		    # code is to create a mask where the region file is. This mask will be used
		    # in conjuction with the data aquired from the event files to use only 
		    # pixels that are declared 'ACTIVE'.
############ Here needs to be redesigned to pull the expmaps from the file:
# Simple approach, pull the file as an array, if it doesn't have it, create it.
#### First part is still the same, it looks at the excl.reg to exclude those PI counts
		    grid_pi = np.ones((1000,1000))
		    if dets == 'A':
			    exp_map_A = np.ones((360,360))
		    else: 
			    exp_map_B = np.ones((360,360))
	

##################### End of redesign
                    if not os.path.isfile(event): continue
		    if os.path.isfile(event):
			if dets == 'A': em = np.copy(edgemaskA)
			else: em = np.copy(edgemaskB)
			em = em.astype('float64')
			with fits.open(event) as hdul:
                                PI = hdul[1].data['PI']
				EXP = hdul[1].header['EXPOSURE']
				RA = hdul[0].header['RA_PNT']
                                DEC = hdul[0].header['DEC_PNT']
                                X = hdul[1].data['X']
                                Y = hdul[1].data['Y']
                                DET1_X = hdul[1].data['DET1X']
                                DET1_Y = hdul[1].data['DET1Y']

#			if j == 2:					  # this is the same as last comment
#				idx_cut, new_exp = cut_time(event, 300) # here too, the 300 is the delta of time
#				EXP = None; EXP = new_exp

			c = sc(ra = RA*u.degree, dec = DEC*u.degree, frame='fk5')
			if abs(c.galactic.b.deg) >= BAD_PNT:
				expmap = np.zeros((360,360))
				
#				if j == 2:
#					idx_XY = (X != -1)*(Y != -1) 
#					idx_good = (idx_XY == 1)*(idx_cut == 1) # fix this #########
#				else: 

				idx_good = (X > 0)*(Y > 0)

				for m in range(len(PI[idx_good])): #NEW and futher [idx_good] index was added to other indicies
				    elow_int = int(round((elow-1.6)/0.04))
				    ehigh_int = int(round((ehigh-1.6)/0.04))
				    if elow_int <= PI[idx_good][m] < ehigh_int:  # HERE is the index change PLUS -- >
# 					if (X[m] == -1) or (Y[m] == -1) or (em[DET1_Y[m]-1,DET1_X[m]-1] == 0) or (grid_pi[X[m]-1,Y[m]-1] == 0):
					if (em[DET1_Y[idx_good][m]-1,DET1_X[idx_good][m] -1] == 0) or (grid_pi[Y[idx_good][m]-1,X[idx_good][m]-1] == 0):  # HERE # NOTICE !!!!!!!!!! IMPORTANT !!!!!!!!!!!!!! for some reason the indexing of grid_pi was X,Y when in python arrays it is indexed Y,X. I can not find a reason as to why I had it this way, it appears as if it was just a mistake, and not some reversal somehow. 
					    continue
					else:
					    if dets == 'A':
						imgA[j][DET1_Y[idx_good][m]-1,DET1_X[idx_good][m]-1] += 1  # HERE
					    else:
						imgB[j][DET1_Y[idx_good][m]-1,DET1_X[idx_good][m]-1] += 1  # HERE

			
			if dets == 'A':
				expA[j] += exp_map_A*EXP*em
				
				fiA = homedir+'/DataA_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
  				exA = homedir+'/ExpA_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
				noA = homedir+'/NormA_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
				update(fiA,imgA[j],obs[i],EXP)
				update(exA,expA[j],obs[i],EXP)
				normA[j] = np.divide(imgA[j],expA[j],out=np.zeros_like(imgA[j]),where=expA[j]!=0)*expA[j].max()
				update(noA,normA[j],obs[i],EXP)
				
			else:
                                expB[j] += exp_map_B*EXP*em

                                fiB = homedir+'/DataB_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
                                exB = homedir+'/ExpB_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
                                noB = homedir+'/NormB_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
                                update(fiB,imgB[j],obs[i],EXP)
                                update(exB,expB[j],obs[i],EXP)
				normB[j] = np.divide(imgB[j],expB[j],out=np.zeros_like(imgB[j]),where=expB[j]!=0)*expB[j].max()
                                update(noB,normB[j],obs[i],EXP)

			em = None
			obs_count += 1


except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(str(exc_type) + '\n' + str(e) +'\n' +str(exc_tb.tb_lineno)+ '\n')
