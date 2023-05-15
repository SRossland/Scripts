#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# syntax: ./image_exp_map.py obsdir obsids mode DET 

# Created by: Steven P. Rossland @ University of Utah

import numpy as np, os, sys, datetime, math
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord as sc
import matplotlib.pyplot as plt
#################################
# UPDATE 10/15/2020: This code is an altered version of image_stack.py. This code will create individual
# images of the exp for all observations given, it will even create ones for obs without the excl.reg file
# Files will be stored in a directory in the Xray directory under Dan's working site.

# This program works as intended.  It creates 6 images, though I am still working on the norm one.
# As a way to check my work, I am applying a line that will let you skip the excl.reg 
# obs so you can see just blank sky. I can then compare the two images used and see what is going on

#################################
obsdir = sys.argv[1]
obslist = sys.argv[2]
mode = sys.argv[3]

if len(sys.argv) == 4:
	det = ['A','B']
elif len(sys.argv) != 5:
	print('Make sure you provide only the information needed!')
	sys.exit()
else: 
	det = [sys.argv[4]]
##################################

NuSTAR_CXB_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/CXB/'+mode+'/expmap/'
mask_dir = '/uufs/astro.utah.edu/common/home/u1019304/temp'
with fits.open(mask_dir+'/fullmaskA_final.fits') as mAhdu:
    edgemaskA = mAhdu[0].data

with fits.open(mask_dir+'/fullmaskB_final.fits') as mBhdu:
    edgemaskB = mBhdu[0].data



expAfull = np.zeros((360,360)); expAsun = np.zeros((360,360)); expAnosun = np.zeros((360,360))
expBfull = np.zeros((360,360)); expBsun = np.zeros((360,360)); expBnosun = np.zeros((360,360))

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

# Pull the files
# create the expomap (remember it needs to be flipped over the y-axis and rotated)
# map to skycoords
# exclude excl.reg
# map back to det1
# sum and mask by edgemask 

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

def get_hk(fi):
  with fits.open(fi) as hdul:
    time = hdul[1].data['TIME']
    livetime = hdul[1].data['LIVETIME']
  return time, livetime

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

bad_obs = []

def robust_expomap(g_start, g_stop, time, idx, X, Y, roll, roll_time, x_ex, y_ex, rad_ex, obs, hk_time, liv_time):
  # This is a robust version of creating the expo_map when the difference between the
  # gti time gap and the total time accouted for in the det1 times are greater than
  # 1.2 seconds.  This will identify where the extra time will be added onto the map
  # and how much time that needs to be added. Now there is some limit to how much time 
  # can be added.
  # This will be done by the individual gti gaps, so it will return the expo_map to be added
  # to the overall expo_map and let the other code keep running.
  
  # First adjust the time to the start if any exist:
  print('ROBUST')
  exp = np.zeros((360,360))
  diff_start = time[idx[0]]-g_start; diff_stop = g_stop-time[idx[-1]]
  exp_target = g_stop-g_start; exp_total = 0; dt_mask = np.ones((360,360))
  if diff_start >= 0.2:
    for j in range(len(x_ex)):
      #temp_point = rotatePoint((X[0],Y[0]),(x_ex[j],y_ex[j]),-roll[0])
      offset_x = x_ex[j]-(X[0]-10)
      offset_y = y_ex[j]-(Y[0]-350)
      #offset_x = temp_point[0]-9
      #offset_y = temp_point[0]-349
      temp_point = rotatePoint((10,350),(offset_x,offset_y),-roll[0])
      # Here we flip the point over the y-axis, this only need to be done for X coords
      off_y = temp_point[1]-1
      offset_x_flip = 360-temp_point[0]-1
      # This should be the coords of the point in DET1
      mask = circ_mask(360,360,offset_x_flip,off_y,rad_ex[j])
      dt_mask[mask] = 0
    exp += dt_mask*diff_start
    exp_total += diff_start
  for m, ix in enumerate(idx):
  #for i in range(len(idx)):
    if exp_target <= exp_total: return exp
    if int(time[ix]) >= roll_time[-1]:                            #####
      offset_roll = offset_roll
    else:
      offset_roll = np.where(roll_time >= int(time[ix]))[0][0]
    if ix >= len(time)-1:
	dt = g_stop - time[len(time)-1]
	idx_hk = (np.abs(hk_time-time[len(time)-1])).argmin()
    else:
	dt = time[ix+1]-time[ix]
	idx_hk = (np.abs(hk_time-time[ix])).argmin()              ####

    #dt = time[idx[i+1]]-time[idx[i]]
    #if exp_target <= exp_total: return exp
    exp_total += dt
    dt *= liv_time[idx_hk]
    for j in range(len(x_ex)):
      #temp_point = rotatePoint((X[i],Y[i]),(x_ex[j],y_ex[j]),-roll[i])
      offset_x = x_ex[j]-(X[ix]-10)
      offset_y = y_ex[j]-(Y[ix]-350)
      temp_point = rotatePoint((10,350),(offset_x,offset_y),-roll[int(offset_roll)])
      #offset_x = temp_point[0]-9
      #offset_y = temp_point[0]-349
      # Here we flip the point over the y-axis, this only need to be done for X coords
      #offset_x_flip = 359-offset_x
      off_y = temp_point[1]-1
      offset_x_flip = 360-temp_point[0]-1
      # This should be the coords of the point in DET1
      mask = circ_mask(360,360,offset_x_flip,off_y,rad_ex[j])
      dt_mask[mask] = 0
    exp += dt_mask*dt
    exp_total += dt
  remainder = exp_target-exp_total
  if remainder > 20.0: 
    print("Robust failer with time remaining of: "+str(remainder))
    print("Ended on observation: "+str(obs))
    bad_obs.append(obs)
  exp += dt_mask*(exp_target-exp_total)
  return exp


def create_expomap(obs, obs_dir, det, j_idx, obs_event):
  print(obs)
  exp_grid = np.zeros((360,360))
  repeat = 0
  failed_array = np.ones((3,3))*-42.0
  # It is assumed that a file (excl.reg) exists in a subdirectory named event_defcl for the obs
  # Do not forget the translation needed between NuSTAR and python arrays
  #obs_hk = obs_event.replace(mo+'_cl.evt','_fpm.hk')##############
  if j_idx == 0:
      obs_event_det = obs_event.replace('01_cl.evt','_det1.fits')
      obs_hk = obs_event.replace('01_cl.evt','_fpm.hk')  ###############
      obs_event_attorb = obs_event.replace('01_cl.evt','.attorb')
      obs_gti = obs_dir+'/'+obs+'/nu'+obs+det+'01_usrgti.fits'
      #obs_gti = obs_event.replace('cl.evt','gti.fits')
  if j_idx != 0:
      root_event,crap = obs_event.split('sep')
      obs_hk = root_event+'cl/nu'+obs+det+'_fpm.hk' ##############
      obs_event_det = root_event+'cl/nu'+obs+det+'_det1.fits'
      obs_event_attorb = root_event+'cl/nu'+obs+det+'.attorb'
      obs_gti = obs_event.replace('_fullevts_','01_gti_')
  obs_excl = obs_dir+'/'+obs+'/event_defcl/excl.reg'

  if (os.path.isfile(obs_gti) == False) or (os.path.isfile(obs_event_det) == False) or \
     (os.path.isfile(obs_event_attorb) == False) or (os.path.isfile(obs_excl) == False):
	return "NO FILE" 

  hk_time, livetime = get_hk(obs_hk)            ################
  gti_start,gti_stop = get_gti(obs_gti)
  #if j_idx == 2:
#	gti_start += 300
#	gti_stop -= 300
  time, X, Y = get_detpos(obs_event_det)
  roll, roll_time = get_roll(obs_event_attorb)
  x_ex, y_ex, rad_ex = getExcl(obs_excl)

  # This should have gathered all of the information, now to get checks and balances
  for i in range(len(gti_start)):
    tick = 0
    dt_o = gti_stop[i] - gti_start[i]
    idx_array = np.where(np.logical_and(time>=gti_start[i],time<=gti_stop[i]))
    idx = idx_array[0]
    # get average time gap in det recordings
    if len(idx) > 1: avg_dt = np.average(np.diff(time[idx]))
    elif len(idx) == 1: avg_dt = 1
    else: continue
    # check time gaps
    time_gap = time[idx[-1]]-time[idx[0]]
    if dt_o > time_gap:
	if np.abs(dt_o - time_gap) > 1.2:
	  print("Large time difference between GTI and det1 recording.")
          print(" delta time: "+str(dt_o-time_gap))
          exp_map = robust_expomap(gti_start[i],gti_stop[i],time,idx,X,Y,roll,roll_time,x_ex,y_ex,rad_ex,obs,hk_time,livetime)
          exp_grid += exp_map
          continue
        else:
          repeat = int(dt_o - time_gap)/avg_dt
    # Find the translated point of excl in det1 coords (remember, the 350,350 point is flipped 
    # around the y axis.
    # so the 350,350 point is ref[349,349] in a python grid)
    # But keep in mind, the values pulled from NuSTAR are absolute (i.e. X[i] = 426.## is the actual
    # physical pixel so it would be X[i]-1 in the grid)
    
    # Find the new position
    # Check for roll difference:
    for m,ix in enumerate(idx):
      if int(time[ix]) >= roll_time[-1]: 
	offset_roll = offset_roll
      else:
	offset_roll = np.where(roll_time >= int(time[ix]))[0][0]
      #elif (time[ix] != roll_time[ix]) and (time[ix] < roll_time[-1]):
        #offr = np.where(roll_time == int(time[ix]))[0] 
        ##offr = np.where(np.logical_and(roll_time >= time[idx[m]]-0.2, roll_time <= time[idx[m]]+0.2))
        #offset_roll = offr[0]-ix
      ##elif ix >= len(roll_time): offset_roll = offset_roll
      #else: offset_roll = 0
      pox,poy,por = [],[],[]
      if ix >= len(time)-1:
	dt = gti_stop[i] - time[len(time)-1]
        idx_hk = (np.abs(hk_time-time[len(time)-1])).argmin() ############
      else:
        dt = time[ix+1]-time[ix] 
      # Here we do livetime correct:
        idx_hk = (np.abs(hk_time - time[ix])).argmin()        ############
      dt *= livetime[idx_hk]                                  ############
      dt_mask = np.ones((360,360))
      tick += dt
      if tick > dt_o: continue # did we exceed gti time?
      for j in range(len(x_ex)):
        #temp_point = rotatePoint((X[ix],Y[ix]),(x_ex[j],y_ex[j]),-roll[int(offset_roll)])
	off_x = x_ex[j]-(X[ix]-10)
        #offset_x = temp_point[0]-(X[ix]-10)
	offset_y = y_ex[j]-(Y[ix]-350)
        #offset_y = (temp_point[1]-(Y[ix]-350))-1
      # Here we flip the point over the y-axis, this only need to be done for X coords
        #offset_x_flip = (360-offset_x)-1
	temp_point = rotatePoint((10,350),(off_x,offset_y),-roll[int(offset_roll)])
      # This should be the coords of the point in DET1
	off_y = temp_point[1]-1
	off_x_flip = 360-temp_point[0]-1
    # Here we sill make a mask of where the excl.reg is and add the time to the remaining pixels
        mask = circ_mask(360,360,off_x_flip,off_y,rad_ex[j])
        dt_mask[mask] = 0
      exp_grid += dt_mask*dt
    # Now I need to ensure we have all the time accounted for by using the last mask and remaining time
    if repeat >= 0: 
      dt_r = dt_o - tick
      exp_grid += dt_mask*dt_r
  return exp_grid  

def update(fil, ary, obsid, ex):
  fits.writeto(fil,np.zeros((360,360)))
  with fits.open(fil,mode='update') as hdul:
    hdul[0].data = ary
    hdul[0].header['EXPOSURE'] = ex
    hdul[0].header['COMMENT'] = obsid
    hdul.flush()

#################################
if mode == '01':
  BAD_PNT = 10.0
  k = 0
elif mode == '02':
  BAD_PNT = 0.0
  k = 1
else:
  print('Need a mode type (01 or 02):')
  sys.exit()

times = ['full','sun','nosun']
soc = ['NOOCC','OCC']

expA = [expAfull,expAsun,expAnosun]
expB = [expBfull,expBsun,expBnosun]

eAf,eAs,eAn = 0,0,0
eBf,eBs,eBn = 0,0,0

obs_count = 0
###

#Here we want to read in those files to the arrays to the loop doesn't have to change that much
# I know this is ugly, i'll try to clean it up later


try:

	obs = np.genfromtxt(obslist,dtype="|U11") # was named obs_orig
	#print(len(obs))
	for dets in det:
	  for j in range(3):
	    for i in range(len(obs)):
		    # Here, check if the exp map already exists:
		    if os.path.isfile(NuSTAR_CXB_dir+'nu'+obs[i]+dets+'_'+times[j]):
			print("Exp map already exists for "+obs[i])
			continue
		
		    if j == 0: 
		        event = obsdir+'/'+obs[i]+'/event_cl/nu'+obs[i]+dets+mode+'_cl.evt'
		    else:
			event = obsdir+'/'+obs[i]+'/event_sep_cl/'+soc[k]+'/nu'+obs[i]+dets+mode+'_fullevts_'+times[j].upper()+'.fits'
		    # Check for source excl.reg file if mode is 01. The idea of this
		    # code is to create a mask where the region file is. This mask will be used
		    # in conjuction with the data aquired from the event files to use only 
		    # pixels that are declared 'ACTIVE'.
		    if mode == '01':
			reg_check = obsdir+'/'+obs[i]+'/event_defcl/excl.reg'
			if os.path.isfile(reg_check):
				print('EXISTS')

				exp_correct = False
				if dets == 'A':
					exp_map_A = create_expomap(obs[i],obsdir,dets,j,event)
					print('ExpA done')
					if (exp_map_A == -42.0).all(): 
						print("Missing a file for exp corr. A")
						continue
				else:
					exp_map_B = create_expomap(obs[i],obsdir,dets,j,event)
					print('ExpB done')
					if (exp_map_B == -42.0).all(): 
						print("Missing a file for exp corr. B")
						continue
				
			else:
				exp_correct = True
				if dets == 'A':
					exp_map_A = np.ones((360,360))
				else:
					exp_map_B = np.ones((360,360))

		    if os.path.isfile(event):
			if dets == 'A': em = np.copy(edgemaskA)
			else: em = np.copy(edgemaskB)
			em = em.astype('float64')
			#if j == 2:
			#	nosun_gti_event = event.replace('_fullevts_','01_gti_')
			#	gti_start_E,gti_stop_E = get_gti(nosun_gti_event)
			#	with fits.open(event) as hdul:
			#		RA = hdul[0].header['RA_PNT']
			#		DEC = hdul[0].header['DEC_PNT']
			#	gti_start_E += 300
			#	gti_stop_E -= 300
			#	EXP = np.sum(gti_stop_E-gti_start_E)
			#else:
			with fits.open(event) as hdul:
				EXP = hdul[1].header['EXPOSURE']
				RA = hdul[0].header['RA_PNT']
                                DEC = hdul[0].header['DEC_PNT']
			c = sc(ra = RA*u.degree, dec = DEC*u.degree, frame='fk5')
			if abs(c.galactic.b.deg) < BAD_PNT:
				continue

			
			if dets == 'A':
				if exp_correct:
					#expA[j] += exp_map_A*EXP*em
					expmap = exp_map_A*EXP*em
				else:
					#expA[j] += exp_map_A*em
					expmap = exp_map_A*em
   # This is the old file position:
   				exA = NuSTAR_CXB_dir+'nu'+obs[i]+dets+'_'+times[j]
   # This will be a test file position pointed at my home directory:
#				exA = '/uufs/astro.utah.edu/common/home/u1019304/nu'+obs[i]+dets+'_'+times[j]
				expmap = np.where(expmap<0,0.,expmap)
				expmap = (expmap/np.max(expmap))*EXP
				update(exA,expmap,obs[i],EXP)
				
			else:
				if exp_correct:
                                        #expB[j] += exp_map_B*EXP*em
					expmap = exp_map_B*EXP*em
                                else:
                                        #expB[j] += exp_map_B*em
					expmap = exp_map_B*em
#				exB = '/uufs/astro.utah.edu/common/home/u1019304/nu'+obs[i]+dets+'_'+times[j]
				exB = NuSTAR_CXB_dir+'nu'+obs[i]+dets+'_'+times[j]
				expmap = np.where(expmap<0,0.,expmap)
				expmap = (expmap/np.max(expmap))*EXP
                                update(exB,expmap,obs[i],EXP)
			em = None; expmap = None
			obs_count += 1
except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(str(exc_type) + '\n' + str(e) +'\n' +str(exc_tb.tb_lineno)+ '\n')
