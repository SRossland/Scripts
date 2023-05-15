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

# To change the type of file read in (full, sun, or nosun) see lines: 338 & 453
# These two lines change the type of file read in by either file name (453) or 
# by classified type (338, full, sun, nosun, etc.)
# Other examples are commented out so it's mostly a preference of which files to use
#################################
obsdir = sys.argv[1]
obslist = sys.argv[2]
mode = sys.argv[3]
elow = float(sys.argv[4])
ehigh = float(sys.argv[5])
detstodo = sys.argv[6]
ski = sys.argv[7]
homedir = sys.argv[8]

skip = ski.lower() in ['true','yes','1']

if detstodo == 'BOTH':
	det = ['A','B']
else: 
	det = [sys.argv[6]]
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

bad_obs = []

def robust_expomap(g_start, g_stop, time, idx, X, Y, roll, roll_time, x_ex, y_ex, rad_ex, obs):
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
  for i in range(len(idx)):
    dt = time[idx[i+1]]-time[idx[i]]
    if exp_target <= exp_total: return exp
    exp_total += dt
    for j in range(len(x_ex)):
      #temp_point = rotatePoint((X[i],Y[i]),(x_ex[j],y_ex[j]),-roll[i])
      offset_x = x_ex[j]-(X[idx[i]]-10)
      offset_y = y_ex[j]-(Y[idx[i]]-350)
      temp_point = rotatePoint((10,350),(offset_x,offset_y),-roll[idx[i]])
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
  if j_idx == 0:
      obs_event_det = obs_event.replace('01_cl.evt','_det1.fits')
      obs_event_attorb = obs_event.replace('01_cl.evt','.attorb')
      #obs_gti = obs_dir+'/'+obs+'/nu'+obs+det+'01_usrgti.fits'
      obs_gti = obs_event.replace('cl.evt','gti.fits')
  if j_idx != 0:
      root_event,crap = obs_event.split('sep')
      obs_event_det = root_event+'cl/nu'+obs+det+'_det1.fits'
      obs_event_attorb = root_event+'cl/nu'+obs+det+'.attorb'
      obs_gti = obs_event.replace('_fullevts_','01_gti_')
  obs_excl = obs_dir+'/'+obs+'/event_defcl/excl.reg'
  
  if (os.path.isfile(obs_gti) == False) or (os.path.isfile(obs_event_det) == False) or \
     (os.path.isfile(obs_event_attorb) == False) or (os.path.isfile(obs_excl) == False):
	return failed_array
  
  gti_start,gti_stop = get_gti(obs_gti)
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
          exp_map = robust_expomap(gti_start[i],gti_stop[i],time,idx,X,Y,roll,roll_time,x_ex,y_ex,rad_ex,obs)
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
	dt = gti_stop[i] - time[len(time)-1] #  Here, look here
      else:
        dt = time[ix+1]-time[ix] # A weird bit of code fixed here (; dt_mask = np.ones((360,360))) didn't need to be in the else
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
  with fits.open(fil,mode='update') as hdul:
    hdul[0].data = ary
    hdul[0].header['EXPOSURE'] += ex
    hdul[0].header['COMMENT'] = obsid
    hdul.flush()

def cut_time(obs_event, delta):
	with fits.open(obs_event) as hdul:
		timetocut = hdul[1].data['TIME']
	gti_file = obs_event.replace('_fullevts_','01_gti_')
	with fits.open(gti_file) as hdul:
		start = hdul[1].data['START']
		stop = hdul[1].data['STOP']
	start += delta
	stop -= delta
	new_exp = sum(stop-start)# trying to add
	idx_cut = np.zeros(len(timetocut), dtype=bool)
	for i in range(len(start)):
		idx_cut += (timetocut>=start[i])*(timetocut<=stop[i])
	return idx_cut, new_exp

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

#times = ['full','sun','nosun']  # All times have been changed, to go to original code, they need the index [j] added to them
#times = 'sun'
times = 'full'
soc = ['NOOCC','OCC']

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
            j = 0
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
           	    event = obsdir+'/'+obs[i]+'/event_cl/nu'+obs[i]+dets+mode+'_cl.evt'

	#	    else:      #The below line was indented   <--------********
	#	    event = obsdir+'/'+obs[i]+'/event_sep_cl/'+soc[k]+'/nu'+obs[i]+dets+mode+'_fullevts_'+times.upper()+'.fits'
		    # Check for source excl.reg file if mode is 01. The idea of this
		    # code is to create a mask where the region file is. This mask will be used
		    # in conjuction with the data aquired from the event files to use only 
		    # pixels that are declared 'ACTIVE'.
############ Here needs to be redesigned to pull the expmaps from the file:
# Simple approach, pull the file as an array, if it doesn't have it, create it.
#### First part is still the same, it looks at the excl.reg to exclude those PI counts
		    if mode == '01':
			reg_check = obsdir+'/'+obs[i]+'/event_defcl/excl.reg'
			if skip:
				if os.path.isfile(reg_check):
					continue
			if os.path.isfile(reg_check):
				print('EXISTS')
				srcfile = open(reg_check,'r+')
                                x_vec,y_vec,rad,data = [],[],[],[]
                                for line in srcfile.readlines():
                                        data.append(line.rstrip())
                                srcfile.close()
                                for q in range(len(data[3:])):
                                        xv, yv, radv = data[3+q].split(',')
                                        x_vec.append(xv[7:])
                                        y_vec.append(yv)
                                        rad.append(radv[:-1])
                                [float(t) for t in x_vec]
                                [float(t) for t in y_vec]
                                [float(t) for t in rad]

                                grid_pi = np.ones((1000,1000))
                                for q in range(len(x_vec)):
                                        mask = create_circ_mask(1000,1000,float(x_vec[q])-1,float(y_vec[q])-1,float(rad[q]))
                                        grid_pi[mask] = 0
				exp_correct = False
				nustar_cxb = obsdir.replace('OBS','CXB/'+mode+'/expmap/nu'+obs[i]+dets+'_'+times)
				if dets == 'A':
#### Here is where I am creating the exp map:
# adding an if statement, if a file exists, read it, if not, create it.
				  
				  if os.path.isfile(nustar_cxb):		#  added lines
					with fits.open(nustar_cxb) as hdul:     #
						exp_map_A = hdul[0].data        #
				  else:	
					print(obs[i]+' '+dets+' has no exp map done, needs to be done')
					continue
					#
					#exp_map_A = create_expomap(obs[i],obsdir,dets,j,event)
					#print('ExpA done')
					#if (exp_map_A == -42.0).all(): 
					#	print("Missing a file for exp corr. A")
					#	continue
					#fits.writeto(nustar_cxb,np.zeros((360,360)))		#
					#with fits.open(nustar_cxb,mode='update') as hdul:	#
					#	hdul[0].data = exp_map_A			#
					#	hdul[0].header['EXPOSURE'] = get_exp(event)	#
					#	hdul[0].header['COMMENT'] = obs[i]		#
					#	hdul.flush()					#
				else:
				  if os.path.isfile(nustar_cxb):		#
					with fits.open(nustar_cxb) as hdul:	#
						exp_map_B = hdul[0].data	#
				  else:	
					print(obs[i]+' '+dets+' has no exp map done, needs to be done')
					#
					#exp_map_B = create_expomap(obs[i],obsdir,dets,j,event)
					#print('ExpB done')
					#if (exp_map_B == -42.0).all(): 
					#	print("Missing a file for exp corr. B")
					#	continue
					#fits.writeto(nustar_cxb,np.zeros((360,360)))            #
                                        #with fits.open(nustar_cxb,mode='update') as hdul:       #
                                        #        hdul[0].data = exp_map_B                        #
                                        #        hdul[0].header['EXPOSURE'] = get_exp(event)     #
                                        #        hdul[0].header['COMMENT'] = obs[i]              #
                                        #        hdul.flush()
			
			else: 
				grid_pi = np.ones((1000,1000))
				exp_correct = True
				if dets == 'A':
					exp_map_A = np.ones((360,360))
				else: 
					exp_map_B = np.ones((360,360))
	
		    else:
				grid_pi = np.ones((1000,1000))
				exp_correct = True
				if dets == 'A':
					exp_map_A = np.ones((360,360))
				else:
					exp_map_B = np.ones((360,360))

##################### End of redesign
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
                                GRADE = hdul[1].data['GRADE']

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

				idx_good = (X > 0)*(Y > 0)*(GRADE <= 26)

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
				if exp_correct:
					expA[j] += exp_map_A*EXP*em
				else:
					expA[j] += exp_map_A*em
				exptimeA[j] += EXP
				
				fiA = homedir+'/DataA_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
  				exA = homedir+'/ExpA_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
				noA = homedir+'/NormA_'+mode+'_'+times+'_'+str(elow)+'_'+str(ehigh)+'keV.fits'
				update(fiA,imgA[j],obs[i],EXP)
				update(exA,expA[j],obs[i],EXP)
				normA[j] = np.divide(imgA[j],expA[j],out=np.zeros_like(imgA[j]),where=expA[j]!=0)*expA[j].max()
				update(noA,normA[j],obs[i],EXP)
				
			else:
				if exp_correct:
                                        expB[j] += exp_map_B*EXP*em
                                else:
                                        expB[j] += exp_map_B*em

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
