#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# syntax: ./spec_stack.py

# Purpose:  The purpose of this script is to stack together, the counted data
# from a given list of observations from the OBSID#AorB_01.evt and 02 files and
# with the total exposure time into a .pha file to be read into XSPEC
# This will be done for full events, no sun, and sun tagged events.

# This program is to be run with the file associated with the stack, in the absense of such file
# the program will create a new one in the programed directory ~line 51

# Created by: Steven P. Rossland, June 2018

# So I started to write the code for stacking all of the stacks at once
# but without a clear file idea of the event files needed for sun or occ
# I have decided to wait on finishing this part until I get some more info
# Ideally I will have the code checking for the existance for each of the files
# and writing what obsids are stacked into each spec.  maybe created 6 or 12 lists 
# or one large file with a boolan table showing which obs are in which stack

# 3rd revamp:
# Adding a conditional option that will limit the exposure time of the stack
# This will create a new stack limited purely by that given exposure on the 
# idea that you can look at similar rates for two stacks so one is not 
# statistically dominate.

# 4th:
# I am now clipping out edge pixels with the edgemask in the auxil dire
# of the NuSTAR dir in my home dir. So far it is just implemented with 
# the 02's in mind.  To do this for the 01's, it will take some consideration
# due to the excl.reg file that may be used.
#
# 5th:
# Added an option to limit the stacks by full, nosun, sun
# also included an option to make an image with an added 
# limit of energy

import numpy as np, os, sys, datetime
from astropy.io import fits
# The following two lines are only for converting coords
from astropy import units as u
from astropy.coordinates import SkyCoord as sc

#def main():
obsdir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS'
obslist = raw_input("OBSID list (with file extention):  ")
dirlist = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/SolarData'

############## 4th
# Create a mask for clipping out the edge pixels that contribute to higher counts than reality

#mask_dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/auxil'
mask_dir = '/uufs/astro.utah.edu/common/home/u1019304/temp'
with fits.open(mask_dir+'/fullmaskA_final2.fits') as mAhdu:
    edgemaskA = mAhdu[0].data
    
with fits.open(mask_dir+'/fullmaskB_final2.fits') as mBhdu:
    edgemaskB = mBhdu[0].data
############## end 4th
############## 5th
# The edition of version of Sun, nosun, full and if you want an image.
print('Choose:  [1] all')
print('         [2] full')
print('         [3] nosun')
print('         [4] sun')
ith = raw_input(':')

if ith == '1': i_choi = range(3)
elif ith == '2': i_choi = [0]
elif ith == '3': i_choi = [1]
elif ith == '4': i_choi = [2]

img_choi = raw_input('Create an image with the stack? [Y/N] *caps are important')
if img_choi == 'Y':
	el_c = float(eval(raw_input('Lower energy limit (in energy): ')))
	hl_c = float(eval(raw_input('High energy limit (in energy): ')))
	el_choi = (el_c-1.6)/0.04
	hl_choi = (hl_c-1.6)/0.04
	TestAfull = np.zeros((360,360))
	TestBfull = np.copy(TestAfull)
	TestAnosun = np.copy(TestAfull)
	TestBnosun = np.copy(TestAfull)
	TestAsun = np.copy(TestAfull)
	TestBsun = np.copy(TestAfull)
	imgApack = [TestAfull,TestAnosun,TestAsun]
	imgBpack = [TestBfull,TestBnosun,TestBsun]
	imgAnames = ['DataA_FULL','DataA_NOSUN','DataA_SUN']
	imgBnames = ['DataB_FULL','DataB_NOSUN','DataB_SUN']
############## end 5th
# Asks if 10.0 degrees is still a good limit for exclusion around the galactic plane

BAD_PNT = 10.0
check_pnt = raw_input("Current GAL. Coord. latitude limit is set at 10.0 degrees, would you like to change?  [y,n]")
if check_pnt == 'y':
	PNT = raw_input("What is your new limit (in Gal. Coords. in degrees)?  ")
	BAD_PNT = float(PNT)

limit = raw_input("Would you like to limit the exposure time to another stack or obs? [y,n]")
if limit == 'y':
	explimit_file_A = raw_input("limiting file: (must be located in spec directory):")
	det = raw_input("Detector? (capital A or B)")
	print("Need some more info")
	print("  Is this: [1] Full spectrum (with sun and no sun combined)")
  	print("           [2] No Sun spectrum")
	print("           [3] Sun spectrum")
	lim_spec = int(raw_input("spec:"))
	print("  Is this: [1] 01 file")
	print("           [2] 02 file")
	lim_evt = int(raw_input("evt:"))
	if os.path.isfile(dirlist+'/spec/'+explimit_file_A) == False:
		print("No valid limiting file was given, Exiting!")
		sys.exit()

# To STACK WHAHAHAHAHA
specs = ['full','nosun','sun']
evt = ['01','02']
soc = ['NOOCC','OCC']
###################################
def create_circ_mask(h, w, centerx, centery,radius):
        Y, X = np.ogrid[:h,:w]
        dist_from_center = np.sqrt((X-centerx)**2 + (Y-centery)**2)
        mask = dist_from_center <= radius
        return mask
###################################
try:
############################################################################
 if limit == 'y':
 
  	lim_spec = specs[lim_spec-1]
  	soc = soc[lim_evt-1]
	lim_evt = evt[lim_evt-1]

  	A_limit = explimit_file_A.split('.')[0]
	specA = dirlist+'/spec/spec'+det+'_'+lim_evt+'_'+lim_spec+'_explimby_'+ A_limit+'.pha'
	if os.path.isfile(specA) == True:
		print("Spectrum has already been done, either delete it or change the name!")
		sys.exit()
	os.system("cp "+dirlist+"/spec/specA_orig.pha "+specA)

	with fits.open(explimit_file_A) as hdulA:
                expA_limit = hdulA[1].header['EXPOSURE']
	with fits.open(specA) as hdulA:
                countsA = hdulA[1].data['COUNTS']
                expA = hdulA[1].header['EXPOSURE']

	q = 0

	f = open(dirlist+'/logs/'+obslist,'r+')
        for line in f.readlines():

                obsid = line[0:11]
                print(obsid)

		if obsid[0] == '7' or obsid[0] == '2':
                        print('Declined series OBS')
                        continue

                if len(line) == 0:
                        break
		if expA > expA_limit:
			break

		q += 1

		if lim_spec == 'full':
                        eventA = obsdir+'/'+obsid+'/event_cl/nu'+obsid+det+lim_evt+'_cl.evt'
                elif lim_spec == 'nosun':
                        eventA = obsdir+'/'+obsid+'/event_sep_cl/'+soc+'/nu'+obsid+det+lim_evt+'_fullevts_NOSUN.fits'
                elif lim_spec == 'sun':
                        eventA = obsdir+'/'+obsid+'/event_sep_cl/'+soc+'/nu'+obsid+det+lim_evt+'_fullevts_SUN.fits'

		if os.path.isfile(eventA) == True:
                        with fits.open(eventA) as hdulA:
                                PIA = hdulA[1].data['PI']
                                expA_s = hdulA[1].header['EXPOSURE']
                                RA_A = hdulA[0].header['RA_PNT']
                                DEC_A = hdulA[0].header['DEC_PNT']
			c = sc(ra = RA_A*u.degree, dec = DEC_A*u.degree, frame='fk5')
                        if abs(c.galactic.b.deg) >= BAD_PNT:
                                for k in range(len(PIA)):
                                        countsA[PIA[k]] += 1
                                expA = expA + expA_s
                                print("Done "+obsid+"A")
                        else:
                                print('To close to galactic equator')
	f.close()

	with fits.open(specA, mode = 'update') as hdulA:
                hdulA[1].data['COUNTS'] = countsA
                hdulA[0].header['EXPOSURE'] = expA
                hdulA[1].header['EXPOSURE'] = expA
                hdulA[2].header['EXPOSURE'] = expA
                hdulA[3].header['EXPOSURE'] = expA
                hdulA.flush()
	print("Total exposure for stack:   "+str(expA))
	print("Exposure limited value:     "+str(expA_limit))
	print("Detector:        "+det)
	print("Number of OBS:   "+str(q))

#############################################################################
 else:
  src_c = raw_input('Is there excl.reg files?  ["y","n"]:  ')
  sub_dir = raw_input('what is the sub directory?:     ')
  for i in i_choi:
#    for j in range(2):
    for j in [1]: # This is done to only stack the 02's from the list 
#    for j in [0]:
	  
	specA = dirlist+'/spec/'+sub_dir+'/specA_'+evt[j]+'_'+specs[i]+'.pha'
	specB = dirlist+'/spec/'+sub_dir+'/specB_'+evt[j]+'_'+specs[i]+'.pha'	
	specA0 = dirlist+'/spec/'+sub_dir+'/specA_'+evt[j]+'_'+specs[i]+'_det0.pha'
	specA1 = dirlist+'/spec/'+sub_dir+'/specA_'+evt[j]+'_'+specs[i]+'_det1.pha'
	specA2 = dirlist+'/spec/'+sub_dir+'/specA_'+evt[j]+'_'+specs[i]+'_det2.pha'
	specA3 = dirlist+'/spec/'+sub_dir+'/specA_'+evt[j]+'_'+specs[i]+'_det3.pha'
	specB0 = dirlist+'/spec/'+sub_dir+'/specB_'+evt[j]+'_'+specs[i]+'_det0.pha'	
	specB1 = dirlist+'/spec/'+sub_dir+'/specB_'+evt[j]+'_'+specs[i]+'_det1.pha'
	specB2 = dirlist+'/spec/'+sub_dir+'/specB_'+evt[j]+'_'+specs[i]+'_det2.pha'
	specB3 = dirlist+'/spec/'+sub_dir+'/specB_'+evt[j]+'_'+specs[i]+'_det3.pha'

	stackedlistA = dirlist+'/logs/specA'+evt[j]+specs[i]+sub_dir+'.dat'
	stackedlistB = dirlist+'/logs/specB'+evt[j]+specs[i]+sub_dir+'.dat'

	oldobsA, oldobsB = [],[]
	if os.path.isfile(stackedlistA) == True:
        	g = open(stackedlistA, 'r+')
        	for line in g.readlines():
               	  ob = line[0:11]
		  if ob == '###########':
		    print('next set')
		    continue
               	  if len(line) == 0:
                    break
                  oldobsA.append(ob)
        	g.close()

	if os.path.isfile(stackedlistB) == True:
                g = open(stackedlistB, 'r+')
                for line in g.readlines():
                  ob = line[0:11]
		  if ob == '###########':
		    continue
                  if len(line) == 0:
                    break
                  oldobsB.append(ob)
                g.close()

	with fits.open(specA) as hdulA:
		countsA = hdulA[1].data['COUNTS']
		expA = hdulA[1].header['EXPOSURE']
	with fits.open(specB) as hdulB:
		countsB = hdulB[1].data['COUNTS']
		expB = hdulB[1].header['EXPOSURE']
	with fits.open(specA0) as hdulA:
                countsA0 = hdulA[1].data['COUNTS']
                expA0 = hdulA[1].header['EXPOSURE']
        with fits.open(specA1) as hdulA:
                countsA1 = hdulA[1].data['COUNTS']
                expA1 = hdulA[1].header['EXPOSURE']
        with fits.open(specA2) as hdulA:
                countsA2 = hdulA[1].data['COUNTS']
                expA2 = hdulA[1].header['EXPOSURE']
        with fits.open(specA3) as hdulA:
                countsA3 = hdulA[1].data['COUNTS']
                expA3 = hdulA[1].header['EXPOSURE']	
	with fits.open(specB0) as hdulB:
                countsB0 = hdulB[1].data['COUNTS']
                expB0 = hdulB[1].header['EXPOSURE']
        with fits.open(specB1) as hdulB:
                countsB1 = hdulB[1].data['COUNTS']
                expB1 = hdulB[1].header['EXPOSURE']
        with fits.open(specB2) as hdulB:
                countsB2 = hdulB[1].data['COUNTS']
                expB2 = hdulB[1].header['EXPOSURE']
        with fits.open(specB3) as hdulB:
                countsB3 = hdulB[1].data['COUNTS']
                expB3 = hdulB[1].header['EXPOSURE']
	detAs = [countsA0,countsA1,countsA2,countsA3]
	detBs = [countsB0,countsB1,countsB2,countsB3]
	print(expA) #########
	print(expB) #########

##################################################
#  The stack  #
##### added image check for counts:

	f = open(dirlist+'/text_files/'+obslist,'r+')
	for line in f.readlines():
	
		obsid = line[0:11]
		print(obsid)
# The next lines are to see if the obs is a cluster (they will not 
# help the background measurement, and will most likely skew it)
		if j == 0:
		    if obsid[0] == '7' or obsid[0] == '2':
			print('Declined series OBS')
			continue

		if j == 1:
		    BAD_PNT = 0.0

		if len(line) == 0:
			break
		newA = 'True'
		newB = 'True'
	#newlist = oldobs + obsid
		
		if i == 0:
			eventA = obsdir+'/'+obsid+'/event_cl/nu'+obsid+'A'+evt[j]+'_cl.evt'
			eventB = obsdir+'/'+obsid+'/event_cl/nu'+obsid+'B'+evt[j]+'_cl.evt'
			#if len(set(oldobsA)&set([obsid])) > 0:
			if oldobsA.count(obsid) > 0:
			#if len(oldobsA) == len(set(oldobsA + [obsid])):
                        	print("duplicate obsid:  "+str(obsid))
                        	newA = 'False'
                	#if len(set(oldobsB)&set([obsid])) > 0:
			if oldobsB.count(obsid) > 0:
			#if len(oldobsB) == len(set(oldobsB + [obsid])):
                        	print("duplicate obsid:  "+str(obsid))
                        	newB = 'False'
			print('full events')
                        print('{0} {1}'.format(i, j))
                        print('{0} {1}'.format(newA, newB))
		elif i == 1:
			eventA = obsdir+'/'+obsid+'/event_sep_cl/'+soc[j]+'/nu'+obsid+'A'+evt[j]+'_fullevts_NOSUN.fits'
			eventB = obsdir+'/'+obsid+'/event_sep_cl/'+soc[j]+'/nu'+obsid+'B'+evt[j]+'_fullevts_NOSUN.fits'
			print('no sun events')
                        if len(set(oldobsA)&set([obsid])) > 0:
		     	#if len(oldobsA) == len(set(oldobsA + [obsid])):
                                print("duplicate obsid:  "+str(obsid))
                                newA = 'False'
			if len(set(oldobsB)&set([obsid])) > 0:
                        #if len(oldobsB) == len(set(oldobsB + [obsid])):
                                print("duplicate obsid:  "+str(obsid))
                                newB = 'False'
                        print('{0} {1}'.format(i, j))
                        print('{0} {1}'.format(newA, newB))
		elif i == 2:
			eventA = obsdir+'/'+obsid+'/event_sep_cl/'+soc[j]+'/nu'+obsid+'A'+evt[j]+'_fullevts_SUN.fits'
			eventB = obsdir+'/'+obsid+'/event_sep_cl/'+soc[j]+'/nu'+obsid+'B'+evt[j]+'_fullevts_SUN.fits'
			print('sun events')
			if len(set(oldobsA)&set([obsid])) > 0:
                        #if len(oldobsA) == len(set(oldobsA + [obsid])):
                                print("duplicate obsid:  "+str(obsid))
                                newA = 'False'
			if len(set(oldobsB)&set([obsid])) > 0:
                        #if len(oldobsB) == len(set(oldobsB + [obsid])):
                                print("duplicate obsid:  "+str(obsid))
                                newB = 'False'
                        print('{0} {1}'.format(i, j))
                        print('{0} {1}'.format(newA, newB))

		print('{0} {1}'.format(eventA,eventB))
# Here we will create a mask around the source if it exists.
		if src_c == 'y':
			print('doing src')
			src_dir = obsdir+'/'+obsid+'/event_defcl/excl.reg'
			if os.path.isfile(src_dir) == True:
				srcfile = open(src_dir,'r+')
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
			else:
				grid_pi = np.ones((1000,1000))
		if os.path.isfile(eventA) and newA == 'True':
			with fits.open(eventA) as hdulA:
                               	PIA = hdulA[1].data['PI']
                               	expA_s = hdulA[1].header['EXPOSURE']
				detidA = hdulA[1].data['DET_ID']
				RA_A = hdulA[0].header['RA_PNT']
				DEC_A = hdulA[0].header['DEC_PNT']
				X_A = hdulA[1].data['X']
				Y_A = hdulA[1].data['Y']
				DET1_XA = hdulA[1].data['DET1X']
				DET1_YA = hdulA[1].data['DET1Y']
			print("A")

# The following lines are to see from the recorded pointing if the background would be contaminated by
# the GRXE, This is repeated for B because there is no guarantee that
# both event files exist.  Yes I could save lines of code by being clever,
# however without rewriting and until I get a better understanding of how
# to use A or B kind of code, this is the way I went.

			c = sc(ra = RA_A*u.degree, dec = DEC_A*u.degree, frame='fk5')
			if abs(c.galactic.b.deg) >= BAD_PNT:
				for k in range(len(PIA)):
				    if src_c == 'y':
				      #if X_A[k] == -1:  # From here, i replace the if with (edgemaskB[DET1_YB[k]-1,DET1_XB[k]-1] == 0): added *4th*
				      if (X_A[k] == -1) or (edgemaskA[DET1_YA[k]-1,DET1_XA[k]-1] == 0):
					continue
				      else:
					countsA[PIA[k]-1] += int(grid_pi[Y_A[k]-1,X_A[k]-1])
					detAs[detidA[k]][PIA[k]-1] += int(grid_pi[Y_A[k]-1,X_A[k]-1])	
				    else:
				      if (X_A[k] == -1) or (edgemaskA[DET1_YA[k]-1,DET1_XA[k]-1] == 0):
				        continue
				      else:
                                	countsA[PIA[k]-1] += 1
					detAs[detidA[k]][PIA[k]-1] += 1
					###### Add image here
					if img_choi == 'Y':
						if (PIA[k] >= el_choi) and (PIA[k] <= hl_choi):
							imgApack[i][DET1_YA[k]-1,DET1_XA[k]-1] += 1
					####################
				expA = expA + expA_s
                        	print("Done "+obsid+"A")
			else:
				print('To close to galactic equator')
				newA = 'False'
		#else:
			#newA = str(os.path.isfile(eventA))

# Should be noted that the test for GRXE is put before the counts addition, and 
# that if either A or B point in the direction of BAD_PNT or less the code goes
# to the next loop iteration under the assumption that A and B always point in 
# the same direction to a nominal degree.

		if os.path.isfile(eventB) and newB == 'True':
    			with fits.open(eventB) as hdulB:
                               	PIB = hdulB[1].data['PI']
                               	expB_s = hdulB[1].header['EXPOSURE']
				detidB = hdulB[1].data['DET_ID']
				RA_B = hdulB[0].header['RA_PNT']
				DEC_B = hdulB[0].header['DEC_PNT'] 
				X_B = hdulB[1].data['X']
				Y_B = hdulB[1].data['Y']
				DET1_XB = hdulB[1].data['DET1X']
				DET1_YB = hdulB[1].data['DET1Y']
			print("B")      

			c = sc(ra = RA_B*u.degree, dec = DEC_B*u.degree, frame='fk5')
			if abs(c.galactic.b.deg) >= BAD_PNT:                   
				for k in range(len(PIB)):
				    if src_c == 'y':
                                      if (X_B[k] == -1) or (edgemaskB[DET1_YB[k]-1,DET1_XB[k]-1] == 0):
                                        continue
                                      else:
                                        countsB[PIB[k]-1] += int(grid_pi[Y_B[k]-1,X_B[k]-1])
                                        detBs[detidB[k]][PIB[k]-1] += int(grid_pi[Y_B[k]-1,X_B[k]-1])
                                    else:
				      if (X_B[k] == -1) or (edgemaskB[DET1_YB[k]-1,DET1_XB[k]-1] == 0):
					continue
				      else:
                               		countsB[PIB[k]-1] += 1
                                        detBs[detidB[k]][PIB[k]-1] += 1
					####### Add image here
					if img_choi == 'Y':
						if (PIB[k] >= el_choi) and (PIB[k] <= hl_choi):
							imgBpack[i][DET1_YB[k]-1,DET1_XB[k]-1] += 1
					######################
				expB = expB + expB_s
                        	print("Done "+obsid+"B")
			else:
				print('To close to galactic equator')
				newB = 'False'
		#else:
			#newB = str(os.path.isfile(eventB))

		if newA == 'True':
		  with open(stackedlistA, 'a+') as blah:
			blah.write(obsid+'\n')
		if newB == 'True':
		  with open(stackedlistB, 'a+') as blah:
                        blah.write(obsid+'\n')
	f.close()
	print(expA) #########
	print(expB) #########

	# This part is to rewrite the spectrum into the file
	with fits.open(specA, mode = 'update') as hdulA:
		hdulA[1].data['COUNTS'] = countsA
		hdulA[0].header['EXPOSURE'] = expA
		hdulA[1].header['EXPOSURE'] = expA
		hdulA[2].header['EXPOSURE'] = expA
		hdulA[3].header['EXPOSURE'] = expA
		hdulA.flush()

	with fits.open(specB, mode = 'update') as hdulB: 
		hdulB[1].data['COUNTS'] = countsB
		hdulB[0].header['EXPOSURE'] = expB
		hdulB[1].header['EXPOSURE'] = expB
		hdulB[2].header['EXPOSURE'] = expB
		hdulB[3].header['EXPOSURE'] = expB
		hdulB.flush()
	
	fitA = [specA0,specA1,specA2,specA3]
	fitB = [specB0,specB1,specB2,specB3]
	for w in range(len(fitA)):
		with fits.open(fitA[w], mode = 'update') as hdulA:
			hdulA[1].data['COUNTS'] = detAs[w]
                	hdulA[0].header['EXPOSURE'] = expA
                	hdulA[1].header['EXPOSURE'] = expA
                	hdulA[2].header['EXPOSURE'] = expA
                	hdulA[3].header['EXPOSURE'] = expA
                	hdulA.flush()

        for w in range(len(fitA)):
                with fits.open(fitB[w], mode = 'update') as hdulB:
                        hdulB[1].data['COUNTS'] = detBs[w]
                        hdulB[0].header['EXPOSURE'] = expB
                        hdulB[1].header['EXPOSURE'] = expB
                        hdulB[2].header['EXPOSURE'] = expB
                        hdulB[3].header['EXPOSURE'] = expB
                        hdulB.flush()

	if img_choi == 'Y':	
		TestA = fits.PrimaryHDU(imgApack[i])
		TestB = fits.PrimaryHDU(imgBpack[i])

		TestA.writeto(dirlist+'/spec/'+sub_dir+'/'+imgAnames[i]+'_'+evt[j]+'_'+str(int(el_c))+'_'+str(int(hl_c))+'keV.fits')
		TestB.writeto(dirlist+'/spec/'+sub_dir+'/'+imgBnames[i]+'_'+evt[j]+'_'+str(int(el_c))+'_'+str(int(hl_c))+'keV.fits')

except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(str(exc_type) + '\n' + str(e) +'\n' +str(exc_tb.tb_lineno))

###################################################
#if __name__ == "__main__":
#	main()
