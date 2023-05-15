#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./excl_02.py

# Created by: Steven P. Rossland, 2018

# This script is meant to be ran from the ../xray/NuSTAR/ directory
# you HAVE to run this from the above directory, xspec does not work well when you don't
# The list given is read in from /uufs/astro05.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs
#
# THIS IS THE BATCH VERSION, IT IS WRITTEN TO RUN ON THE KINGSPEAK SERVER AND WILL
# SUBMIT A JOB TO RUN ON THE KINGSPEAK SCRATCH DRIVE  

import sys, shutil, os, string, time, datetime, pidly, numpy as np
    
def main():
    dirhome = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs'
# Give working dir
    dirnow = '/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/xray/NuSTAR'  #it was os.getcwd(), however I changed it to the nustar/obs dir in the XRAY directory
    scratch = '/scratch/kingspeak/serial/u1019304/NuSTAR'

    dir = dirnow+'/OBS'
#####
# THE OBSIDS ARE GIVEN IN THE SYS ARGV ARGUMENT OF THE BATCH SUBMISSION
    obsid = sys.argv[1][0:11]
    if len(obsid) != 11:
	print('Not a valid obsid length')
	sys.exit()

    # Here I am trying to catch if a list has been run before:
#############I need to rewrite this to use the obs.txt file ####
#    if os.path.isfile(dirnow+'/LogFile/Running_02.txt') == True:
#	obs_test = []
#	file_obs_test = dirnow+'/LogFile/Running_02.txt'
#	f = open(file_obs_test,'r+')
#	for line in f.readlines():
#	    l = line.rstrip()
#	    obs_test.append(l)
#	f.close()
#	obs_test.append(mlist)
#	if obs_test.count(obsid) > 0:
#	    print('Duplicate List given: Check running list')
#	    sys.exit()    
######################
    try:

# From here all files will be moved to the scratch drive and processed
	    
	    os.chdir(scratch)
#	    os.system("mkdir -p NuSTAR/OBS")
#	    os.chdir("NuSTAR")
#	    os.system("cp -r -p "+dir+"/"+obsid+" ./OBS")
	    dirs = scratch+'/OBS'
#	    dirs = dirs+'/OBS'
	    
	    os.system("gunzip -d -f "+dirs+"/"+obsid+"/event_uf/*.gz")

            sttime = str(datetime.datetime.now())
            print(obsid +"--- Start: "+ sttime + '\n')

            cldir = dirs+'/'+obsid+'/event_defcl'

            if os.path.isfile(cldir+'/nu'+obsid+'A02_cl.evt') and os.path.isfile(cldir+'/nu'+obsid+'B02_cl.evt') == True:
	        ab = ['A','B']
	    elif os.path.isfile(cldir+'/nu'+obsid+'A02_cl.evt') == True:
		ab = 'A'
	    elif os.path.isfile(cldir+'/nu'+obsid+'B02_cl.evt') == True:
		ab = 'B'
	    else:
		sys.exit() #continue  

 
            idl = pidly.IDL('/uufs/chpc.utah.edu/sys/pkg/idl/8.4/idl84/bin/idl') #


            idl("cd, current=dir")
	    idl("dir=dir+'/OBS/'")
            idl("obsid='"+obsid+"'")
            idl("imname='im3to30keV_02.fits'")
            idl("cldir = '"+cldir+"/'")
	    for i in range(len(ab)):
		idl("mkimgs_02,cldir,obsid,'"+ab[i]+"',3,30")
	    if len(ab) == 2:
            	idl("fits_read,cldir+'imA3to30keV_02.fits',im1,h")
            	idl("fits_read,cldir+'imB3to30keV_02.fits',im2")
            	idl("im=im1+im2")
	    else:
		idl("fits_read,cldir+'im"+ab[0]+"3to30keV_02.fits',im1,h")
		idl("im=im1")
            idl("fits_write,cldir+imname,im,h")
            idl("tbin = 100")
	    
##################################################################
            for i in range(len(ab)):
	    	idl("lcfilter_02, dir, obsid, '"+ab[i]+"', 50, 150, imname, tbin")
            	idl("lcfilter_02, dir, obsid, '"+ab[i]+"', 50, 150, imname, tbin, /usr")
	    	idl("lcfilter_02, dir, obsid, '"+ab[i]+"', 3, 7, imname, tbin, /usr")

	    for dets in ['A','B']:
            	os.system("mv "+cldir+"/nu"+obsid+dets+"02_usrgti.fits "+dirs+"/"+obsid+"/")

	    for i in range(len(ab)):
	    	os.system("/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/Scripts/run_pipe_usrgti_notstrict_02.sh "+dirs+"/"+obsid+" "+ab[i]) 


            cldir = dirs+'/'+obsid+'/event_cl'

 	    idl("cd, current=dir")
            idl("obsid='"+obsid+"'")
            idl("imname='im3to30keV_02.fits'")
            idl("cldir = '"+cldir+"/'")
            for i in range(len(ab)):
                idl("mkimgs_02,cldir,obsid,'"+ab[i]+"',3,30")
            if len(ab) == 2:
                idl("fits_read,cldir+'imA3to30keV_02.fits',im1,h")
                idl("fits_read,cldir+'imB3to30keV_02.fits',im2")
                idl("im=im1+im2")
            else:
                idl("fits_read,cldir+'im"+ab[0]+"3to30keV_02.fits',im1,h")
                idl("im=im1")
            idl("fits_write,cldir+imname,im,h")
            idl.close()

####Insert here, elevcorr.py


	    endtime = str(datetime.datetime.now()) 
	    print("    FINISHED:  "+endtime+'\n')
	    os.system('mv '+dirhome+'/OBS_ID/'+obsid+'.obs '+dirhome+'/OBS_ID_DONE/')

# The next 3 lines are to copy only the changed files back over to the origninal source    	    
#	    os.system("cp -r -u ./OBS/"+obsid+" "+dir)
#	    os.chdir(scratch)
#	    shutil.rmtree(scratch+"/NuSTAR/OBS/"+obsid)
#	    os.chdir(dirnow)
    except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
	    #with open(wfile,'a+') as obsfile:
            print('\n'+"="*41+" ERRORS "+"="*41+ '\n')
            print(str(exc_type) + '\n' + str(e) + '\n' + str(exc_tb.tb_lineno))
            print('\n' + "="*90 + '\n')
	    sys.exit()  # To run as a batch, a sys.exit() is needed over the continue it had

"""Need to convert my writing of errors and such to logging"""
if __name__ == "__main__":
    main()

