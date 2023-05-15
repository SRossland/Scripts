#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./excl_02.py

# Created by: Steven P. Rossland, 2018

# This script is meant to be ran from the ../xray/NuSTAR/ directory
# The list given is read in from /uufs/astro05.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs
#  

import sys, shutil, os, string, time, datetime, pidly, numpy as np
from photutils import find_peaks
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
    
def main():
    dirhome = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs'
# Give working dir
    dirnow = os.getcwd()  #
    print("Which sratch drive:")
    print("  [1] astro05 local ")
    print("  [2] general/lustre")
    print("  [3] kingspeak/serial")
    print("  [4] lonepeak/serial")

    choice = eval(raw_input("    Choose 1-4:"))
# Input server
    while choice not in [1,2,3,4]:
        print("That is not a choice, program exiting")
	sys.exit()
    if choice == 1:
	waiting = raw_input("Are you logged into astro05?:  [y/n]")
        if waiting == 'n':
            print("Please log into astro05")
            sys.exit()
    elif choice != 1:
	waiting = raw_input("Are you logged into notchpeak?:  [y/n]")
        if waiting == 'n':
            print("Please log into notchpeak")
	    sys.exit()
    if choice == 1:
    	scratch = '/scratch/local'
    elif choice == 2:
        scratch = '/scratch/general/lustre'
    elif choice == 3:
        scratch = '/scratch/kingspeak/serial'
    elif choice == 4:
   	scratch = '/scratch/lonepeak/serial'

    dir = dirnow+'/OBS'
    wlist = raw_input("What download list does this correspond to? (found in ~/NuSTAR/Scripts/logs (=/= '*.dat'):  ")
	# Give obsid list
    mlist, crap = wlist.split(".")
    
    os.system("mkdir -p "+scratch+"/NuSTAR/OBS")

    # Here I am trying to catch if a list has been run before:

    if os.path.isfile(dirnow+'/LogFile/Running_02.txt') == True:
	obs_test = []
	file_obs_test = dirnow+'/LogFile/Running_02.txt'
	f = open(file_obs_test,'r+')
	for line in f.readlines():
	    l = line.rstrip()
	    obs_test.append(l)
	f.close()
	obs_test.append(mlist)
	if len(obs_test) != len(set(obs_test)) == True:
	    print('Duplicate List given: Check running list')
	    sys.exit()    

    with open(dirnow+'/LogFile/Running_02.txt','a+') as run:
	run.write(mlist+'\n')    

    # Open the file to get the values of the background for both A and B

#    fbgd = np.loadtxt(dirnow+'/LogFile/avg_level_3_30keV.dat') #
#    A_norm_bgd = fbgd[0]
#    B_norm_bgd = fbgd[1]
    
    
#    fPSF = np.loadtxt(dirnow+'/LogFile/prac_PSF.dat') #
#    fitx = fPSF[:,0]
#    fity = fPSF[:,1]

    f = open(dirhome+'/'+mlist+'.dat','r+') #
        
    for line in f.readlines():
        try:
            obsid = line[0:11]
            
	    if len(line) == 0:
                break
# This section will erase the defcl folder and set the progress back
	#    if os.path.isdir(dir+'/'+obsid+'/event_defcl') == True:
        #        if os.path.isdir(dir+'/'+obsid+'/event_cl') == True:
        #            os.system("rm -r -f "+dir+'/'+obsid+'/event_cl')
        #        os.system("mv "+dir+'/'+obsid+"/event_defcl "+dir+'/'+obsid+"/event_cl")


            with open(dirnow+'/LogFile/Errors_02.txt','a+') as E:
		E.write('\n'+obsid+'\n')

 	    with open(dirnow+'/LogFile/ObsProc'+mlist+'.txt','a+') as O:
		O.write(obsid+'\n') 	    

# From here all files will be moved to the scratch drive and processed
	    
	    os.chdir(scratch)
	    os.system("mkdir -p NuSTAR/OBS")
	    os.chdir("NuSTAR")
	    os.system("cp -r -p "+dir+"/"+obsid+" ./OBS")
	    dirs = scratch+'/NuSTAR'
	    dirs = dirs+'/OBS'
	    
	    os.system("gunzip -d -f "+dirs+"/"+obsid+"/event_uf/*.gz")

	    wfile = dirs+'/'+obsid+'/Proc_log.txt'
    
            sttime = str(datetime.datetime.now())
            with open(wfile,'a+') as obsfile:
            	obsfile.write(obsid +"--- Start: "+ sttime + '\n')

           # excltrue = (['y','n'])
           # etrue = excltrue[1]
            cldir = dirs+'/'+obsid+'/event_defcl'

            if os.path.isfile(cldir+'/nu'+obsid+'A02_cl.evt') and os.path.isfile(cldir+'/nu'+obsid+'B02_cl.evt') == True:
	        ab = ['A','B']
	    elif os.path.isfile(cldir+'/nu'+obsid+'A02_cl.evt') == True:
		ab = 'A'
	    elif os.path.isfile(cldir+'/nu'+obsid+'B02_cl.evt') == True:
		ab = 'B'
	    else:
		continue  

	   # os.system("gunzip -d "+cldir+"/*.gz")
 
            idl = pidly.IDL('/uufs/chpc.utah.edu/sys/pkg/idl/8.4/idl84/bin/idl') #

           # if os.path.isfile(cldir+'/excl.reg'):
               # os.system("rm -f "+cldir+'/excl.reg')

            idl("cd, current=dir")
	    idl("dir=dir+'/OBS/'")
            idl("obsid='"+obsid+"'")
            idl("imname='im3to30keV_02.fits'")
            idl("cldir = '"+cldir+"/'")
	    for i in range(len(ab)):
		idl("mkimgs_02,cldir,obsid,'"+ab[i]+"',3,30")
		    #idl("mkimgs_02,cldir,obsid,'B',3,30")
	    	#os.system("mv "+cldir+"/im"+ab[i]+"3to30keV.fits "+cldir+"/im"+ab[i]+"3t030keV_02.fits")
            	#os.system("mv "+cldir+"/imB3to30keV.fits imB3to30keV_02.fits")
	    if len(ab) == 2:
            	idl("fits_read,cldir+'imA3to30keV_02.fits',im1,h")
            	idl("fits_read,cldir+'imB3to30keV_02.fits',im2")
            	idl("im=im1+im2")
	    else:
		idl("fits_read,cldir+'im"+ab[0]+"3to30keV_02.fits',im1,h")
		idl("im=im1")
            idl("fits_write,cldir+imname,im,h")
            idl("tbin = 100")
	    
	    with open(wfile,'a+') as obsfile:
		obsfile.write("     A, B, combined image created for occ"+'\n')            
     # All the following is for excl.reg finding
#################################################################
	   # tbl, exp = peaks(cldir+'/im3to30keV.fits',A_norm_bgd)
           # tblcopy = tbl.copy()
            

           # A_bgd = A_norm_bgd * exp

           # if len(tbl) >= 1:
                
	#	with open(wfile, 'a+') as obsfile:
	#	    obsfile.write("	excl.reg attempt"+'\n')

        #        srcarray = createSrc(tblcopy, fitx, fity, A_bgd, wfile)
        #        etrue = excltrue[0]
        #        srcarray = np.reshape(srcarray,(len(srcarray)/3,3))
                    
        #    else:
        #        etrue = excltrue[1]
                #os.system("rm -f "+cldir+'excl.reg')
                        #Need to erase all values in the tbl associated with mask



        #    if etrue == 'y':
        #        with open(cldir+'/excl.reg','a+') as o:
        #            o.write('# Region file format: DS9 version 4.1\n')
        #            o.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n physical\n')
        #        for i in range(len(srcarray)):
        #            with open(cldir+'/excl.reg','a+') as o:
        #                o.write("circle("+str(srcarray[i,0])+","+str(srcarray[i,1])+","+str(srcarray[i,2])+")"+'\n')
	#	with open(wfile,'a+') as obsfile:
	#	    obsfile.write("    excl.reg file created"+'\n')
##################################################################
            for i in range(len(ab)):
	    	idl("lcfilter_02, dir, obsid, '"+ab[i]+"', 50, 150, imname, tbin")
            	idl("lcfilter_02, dir, obsid, '"+ab[i]+"', 50, 150, imname, tbin, /usr")
           # if etrue == 'y':
           #     idl("lcfilter, dir, obsid, 'A', 3, 20, imname, tbin, /usr, /excl")
	   # else:
	    	idl("lcfilter_02, dir, obsid, '"+ab[i]+"', 3, 7, imname, tbin, /usr")

            #idl("lcfilter_02, dir, obsid, 'B', 50, 160, imname, tbin")
            #idl("lcfilter_02, dir, obsid, 'B', 50, 160, imname, tbin, /usr")
           # if etrue == 'y':
           #     idl("lcfilter, dir, obsid, 'B', 3, 20, imname, tbin, /usr, /excl")
	   # else:
	    #idl("lcfilter_02, dir, obsid, 'B', 3, 20, imname, tbin, /usr")

	    with open(wfile,'a+') as obsfile:
  		obsfile.write("    GTI_02 file created"+'\n')
	    for dets in ['A','B']:
            	os.system("mv "+cldir+"/nu"+obsid+dets+"02_usrgti.fits "+dirs+"/"+obsid+"/")
           # os.system("mv "+cldir+" "+dir+"/"+obsid+"/event_defcl")
	    for i in range(len(ab)):
		os.system("source /uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/xray/caldb/heasoft/CALDB/software/tools/caldbinit.sh")
	    	os.system("/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/Scripts/run_pipe_usrgti_notstrict_02.sh "+dirs+"/"+obsid+" "+ab[i]) #
            #os.system("/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/Scripts/run_pipe_usrgti_notstrict.sh "+dir+"/"+obsid+" B") #

	    with open(wfile, 'a+') as obsfile:
		obsfile.write("    reprocessed w/GTI_02 file"+'\n')

# Now everything should be in the event_cl file

            cldir = dirs+'/'+obsid+'/event_cl'

#	    idl("cd, current=dir")
            idl("obsid='"+obsid+"'")
            idl("imname='im3to30keV_02.fits'")
            idl("cldir = '"+cldir+"/'")
            for i in range(len(ab)):
                idl("mkimgs_02,cldir,obsid,'"+ab[i]+"',3,30")
                #idl("mkimgs_02,cldir,obsid,'B',3,30")
                #os.system("mv "+cldir+"/im"+ab[i]+"3to30keV.fits "+cldir+"/im"+ab[i]+"3t030keV_02.fits")
                #os.system("mv "+cldir+"/imB3to30keV.fits imB3to30keV_02.fits")
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

            for j in ['A','B']:
	        os.system("sepsunocc2.py "+dirs+" "+obsid+" "+j)

	    endtime = str(datetime.datetime.now()) 
	    with open(wfile, 'a+') as obsfile:
		obsfile.write("    FINISHED:  "+endtime+'\n')

# The next 3 lines are to copy only the changed files back over to the origninal source    	    
	    os.system("cp -r -u ./OBS/"+obsid+" "+dir)
	    os.chdir(scratch)
	    shutil.rmtree(scratch+"/NuSTAR/OBS/"+obsid)
 	    os.chdir(dirnow)
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
	    with open(wfile,'a+') as obsfile:
                   obsfile.write('\n'+"="*41+" ERRORS "+"="*41+ '\n')
                   obsfile.write(str(exc_type) + '\n' + str(e) + '\n' + str(exc_tb.tb_lineno))
                   obsfile.write('\n' + "="*90 + '\n')
	    continue
    print("EXCL_02.PY COMPLETE")
    os.system("date")
    f.close()

"""Need to convert my writing of errors and such to logging"""
if __name__ == "__main__":
    main()

