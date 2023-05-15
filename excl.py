#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./excl.py

# Created by: Steven P. Rossland, 2018

# This script is meant to be ran from the ../xray/NuSTAR/ directory

#  

import sys, shutil, os, string, time, datetime, pidly, numpy as np
from photutils import find_peaks
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
    
def main():
    dirhome = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs'
    dirnow = os.getcwd()  #
    
    print("Which sratch drive:")
    print("  [1] astro05 local ")
    print("  [2] general/lustre")
    print("  [3] kingspeak/serial")
    print("  [4] lonepeak/serial")

    choice = eval(raw_input("    Choose 1-4:"))

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
    
    wlist = raw_input("What download list does this correspond to? (local path from current location with no file extension (=/= '*.dat'):  ")
    mlist, crap = wlist.split(".")

    os.system("mkdir -p "+scratch+"/NuSTAR/OBS")

    # Here I am trying to catch if a list has been run before:

    if os.path.isfile(dirnow+'/LogFile/Running.txt') == True:
	obs_test = []
	file_obs_test = dirnow+'/LogFile/Running.txt'
	f = open(file_obs_test,'r+')
	for line in f.readlines():
	    l = line.rstrip()
	    obs_test.append(l)
	f.close()
	obs_test.append(mlist)
	if len(obs_test) != len(set(obs_test)) == True:
	    print('Duplicate List given: Check running list')
	    sys.exit()    

    with open(dirnow+'/LogFile/Running.txt','a+') as run:
	run.write(mlist+'\n')    

    # Open the file to get the values of the background for both A and B

    fbgd = np.loadtxt(dirnow+'/LogFile/avg_level_3_30keV.dat') #
    A_norm_bgd = fbgd[0]
    B_norm_bgd = fbgd[1]
    
    
    fPSF = np.loadtxt(dirnow+'/LogFile/prac_PSF.dat') #
    fitx = fPSF[:,0]
    fity = fPSF[:,1]

    f = open(dirhome+'/'+mlist+'.dat','r+') #
        
    for line in f.readlines():
        try:
            obsid = line[0:11]
            
	    if len(line) == 0:
                break

	    if os.path.isdir(dir+'/'+obsid+'/event_defcl') == True:
                if os.path.isdir(dir+'/'+obsid+'/event_cl') == True:
                    os.system("rm -r -f "+dir+'/'+obsid+'/event_cl')
               # os.system("mv "+dir+'/'+obsid+"/event_defcl "+dir+'/'+obsid+"/event_cl")


            with open(dirnow+'/LogFile/Errors.txt','a+') as E:
		E.write('\n'+obsid+'\n')

 	    with open(dirnow+'/LogFile/ObsProc'+mlist+'.txt','a+') as O:
		O.write(obsid+'\n') 	    

            os.chdir(scratch)
            os.system("mkdir -p NuSTAR/OBS")
            os.chdir("NuSTAR")
	    os.system("source /uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/xray/caldb/heasoft/CALDB/software/tools/caldbinit.sh")
            os.system("cp -r -p "+dir+"/"+obsid+" ./OBS")
            dirs = scratch+'/NuSTAR'
            dirs = dirs+'/OBS'

# This line is to move the scratch defcl file to the cl version to redo the processing
	    if os.path.isdir(dirs+'/'+obsid+'/event_defcl') == True:
	    	os.system("mv "+dirs+"/"+obsid+"/event_defcl "+dirs+"/"+obsid+"/event_cl")

            os.system("gunzip -d -f "+dirs+"/"+obsid+"/event_uf/*.gz")
 
	    os.system("source /uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/xray/caldb/heasoft/CALDB/software/tools/caldbinit.sh")

	    wfile = dirs+'/'+obsid+'/Proc_log.txt'
    
            sttime = str(datetime.datetime.now())
            with open(wfile,'a+') as obsfile:
            	obsfile.write(obsid +"--- Start: "+ sttime + '\n')

            excltrue = (['y','n'])
            etrue = excltrue[1]
            cldir = dirs+'/'+obsid+'/event_cl'
            
	    os.system("gunzip -d -f "+cldir+"/*.gz")
 
            idl = pidly.IDL('/uufs/chpc.utah.edu/sys/pkg/idl/8.4/idl84/bin/idl') #

            if os.path.isfile(cldir+'/excl.reg'):
                os.system("rm -f "+cldir+'/excl.reg')

            idl("cd, current=dir")
	    idl("dir=dir+'/OBS/'")
            idl("obsid='"+obsid+"'")
            idl("imname='im3to30keV.fits'")
            idl("cldir = '"+cldir+"/'")
            idl("mkimgs,cldir,obsid,'A',3,30")
            idl("mkimgs,cldir,obsid,'B',3,30")
            idl("fits_read,cldir+'imA3to30keV.fits',im1,h")
            idl("fits_read,cldir+'imB3to30keV.fits',im2")
            idl("im=im1+im2")
            idl("fits_write,cldir+imname,im,h")
            idl("tbin = 100")
	    
	    with open(wfile,'a+') as obsfile:
		obsfile.write("     A, B, combined image created"+'\n')            
            tbl, exp = peaks(cldir+'/im3to30keV.fits',A_norm_bgd)
            tblcopy = tbl.copy()
            

            A_bgd = A_norm_bgd * exp

            if len(tbl) >= 1:
                
		with open(wfile, 'a+') as obsfile:
		    obsfile.write("	excl.reg attempt"+'\n')

                srcarray = createSrc(tblcopy, fitx, fity, A_bgd, wfile)
                etrue = excltrue[0]
                srcarray = np.reshape(srcarray,(len(srcarray)/3,3))
                    
            else:
                etrue = excltrue[1]
                #os.system("rm -f "+cldir+'excl.reg')
                        #Need to erase all values in the tbl associated with mask



            if etrue == 'y':
                with open(cldir+'/excl.reg','a+') as o:
                    o.write('# Region file format: DS9 version 4.1\n')
                    o.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n physical\n')
                for i in range(len(srcarray)):
                    with open(cldir+'/excl.reg','a+') as o:
                        o.write("circle("+str(srcarray[i,0])+","+str(srcarray[i,1])+","+str(srcarray[i,2])+")"+'\n')
		with open(wfile,'a+') as obsfile:
		    obsfile.write("    excl.reg file created"+'\n')

            idl("lcfilter, dir, obsid, 'A', 50, 150, imname, tbin")
            idl("lcfilter, dir, obsid, 'A', 50, 150, imname, tbin, /usr")
            if etrue == 'y':
                idl("lcfilter, dir, obsid, 'A', 3, 7, imname, tbin, /usr, /excl")
	    else:
		idl("lcfilter, dir, obsid, 'A', 3, 7, imname, tbin, /usr")

            idl("lcfilter, dir, obsid, 'B', 50, 150, imname, tbin")
            idl("lcfilter, dir, obsid, 'B', 50, 150, imname, tbin, /usr")
            if etrue == 'y':
                idl("lcfilter, dir, obsid, 'B', 3, 7, imname, tbin, /usr, /excl")
	    else:
		idl("lcfilter, dir, obsid, 'B', 3, 7, imname, tbin, /usr")

	    with open(wfile,'a+') as obsfile:
  		obsfile.write("    GTI file created"+'\n')

            os.system("mv "+cldir+"/nu"+obsid+"*01_usrgti.fits "+dirs+"/"+obsid+"/")
            os.system("mv "+cldir+" "+dirs+"/"+obsid+"/event_defcl")
	    os.system("source /uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/xray/caldb/heasoft/CALDB/software/tools/caldbinit.sh")
            os.system("/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/Scripts/run_pipe_usrgti_notstrict.sh "+dirs+"/"+obsid+" A") #
            os.system("/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/Scripts/run_pipe_usrgti_notstrict.sh "+dirs+"/"+obsid+" B") #

	    with open(wfile, 'a+') as obsfile:
		obsfile.write("    reprocessed w/GTI file"+'\n')
           
	    idl("cd, current=dir")
            idl("obsid='"+obsid+"'")
            idl("imname='im3to30keV.fits'")
            idl("cldir = '"+cldir+"/'")
            idl("mkimgs,cldir,obsid,'A',3,30")
            idl("mkimgs,cldir,obsid,'B',3,30")
            idl("fits_read,cldir+'imA3to30keV.fits',im1,h")
            idl("fits_read,cldir+'imB3to30keV.fits',im2")
            idl("im=im1+im2")
            idl("fits_write,cldir+imname,im,h")
            idl.close()
           
	    endtime = str(datetime.datetime.now()) 
	    with open(wfile, 'a+') as obsfile:
		obsfile.write("    FINISHED:  "+endtime+'\n')

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
    print("EXCL.PY COMPLETE")
    f.close()

#########################################################################################################
def peaks(file, A_norm_bgd):
    imgar = fits.open(file)
    dataar = imgar[0].data
    hdr = imgar[0].header
    exp = hdr['Exposure']
    A_bgd = A_norm_bgd*exp
    gauss_kernel = Gaussian2DKernel(stddev = 6)
    dataarray = convolve(dataar,gauss_kernel)
    threshold = A_bgd*5
    tbl = find_peaks(dataarray, threshold, box_size=5.0)
    imgar.close()
    return tbl, exp
#########################################################################################################
def createSrc(tbl, fitx, fity, A_bgd, wfile):
    src = []
    j = 0
    tbl.sort(['peak_value'])
    for i in range(1,len(tbl)+1):
        xpeak, ypeak, peak = tbl['x_peak','y_peak','peak_value'][len(tbl)-i]
        cent = [xpeak,ypeak]
        idx = find_nearest(fity*peak,A_bgd*0.03)
        rad = fitx[idx]
        if rad >= 150:
            with open(wfile,'a+') as obsfile:
		obsfile.write("    ***** RADIUS OF EXCL.REG IS LARGER THAN 150 *****"+'\n')
        if peak-np.sqrt(peak) <= 3*A_bgd:
            j += 1
	    with open(wfile,'a+') as obsfile:
            	obsfile.write("You have "+str(j)+"regions that show low sig."+'\n')
        src = np.append([src],[cent[0],cent[1],rad])
        if rad >= 200:
            return src
    return src
#########################################################################################################
def find_nearest(array,bgd_level):
    idx = (np.abs(array-bgd_level)).argmin()
    return idx
#########################################################################################################
"""Need to convert my writing of errors and such to logging"""
if __name__ == "__main__":
    main()

