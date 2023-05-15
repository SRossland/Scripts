#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./fixblank.py

# Created by: Steven P. Rossland, 2018

# This script is, hopefully, a onetime use program to fix the lack of flare
# filtering for obs that showed no sources in low energies 

import sys, os, string, time, datetime, pidly, numpy as np

# the peakfinding is not utilized because it is assumed that no source is 
# detected. 

def main():
	print("It is assumed you are running this from XRAY/NuSTAR/ directory")
	dirnow = os.getcwd()
	dir = dirnow+'/OBS'
	wlist = raw_input("Input list of OBSID's:  (file ext. needed, i.e. *.dat, also assumed to be located in home directory of u1019304/NuSTAR/Scripts/Scripts/)  ")
	mlist, crap = wlist.split(".")
	
	f = open('/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/Scripts/'+wlist,'r+')
	
	for line in f.readlines():
		try: 
			obsid = line[0:11]
		
			if len(line) == 0:
				break
			
			if os.path.isdir(dir+'/'+obsid+'/event_cl/event_defcl') == True:
				os.system('mv '+dir+'/'+obsid+'/event_cl/event_defcl '+dir+'/'+obsid)
			
			if os.path.isdir(dir+'/'+obsid+'/event_defcl') == True:
				if os.path.isdir(dir+'/'+obsid+'/event_cl') == True:
					os.system("rm -r "+dir+'/'+obsid+'/event_cl')
				os.system('mv '+dir+'/'+obsid+'/event_defcl '+dir+'/'+obsid+'/event_cl')
				
			wfile = dir+'/'+obsid+'/Proc_log.txt'
			sttime = str(datetime.datetime.now())
			with open(wfile,'a+') as obsfile:
				obsfile.write("    reprocess:  Start "+sttime+'\n')

			cldir = dir+'/'+obsid+'/event_cl'

			os.system("gunzip -d "+cldir+"/*.gz")

			idl = pidly.IDL('/uufs/chpc.utah.edu/sys/pkg/idl/8.4/idl84/bin/idl')
# This section is the image step			
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

			idl("lcfilter, dir, obsid, 'A', 50, 160, imname, tbin")
            		idl("lcfilter, dir, obsid, 'A', 50, 160, imname, tbin, /usr")
               		idl("lcfilter, dir, obsid, 'A', 3, 20, imname, tbin, /usr, /blank")

            		idl("lcfilter, dir, obsid, 'B', 50, 160, imname, tbin")
            		idl("lcfilter, dir, obsid, 'B', 50, 160, imname, tbin, /usr")
                	idl("lcfilter, dir, obsid, 'B', 3, 20, imname, tbin, /usr, /blank")

            		with open(wfile,'a+') as obsfile:
                		obsfile.write("    GTI file created"+'\n')

            		os.system("mv "+cldir+"/nu"+obsid+"*01_usrgti.fits "+dir+"/"+obsid+"/")
            		os.system("mv "+cldir+" "+dir+"/"+obsid+"/event_defcl")
            		os.system("/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/Scripts/run_pipe_usrgti_notstrict.sh "+dir+"/"+obsid+" A") #
            		os.system("/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/Scripts/run_pipe_usrgti_notstrict.sh "+dir+"/"+obsid+" B") #

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

		except Exception as e:
            		exc_type, exc_obj, exc_tb = sys.exc_info()
            		fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                   	print('\n'+"="*41+" ERRORS "+"="*41+ '\n')
                   	print(str(exc_type) + '\n' + str(e) + '\n' + str(exc_tb.tb_lineno))
                   	print('\n' + "="*90 + '\n')
            	continue
    	print("EXCL.PY COMPLETE")
    	f.close()
     
if __name__ == "__main__":
	main()
