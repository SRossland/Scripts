#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

import sys, shutil, os, string, time, datetime, pidly, numpy as np
import subprocess

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
    #try:

# From here all files will be moved to the scratch drive and processed
    print('A')
    os.chdir(scratch)
#    os.system("mkdir -p NuSTAR/OBS")
    dirs = scratch+'/OBS'
#    os.chdir("NuSTAR")
#    src_dir = dir+'/'+obsid
#    dst_dir = dirs+'/'+obsid
#    if os.path.isdir(dst_dir) == True:
#	shutil.rmtree(dst_dir)
#    shutil.copytree(src_dir, dst_dir)

# This is a test to see what files are actually there    
#    for root, directories, filenames in os.walk(dst_dir):
#	for directory in directories:
#	    print os.path.join(root, directory)
#	for filename in filenames:
#	    print os.path.join(root, filename)

    #os.system("cp -r -p "+dir+"/"+obsid+" ./OBS")
    #dirs = scratch+'/NuSTAR'
    #dirs = dirs+'/OBS'

    os.system("gunzip -d -f "+dirs+"/"+obsid+"/event_uf/*.gz")
    sttime = str(datetime.datetime.now())
    print(obsid +"--- Start: "+ sttime + '\n')
    clder = dirs+'/'+obsid+'/event_defcl'
    if os.path.isfile(clder+'/nu'+obsid+'A02_cl.evt') and os.path.isfile(clder+'/nu'+obsid+'B02_cl.evt') == True:
        ab = ['A','B']
    elif os.path.isfile(clder+'/nu'+obsid+'A02_cl.evt') == True:
        ab = 'A'
    elif os.path.isfile(clder+'/nu'+obsid+'B02_cl.evt') == True:
        ab = 'B'
    else:
        sys.exit()
    idl = pidly.IDL('/uufs/chpc.utah.edu/sys/pkg/idl/8.4/idl84/bin/idl')
    try:
	    print('B')
	    obs_dir = dirs+'/'+obsid
	    sys.stdout.flush()
            for i in range(len(ab)):
                subprocess.call(['/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/Scripts/run_pipe_usrgti_notstrict_02.sh', obs_dir, ab[i]])

    	    print('E')
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

	    endtime = str(datetime.datetime.now())
            print("    FINISHED:  "+endtime+'\n')

# The next 3 lines are to copy only the changed files back over to the origninal source             
#            os.system("cp -r -u ./OBS/"+obsid+" "+dir)
#            os.chdir(scratch)
#            shutil.rmtree(scratch+"/NuSTAR/OBS/"+obsid)
#            os.chdir(dirnow)
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

