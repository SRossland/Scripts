#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./bgdstack_02.py filename datatype

# filename = anytype that is readable 

# datatype = sci or occ
# images are stacked them in an effort to find a reasonable value for the standard background value
#
# Update:  rewrote the program to stack either SCIENCE or OCCULTED data with or without sources.
# It requires that you run your source/no-source lists seprately (if implemented, I need to find a
# way to track the movement of the pixel on the detector, i.e. how the movement happens in det coords. 
# compared to that of gal. coords.).  This is done by passing the list name and if it is sciene or
# occulted data in the program call.



import sys, os, pidly, string, datetime, time, numpy as np
from astropy.io import fits
from numpy import *

########following is added in update#######
listfile = sys.argv[1]

usr_arg = sys.argv[2].lower()

if usr_arg not in ['sci','occ']:
        print('Need a valid datatype [sci or occ], exiting')
        sys.exit()

usr = '01' if usr_arg == 'sci' else '02'
###########################################

print("="*80)
print("      This program uses a file in the logs directory named bgdobsid.dat")
print("         If you do not have this file the program will not run")
print("         Ensure you have all the obsid's you want in this list")
print("         The det images will be saved to the logs directory...just fyi")
print("      [1] To create det images ")
print("      [2] To stack those images ")
print("="*80)
choi = raw_input("Are you ready with your list?  [y/n]")

while choi not in ['y']:
    print('Exiting')
    sys.exit()
    
choice = eval(raw_input("What would you like to do? [1 or 2]"))
while choice not in [1,2]:
    print('Exiting')
    sys.exit()

dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs'
obsdir = "/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS"

def main():
    obsid=[]
    file = dir+"/OBSIDDETLOG_"+usr+".txt"
    fwit = "a+"
    with open(file,fwit) as obsfile:
        obsfile.write("# University of Utah NuStar Archive log:  \n"+"# Bgdstack.py run \n"+"="*90+'\n')
    
    
    print("Low Energy: 3"+"\n"+"HighEnergy: 7"+"\n")
    ch = raw_input("Change values?  [y/n]  (lowercase y or n only): ")
    if ch == 'y':
        elow = str(eval(raw_input("Low Energy:  ")))
        ehigh = str(eval(raw_input("High Energy:  ")))
    else:
        elow = '3'
        ehigh = '7'
        
    f = open(dir+'/'+listfile,'r+')
    for line in f.readlines():
        obs = line[0:11]
        obsid.append(obs)
    f.close()
    if choice == 1:    
        if os.path.isdir(dir+'/bgdstackimages'+usr) == False:
            os.system('mkdir '+dir+'/bgdstackimages'+usr)
        sunocc = str(raw_input("Would you like to do sun/nosun seperation?  [y/n]"))
        
        try:
            for i in range(len(obsid)):
                
                if len(obsid[i]) != 11:
                    print("Not a valid NuStar obsid")
                    continue
            
                idl = pidly.IDL('/uufs/chpc.utah.edu/sys/pkg/idl/8.4/idl84/bin/idl')
            	
		occtype = 'NOOCC' if usr_arg == 'sci' else 'OCC'
		
                cldir = obsdir+'/'+obsid[i]+'/event_cl/'
                if sunocc == 'n':
                    cldir = obsdir+'/'+obsid[i]+'/event_cl/'
                    os.system("mkdetimgs_01.py . "+obsid[i]+" "+elow+" "+ehigh+" "+usr)
                elif sunocc == 'y':
                    cldir = obsdir+'/'+obsid[i]+'/event_sep_cl/'
                    for su in ['SUN','NOSUN']:
			dirs = os.getcwd()
                        os.system("mkdetimgs_01.py "+dirs+" "+obsid[i]+" "+elow+" "+ehigh+" "+usr+" "+occtype+" "+su)
                else:
                    print("Not a valid choice for sun or nosun seperation")
                    sys.exit()
#		for det in ['B']:
                for det in ['A','B']:
                    if sunocc == 'n':
                        if os.path.isfile(cldir+'im'+det+elow+'to'+ehigh+'keVDET.fits') == True:
                            os.system("mv "+cldir+'im'+det+elow+'to'+ehigh+'keVDET.fits '+dir+'/bgdstackimages'+usr+'/im'+det+elow+'to'+ehigh+'keVDET'+obsid[i]+'.fits')
                    if sunocc == 'y':
                        for su in ['SUN','NOSUN']:
                            if os.path.isfile(cldir+occtype+'/'+'im'+det+elow+'to'+ehigh+'keV'+su+'DET.fits') == True:
                                os.system("mv "+cldir+occtype+'/'+'im'+det+elow+'to'+ehigh+'keV'+su+'DET.fits '+dir+'/bgdstackimages'+usr+'/im'+det+elow+'to'+ehigh+'keV'+su+'DET'+obsid[i]+'.fits')
                #os.system("mv "+cldir+'imB'+elow+'to'+ehigh+'keVDET.fits '+dir+'/bgdstackimages02/imB'+elow+'to'+ehigh+'keVDET'+obsid[i]+'.fits')
                idl("dir = '"+obsdir+"'")
                idl("obsid = '"+obsid[i]+"'")
                idl("nuskybgd_instrmap, dir, obsid,'A','bgd'")
                idl("nuskybgd_instrmap, dir, obsid,'B','bgd'")
                idl.close()
        
            sys.exit()
        
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(str(exc_type) + '\n' + str(e) + '\n' + str(exc_tb.tb_lineno) + '\n' + fname)

    if choice == 2:
        try:
            stack_image_A = np.zeros([360,360]) # careful of type of array specifty type
            modA = np.zeros([360,360])
            stack_image_B = np.zeros([360,360])
            modB = np.zeros([360,360])

            count_A = 0
            count_B = 0
            expA = 0
            expB = 0
            
            so = str(raw_input("Is this a SUN/NOSUN seperations?  [y/n]"))
            if so == 'y':
                suno = str(raw_input("SUN or NOSUN? (use all caps!)"))
            for i in range(len(obsid)):
                for det in ["A","B"]:
                    cldir = obsdir+'/'+obsid[i]+'/event_cl/'
                    
                    if so == 'y':
                        if os.path.isfile(dir+'/bgdstackimages'+usr+'/im'+det+elow+'to'+ehigh+'keV'+suno+'DET'+obsid[i]+'.fits') == True:
                            fits_image_filename = dir+'/bgdstackimages'+usr+'/im'+det+elow+'to'+ehigh+'keV'+suno+'DET'+obsid[i]+'.fits'
                        else:
			    print(obsid[i]+' does not exist')
                            continue
                    elif so == 'n':
                        fits_image_filename = dir+'/bgdstackimages'+usr+'/im'+det+elow+'to'+ehigh+'keVDET'+obsid[i]+'.fits'
                    if os.path.isfile(fits_image_filename) == False:
                      	  print(obsid[i]+' does not exist')
			  continue
                    hdul = fits.open(fits_image_filename)
                    hdr = hdul[0].header
                    exp = hdr['Exposure']  # The exposure time per observation
        
                    data =hdul[0].data  ###### #####  Counts per pixel in matrix form (360 x 360 array)
                    #hdul.close()
                    
                    instrmap = cldir+'/bgd/newinstrmap'+det+'.fits'
                    fi = fits.open(instrmap)
                    map_data_full = fi[1].data # Holding the instrument data to a variable
                    #fi.close()
                    
                    # To get the average flux: sum(data)/sum(map_data*exp)
                    # This will all be summed over the total observations, however I do need to ensure that any
                    # infinities and NaN's are taken care of
                    
                    map_data = np.float64(map_data_full > 0.)
                    # map_data.astype(np.int) # integrated into the follwing lines
                
                    modC = exp*map_data  #Fix the mask
                    #datanorm = sum(data)/mod
                    
                    if det == "A":
                        stack_image_A = stack_image_A + data*map_data#.astype(np.int)
                        modA = modA + modC
                        expA = expA + exp
                        count_A += 1
                    elif det == "B":
                        stack_image_B = stack_image_B + data*map_data#.astype(np.int)
                        modB = modB + modC
                        expB = expB + exp
                        count_B += 1
		    else:
			print('Error! Detector not known')

                    hdul.close()
		    fi.close()
                    # both count_A and count_B should be the same value, however, this should also catch any
                    # errors the program has made
                    
	    avg_bgd_A_pix = np.divide(stack_image_A,modA,where=modA!=0)
	    avg_bgd_B_pix = np.divide(stack_image_B,modB,where=modB!=0)

    ##        with np.errstate(divide = 'ignore'):
    ##            avg_bgd_A_pix = stack_image_A/modA
    ##            avg_bgd_B_pix = stack_image_B/modB
            
            # find the nan's and inf's
            where_are_nans_A = isnan(avg_bgd_A_pix)
            where_are_nans_B = isnan(avg_bgd_B_pix)
            where_are_infs_A = isinf(avg_bgd_A_pix)
            where_are_infs_B = isinf(avg_bgd_B_pix)
            
            avg_bgd_A_pix[where_are_nans_A] = 0
            avg_bgd_A_pix[where_are_infs_A] = 0
            avg_bgd_B_pix[where_are_nans_B] = 0
            avg_bgd_B_pix[where_are_infs_B] = 0

            hduAtest = fits.PrimaryHDU(stack_image_A)
            hduBtest = fits.PrimaryHDU(stack_image_B)
            # Write a new fits file:
            hduA = fits.PrimaryHDU(avg_bgd_A_pix)
            hduB = fits.PrimaryHDU(avg_bgd_B_pix)
            
            # These should be a fits file that is a displayable image of the stacked data
            # normalized to the number of observations used in the stack
            if so == 'y':
                hduA.writeto(dir+'/bgdstackimages'+usr+'/avg_bgd_A'+elow+'_'+ehigh+suno+'.fits')
                hduB.writeto(dir+'/bgdstackimages'+usr+'/avg_bgd_B'+elow+'_'+ehigh+suno+'.fits')

                hduAtest.writeto(dir+'/bgdstackimages'+usr+'/avg_bgd_A_test'+elow+'_'+ehigh+suno+'.fits')
                hduBtest.writeto(dir+'/bgdstackimages'+usr+'/avg_bgd_B_test'+elow+'_'+ehigh+suno+'.fits')
            else:
                hduA.writeto(dir+'/bgdstackimages'+usr+'/avg_bgd_A'+elow+'_'+ehigh+'.fits')
                hduB.writeto(dir+'/bgdstackimages'+usr+'/avg_bgd_B'+elow+'_'+ehigh+'.fits')

                hduAtest.writeto(dir+'/bgdstackimages'+usr+'/avg_bgd_A_test'+elow+'_'+ehigh+'.fits')
                hduBtest.writeto(dir+'/bgdstackimages'+usr+'/avg_bgd_B_test'+elow+'_'+ehigh+'.fits')

            # Find the value of the whole detector average background level
            
            avg_bgd_A_sum = 0
            avg_bgd_B_sum = 0
            A_c = 0
            B_c = 0
            
            for i in range(360):
                for j in range(360):
                    if avg_bgd_A_pix[i,j] > 0.:
                        avg_bgd_A_sum += avg_bgd_A_pix[i,j]
                        A_c += 1
                    if avg_bgd_B_pix[i,j] > 0.:
                        avg_bgd_B_sum += avg_bgd_B_pix[i,j]
                        B_c += 1
                        
            #avg_bgd_A = avg_bgd_A_sum/A_c
            #avg_bgd_B = avg_bgd_B_sum/B_c
            
            #print(avg_bgd_A)
            #print(avg_bgd_B)
            print(expA)
            print(expB)
    
            if so == 'y':
                with open(dir+'/bgdstackimages'+usr+'/avg_level_'+suno+elow+'_'+ehigh+'keV.dat','a+') as lvl:
            #lvl.write("Average background level from fits files created from bgdobsid.dat"+'\n')
                    lvl.write(str(avg_bgd_A)+'\n')
                    lvl.write(str(avg_bgd_B))
            
                sys.exit()
            else:
                with open(dir+'/bgdstackimages'+usr+'/avg_level_'+elow+'_'+ehigh+'keV.dat','a+') as lvl:
            #lvl.write("Average background level from fits files created from bgdobsid.dat"+'\n')
                    lvl.write(str(avg_bgd_A)+'\n')
                    lvl.write(str(avg_bgd_B))
            
                sys.exit()
        
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(str(exc_type) + '\n' + str(e) + '\n' + str(exc_tb.tb_lineno) + '\n' + fname)


if __name__ == "__main__":
    main()
