#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./start.py

# Created by: Steven P. Rossland, 2017

#"""This will pull all files from the list given and place them in the directory below your
#   parent current working directory (so the one you are in right now). In other words, place
#   your list in the directory above you want to put all your files."""

import sys, traceback, os, string, time, datetime, pidly, logging, numpy as np
from photutils import find_peaks
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.convolution import Gaussian2DKernel, convolve

#mlist=sys.argv[1]

#"""Print a warning to inform the usr of how the environemnt should be setup"""
print("="*70)
print("It is assumed that this program is being ran from the directory where")
print("                all obsID files will be located")
print("       ******If not, you are going to have a bad day******")
print("Files will be retrieved from https://heasarc.gsfc.nasa.gov/FTP/nustar/data/obs")
print("+"*70)
print("   NOTE: YOUR LIST SHOULD BE SAVED SO THE LAST ENTRY FOR THE LINE IS BLANK")
print("+"*70)
print("="*70)
waitforenter=raw_input("Are you in the right DIR? [y,n] -- if not, please fix it!" )
if waitforenter == 'n':
    sys.exit()
print("It's your funeral!")

dir = os.getcwd()
####If the name of this file is changed, it has to be changed in pidly at line ~850 ish for halt and stop


####################################################################################################
"""Main: The main part of this script is an effort to take the whole of the master list given
    to retrieve the files from heasarc and to process them with light curves to fit them witha function
    that can be used to clip the points that are not of importance. From there the file is then reprocessed
    using nupipline with the new gti and then seperated based on header flags (sepsunocc.py) to be
    used at a later date.  Any error is written to the file given above"""

"""Future work:  adjust the src.reg file given (excl.reg) to include all of the source so the background
    is preseverd.
    process the files with getspecarf and getspecnoarf to get spectrum for the fore/background
    nuskybgd_fitab
    create background image
    """
####################################################################################################
####################################################################################################
os.system("clear")
print("Welcome to NuSTAR's: Choose Your Own Adventure!")
print("Please take a moment and select which path your life will take")
print("-"*75)
print("         [1]:  Download a list of observations")
print("         [2]:  Process observations from DLobsid.dat list")
print("-"*75)
choice = eval(raw_input("Which adventure would you like? [1 or 2]:  "))

def main():
    if choice == 1:
	dllog = os.pardir+'/LogFile/DownLog.txt'
	lwit = "a+"
        print("Download list should be located one directory above where you")
        print("want to download to and is your current working directory")
        mlist = raw_input("File list (just the name, not the path):  ")
        wlist, crap = mlist.split(".")
	with open(dllog,lwit) as lwr:
		lwr.write(wlist+"--"+str(datetime.datetime.now())) 
	filepath = os.path.join(os.pardir+'/LogFile',mlist)
        with open(filepath) as file:
            last_line = file.readlines()[-1]
        if last_line == '# fin':
            print("File has already been ran")
            sys.exit()
        f = open(filepath,'r+')
        for line in f.readlines():
            if line[0] == '#':
                continue
            li = line.rstrip()
            if len(li) == 0:
                break
            obsidli = li[-12:-1]
            
            if os.path.isdir(dir+'/'+obsidli) == True:
                print(obsidli+" has been downloaded already")
                continue
            os.system(line)
            with open(os.pardir+"/LogFile/DLobsid"+wlist+".dat","a+") as o:
                o.write(obsidli +'--'+str(datetime.datetime.now())+'\n')
        with open(filepath, 'a+') as end:
            end.write('\n'+'# fin')
	with open(dllog,lwit) as lwr:
		lwr.write("fin --"+str(datetime.datetime.now()))
        f.close()
        sys.exit()
###########################################################################################
    if choice == 2:
        # This file corresponds to the same file given to pidly.py to write to when halt and stop errors
        # are given
        file = os.pardir+"/LogFile/OBSIDLOG.txt"
        fwit = "a+"

        if os.path.isfile(file):
            with open(file,fwit) as obsfile:
                obsfile.write('\n'+"="*90+"\n New list \n")

        with open(file,fwit) as obsfile:
            obsfile.write("#University of Utah Nustar archive \n" \
                      "#File started: "+str(datetime.datetime.now())+" \n#Masterlist: \n \n"+"="*90+'\n')
    #    f = open( dir+'/'+mlist, 'r')

    #    for line in f.readlines():
    #        try:
    #            print("Retrieving: "+line)
    #            index = line.find('-nH')
    #            lineret = line[:index] + '--show-progress '+line[index:]
    #            os.system(lineret)
    #            print("FILE RETRIEVED")

        # """Get the obsid from the input line, this is done by assuming that the last character is a blank space,
        #    i.e. ##/#//#OBSID#####/_ ::: where _ is indicating a blank space"""
    #            getobs=line[-13:]
    #            obsid=getobs[:11]


        #    """Create a timestamp for the log of OBSIDs processed"""
        mlist = raw_input("What download list does this correspond to?:  ")
        f = open(dir+'/DLobsid'+mlist+'.dat','r+')
        for line in f.readlines():
            try:
                obsid = line[0:11]
                if os.path.isfile(dir+'/'+obsid+'/event_defcl') == 'True':
                    with open(file,fwit) as obsfile:
                        obsfile.write("Processed file detetected, assuming no further processing to be taken")
                    continue
            
                sttime = str(datetime.datetime.now())
                with open(file,fwit) as obsfile:
                    obsfile.write(obsid +"---"+ sttime + '\n')
                with open('OTFfile.txt','a+') as OTF:
                    OTF.write('\n'+'#'+obsid +"--"+sttime+"  ")

                excltrue = (['y','n'])
                etrue = excltrue[1]
                cldir = dir+'/'+obsid+'/event_cl'
            
                idl = pidly.IDL('/Users/Steve/Documents/Exelis/idl84/bin/idl')
                
                os.system("gunzip -d "+obsid+"/event_cl/*.gz")
                #os.system("gunzip -d "+obsid+"/event_cl/*.gz")

    #            if os.path.isfile(dir+"/"+obsid+"/event_cl/nu"+obsid+"A01_src.reg"):
    
            #os.system("cp "+dir+"/"+obsid+"/event_cl/nu"+obsid+"A01_src.reg "+dir+"/"+obsid+"/event_cl/excl.reg")
            #      o = open(dir+'/'+obsid+'/event_cl/excl.reg','r+')
            #       lines = o.readlines()
            #       w = []
            #       for i in lines:
            #           j = i.replace(' ','')
            #           w.append(j)
            #       o.seek(0)
            # The following is adjusted to o.seek(0) before the with open
            #       with open(cldir+'/excl.reg',fwit) as o:
            #           o.write('# Region file format: DS9 version 4.1\n')
            #           o.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1                      highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n physical\n')
                    # str1 = ''.join(w)
                    # This is not needed, as far as I can tell, the radius is always 20 pixels at 500,500
                    # for line in str1:
                    #    o.write(str.lower(str1))
                    #o.close()
                

                idl("cd, current=dir")
                idl("obsid='"+obsid+"'")
                idl("imname='im3to30keV.fits'")
                idl("cldir=dir+'/'+obsid+'/event_cl/'")
                idl("mkimgs,cldir,obsid,'A',3,30")
                idl("mkimgs,cldir,obsid,'B',3,30")
                idl("fits_read,cldir+'imA3to30keV.fits',im1,h")
                idl("fits_read,cldir+'imB3to30keV.fits',im2")
                idl("im=im1+im2")
                idl("fits_write,cldir+imname,im,h")
                idl("imsmooth = GAUSS_SMOOTH(im,10)")
                idl("fits_write,cldir+'imsmooth3to30keV.fits',imsmooth,h")
                idl("tbin = 100")
                with open(file,fwit) as obsfile:
                    obsfile.write("    *Images 3 to 30keV    ")

                tbl, imgar, pvalues = peaks(cldir+'/im1to30keV.fits')
                tblcopy = tbl.copy()
    #pcopy = pvalues.copy()

                radius = 20
                srcarray = []
                if np.sum(pvalues) > 100:
                    while np.sum(pvalues) > 50:
                        pcopy = pvalues.copy()
                        p2 = pvalues.copy()
                        src, mask, center  = createSrc(tblcopy, radius, pcopy)
                        if np.sum(src) <= 100:
                            pvalues[mask] = 0
                            radius = 5
                            for i in range(len(tblcopy)):
                                if mask[(tblcopy['x_peak'][i]),(tblcopy['y_peak'][i])] == True:
                                    tblcopy['peak_value'][i]=0
                            pcopy = None
                            p2 = None
                            src = None
                            srctest = None
                            continue
                        radius2 = radius + 5
                        srctest, mask2, center2 = createSrc(tblcopy, radius2, p2)
                        radius += 5
                        if np.sum(src) == np.sum(srctest):
                            etrue = excltrue[0]
                            srcarray = np.append([srcarray],[center[0],center[1],radius])
                            pvalues[mask2] = 0
                            for i in range(len(tblcopy)):
                                if mask2[(tblcopy['x_peak'][i]),(tblcopy['y_peak'][i])] == True:
                                    tblcopy['peak_value'][i]=0
                            radius = 5
                        pcopy = None
                        p2 = None
                else:
                    etrue = excltrue[1]
                    with open(file,fwit) as obsfile:
                        obsfile.write("   ******NO SOURCE NEEDED*******  \n")
                    os.system("rm -f "+cldir+'excl.reg')
                        #Need to erase all values in the tbl associated with mask
                srcarray = np.reshape(srcarray,(len(srcarray)/3,3))

                if etrue == 'y':
                    with open(cldir+'/excl.reg',fwit) as o:
                        o.write('# Region file format: DS9 version 4.1\n')
                        o.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n physical\n')
                    for i in range(len(srcarray)):
                        with open(dir+'/'+obsid+'/event_cl/excl.reg','a+') as o:
                            o.write("circle("+str(srcarray[i,0])+","+str(srcarray[i,1])+","+str(srcarray[i,2])+")"+'\n')
                    with open(file,fwit) as obsfile:
                                obsfile.write("*Region file included \n")
                else:
                    with open(file,fwit) as obsfile:
                        obsfile.write("*No Region file       \n")


                idl("lcfilter, dir, obsid, 'A', 50, 160, imname, tbin")
                idl("lcfilter, dir, obsid, 'A', 50, 160, imname, tbin, /usr")
                if etrue == 'y':
                    idl("lcfilter, dir, obsid, 'A', 3, 20, imname, tbin, /usr, /excl")
                    with open(file,fwit) as obsfile:
                        obsfile.write("    *3 A hist/lc files")
                else:
                    with open(file,fwit) as obsfile:
                        obsfile.write("    *2 A hist/lc files")

                idl("lcfilter, dir, obsid, 'B', 50, 160, imname, tbin")
                idl("lcfilter, dir, obsid, 'B', 50, 160, imname, tbin, /usr")
                if etrue == 'y':
                    idl("lcfilter, dir, obsid, 'B', 3, 20, imname, tbin, /usr, /excl")
                    with open(file,fwit) as obsfile:
                        obsfile.write("    *3 B hist/lc files \n")
                else:
                    with open(file,fwit) as obsfile:
                        obsfile.write("    *2 B hist/lc files \n")

                os.system("mv "+cldir+"/nu"+obsid+"*01_usrgti.fits "+dir+"/"+obsid+"/")
                os.system("mv "+cldir+" "+dir+"/"+obsid+"/event_defcl")
                os.system("/Users/Steve/Documents/Research/Scripts/run_pipe_usrgti_notstrict.sh "+dir+"/"+obsid+" A")
                os.system("/Users/Steve/Documents/Research/Scripts/run_pipe_usrgti_notstrict.sh "+dir+"/"+obsid+" B")
    #idl("reproc, dir, obsid")
                with open(file,fwit) as obsfile:
                    obsfile.write("    *Reprocessed A and B")
                idl("cd, current=dir")
                idl("obsid='"+obsid+"'")
                idl("imname='im3to30keV.fits'")
                idl("cldir=dir+'/'+obsid+'/event_cl/'")
                idl("mkimgs,cldir,obsid,'A',3,30")
                idl("mkimgs,cldir,obsid,'B',3,30")
                idl("fits_read,cldir+'imA3to30keV.fits',im1,h")
                idl("fits_read,cldir+'imB3to30keV.fits',im2")
                idl("im=im1+im2")
                idl("fits_write,cldir+imname,im,h")
                with open(file,fwit) as obsfile:
                    obsfile.write("     *New images made 3 to 30keV \n")
                idl.close()

    #"""System commands to ---rename the given event_cl file, nupipeline, and--- seperate the gti's based on
    #            SUNSHINE and OCCULTED flags"""
                #os.system("mv "+dir+'/'+obsid+"/event_cl "+dir+'/'+obsid+"/event_arch_cl")
                #os.system("run_pipe_notstrict.sh "+obsid)
                os.system("sepsunocc2.py "+dir+" "+obsid+" A")
                os.system("sepsunocc2.py "+dir+" "+obsid+" B")
                with open(file,fwit) as obsfile:
                    obsfile.write("    *event_sep_cl created successfully \n")

                endtime = str(datetime.datetime.now())
                with open(file,fwit) as obsfile:
                    obsfile.write("    Completed successfully:  " +endtime+ '\n')

            except Exception as e:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                with open(file,fwit) as obsfile:
                   obsfile.write('\n'+"="*41+" ERRORS "+"="*41+ '\n')
                   obsfile.write(str(exc_type) + '\n' + str(e) + '\n' + str(exc_tb.tb_lineno))
                   obsfile.write('\n' + "="*90 + '\n')
                continue
            except IDLError as er:
                with open(file,fwit) as obsfile:
                    obsfile.write("="*90+'\n'+"Caught IDL error:  " + str(er)+"\n"+"="*90+'\n')
                idl.close()
                continue

        f.close()


        print("START.PY COMPLETE")
#########################################################################################################
def createCircularMask(center,radi):
    X,Y = np.ogrid[:1000,:1000]
    dist_from_center = np.sqrt((X-center[0])**2 + (Y-center[1])**2)
    mask = dist_from_center <= radi
    return mask

#########################################################################################################
def peaks(file):
    imgar = fits.open(file)
    dataar = imgar[0].data
    gauss_kernel = Gaussian2DKernel(6)
    dataarray = convolve(datar,gauss_kernel)
    datastats = np.where(dataarray > 0)
    mean, median, std = sigma_clipped_stats(dataarray[datastats], sigma=3.0)
    threshold = median + (10*std)
    tbl = find_peaks(dataarray, threshold, box_size=5.0)
    pvalues = np.zeros((1000,1000))
    for i in range(len(tbl)):
        x,y = tbl['x_peak','y_peak'][i]
        pvalues[x,y] = tbl['peak_value'][i]
    return tbl, imgar, pvalues
#########################################################################################################
def createSrc(tbl, rad, pval):
    tbl.sort(['peak_value'])
    xpeak, ypeak = tbl['x_peak','y_peak'][len(tbl)-1]
    cent = [xpeak,ypeak]
    masked = createCircularMask(cent,rad)
    pval[~masked] = 0
    return pval, masked, cent
#########################################################################################################
"""Need to convert my writing of errors and such to logging"""
if __name__ == "__main__":
    main()

