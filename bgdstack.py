#!/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
# syntax: ./bgdstack.py

# Script will d/l a list of background images from heasarc and use idl to read in the fits
# images and stack them in an effort to find a reasonable value for the standard background value

import sys, os, logging, pidly, string, datetime, time, numpy as np

print("="*80)
print("      Choose what you would like to do:")
print("         [1] : Download a list of files")
print("         [2] : Create gti's, reprocess, and make images of those file")
print("         [3] : Stack those images")
print("="*80)
choice = eval(raw_input("Choice:  "))

dir = os.getcwd()
bgd_log = logging.getLogger(__name__)

obstest = []
if os.path.isfile(dir+'/bgdobsid.dat'):
    c = open(dir+'/bgdobsid.dat','r+')
    for line in c.readlines():
        obs = line[0:11]
        obstest.append(obs)
mlist = list
def main():
    if choice == 1:
        mlist = raw_input("File list:  ")
        with open(dir+'/'+mlist) as file:
            last_line = file.readlines()[-1]
        if last_line == '# fin':
            print("File has already been ran")
            sys.exit()
        f = open(dir+'/'+mlist,'r+')
        for line in f.readlines():
            bgd_log.debug("main")
            if line[0] == '#':
                continue
            li = line.rstrip()
            if len(li) == 0:
                break
            obsidli = li[-12:-1]

            if obsidli in obstest:
                print(obsidli+" has been fully downloaded already")
                continue
            index = line.find('-nH')
            lineret = line[:index] + '--show-progress '+line[index:]
            os.system(lineret)
            with open("bgdobsid.dat","a+") as o:
                o.write(obsidli +'--'+str(datetime.datetime.now())+'\n')
        with open(dir+'/'+mlist, 'a+') as end:
            end.write('\n'+'# fin')
        f.close()
        sys.exit()
####################################################################################
    if choice == 2:
        bgd_log.debug("main")
        obsid=[]
        file = dir+"/OBSIDLOG.txt"
        fwit = "a+"
        with open(file,fwit) as obsfile:
            obsfile.write("# University of Utah NuStar Archive log:  \n"+"# Bgdstack.py run \n"+"="*90+'\n')
        
        print("Low Energy: 3"+"\n"+"HighEnergy: 79"+"\n")
        ch = eval(raw_input("Change values?  [y/n]  (lowercase y or n only): "))
        if ch == 'y':
            elow = eval(raw_input("Low Energy:  "))
            ehigh = eval(raw_input("High Energy:  "))
        else:
            elow = 3
            ehigh = 79
        f = open(dir+'/bgdobsid.dat','r+')
        for line in f.readlines():
            obs = line[0:11]
            obsid.append(obs)
        if os.path.isdir(dir+'/bgdstackimages') == 'False':
            os.system("mkdir "+dir+'/bgdstackimages')
        try:
            for i in range(len(obsid)):
                
                if len(obsid[i]) != 11:
                    print("Not a valid NuStar obsid")
                    continue
            
                cldir = dir+'/'+obsid[i]+'/event_cl/'
                if os.path.isdir(dir+"/"+obsid[i]+"/event_defcl"):
                    print(obsid[i]+" has been processed")
                    continue
                
                with open(file,fwit) as obsfile:
                    obsfile.write(obsid[i] +"---"+ sttime + '\n')

                idl = pidly.IDL('/Users/Steve/Documents/Exelis/idl84/bin/idl')
                
                os.system("gunzip -d "+cldir+"*.gz")
                idl("cd, current=dir")
                idl("obsid='"+obsid[i]+"'")
                idl("imname='im"+elow+"to"+ehigh+"keV.fits'")
                idl("cldir=dir+'/'+obsid+'/event_cl/'")
                idl("mkimgs,cldir,obsid,'A',"+elow+","+ehigh)
                idl("mkimgs,cldir,obsid,'B',"+elow+","+ehigh)
                idl("fits_read,cldir+'imA"+elow+"to"+ehigh+"keV.fits',im1,h")
                idl("fits_read,cldir+'imB"+elow+"to"+ehigh+"keV.fits',im2")
                idl("im=im1+im2")
                idl("fits_write,cldir+imname,im,h")
                with open(file,fwit) as obsfile:
                    obsfile.write("     *Images made "+elow+" to "+ehigh+"keV \n")
                idl("tbin = 100")
                
                idl("lcfilter, dir, obsid, 'A', 50, 160, imname, tbin")
                idl("lcfilter, dir, obsid, 'A', 50, 160, imname, tbin, /usr")
                idl("lcfilter, dir, obsid, 'A',  3,  20, imname, tbin, /usr")
                with open(file,fwit) as obsfile:
                    obsfile.write("    *3 A hist/lc files")
                idl("lcfilter, dir, obsid, 'B', 50, 160, imname, tbin")
                idl("lcfilter, dir, obsid, 'B', 50, 160, imname, tbin, /usr")
                idl("lcfilter, dir, obsid, 'B',  3,  20, imname, tbin, /usr")
                with open(file,fwit) as obsfile:
                    obsfile.write("    *3 B hist/lc files \n")
                os.system("mv "+cldir+"/nu"+obsid[i]+"*01_usrgti.fits "+dir+"/"+obsid[i]+"/")
                os.system("mv "+cldir+" "+dir+"/"+obsid[i]+"/event_defcl")
                os.system("/Users/Steve/Documents/Research/Scripts/run_pipe_usrgti_notstrict.sh "+dir+"/"+obsid[i]+" A")
                os.system("/Users/Steve/Documents/Research/Scripts/run_pipe_usrgti_notstrict.sh "+dir+"/"+obsid[i]+" B")

#########For reprocessed data#######
                with open(file,fwit) as obsfile:
                    obsfile.write("    *Reprocessed A and B")
                idl("cd, current=dir")
                idl("obsid='"+obsid[i]+"'")
                idl("imname='im"+elow+"to"+ehigh+"keV.fits'")
                idl("cldir=dir+'/'+obsid+'/event_cl/'")
                idl("mkimgs,cldir,obsid,'A',"+elow+","+ehigh)
                idl("mkimgs,cldir,obsid,'B',"+elow+","+ehigh)
                idl("fits_read,cldir+'imA"+elow+"to"+ehigh+"keV.fits',im1,h")
                idl("fits_read,cldir+'imB"+elow+"to"+ehigh+"keV.fits',im2")
                idl("im=im1+im2")
                idl("fits_write,cldir+imname,im,h")

                with open(file,fwit) as obsfile:
                    obsfile.write("     *New images made "+elow+" to "+ehigh+"keV \n")

                os.system("mkdetimg.py "+dir+" "+obsid[i]+" "+elow+" "+ehigh)
                os.system("mv "+cldir+'imA'+elow+'to'+ehigh+'keV.fits '+dir+'/bgdstackimages')
                idl("nuskybgd_instamap, "+dir+", "+obsid[i]+",'A','bgd'")
                idl("nuskybgd_instamap, "+dir+", "+obsid[i]+",'B','bgd'")

                os.system("sepsunocc2.py "+dir+" "+obsid+" A")
                os.system("sepsunocc2.py "+dir+" "+obsid+" B")
                with open(file,fwit) as obsfile:
                    obsfile.write("    *event_sep_cl created successfully \n")

                endtime = str(datetime.datetime.now())
                with open(file,fwit) as obsfile:
                    obsfile.write("    Completed successfully:  " +endtime+ '\n')


                idl.close()
        except keyboardInterrupt:
            if eval(raw_input("Did you want to quit?  [y/n]:  ")) == 'y':
                sys.exit()

####################################################################################
    if choice == 3:
        print('Still being built asshole, hold your damn horses')
        sys.exit()

if __name__ == "__main__":
    main()
