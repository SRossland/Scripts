#!/usr/bin/python
# syntax: ./backup.py

"""This code is just to ensure that my files are properly backed up and the newest version are held"""

import string, os, sys, datetime, time, subprocess

def main():
  try:
    yesno = raw_input("Do You need to update your current directory? [y/n]:  ")
    if yesno == 'y':
        path1 = raw_input("What directory would you like to copy (please include last '/')? :  ")
        file1 = raw_input("Give the file name:  ")
    elif yesno == 'n':
        path1 = "/Users/Steve/Documents/Research/Scripts/"
        file1 = "Scripts"
    else:
        print("You need to choose y or n!")
        sys.exit()

    print("This will copy all files in this directory and subdirectory to your backup \n" \
          "with a timedate stamp")
    path2 = "/Volumes/USB_Storage/ResearchBackup/Scripts/"
    dirs  = os.listdir(path1)

    os.system("ln -s "+path2+" "+path1)
    stat = 'stat -f "%Sm" -t "%Y%m%d_%H%M%S" '

    list1 = []
    list2 = []
    i = 0
    date = datetime.datetime.today().strftime('%Y%m%d')

    for line in dirs:
        if line.endswith(".py"):
            list1.append(line)
        if line.endswith(".sh"):
            list2.append(line)

    for line in list1:
        
        if os.path.isfile(path2+line):
            check = subprocess.check_output(stat+path1+line, shell=True)
            bcheck= subprocess.check_output(stat+path2+line, shell=True)
            if check != bcheck:
                
                index = line.find('.py')
                linedate = line[:index]+date+line[index:]
                os.system("mv -- "+path2+line+" "+path2+linedate)
                os.system("cp -R "+path1+line+" "+path2+line)
                i += 1
        else:
            os.system("cp -R "+path1+line+" "+path2+line)
            i += 1

    for line in list2:
        if os.path.isfile(path2+line):
            check = subprocess.check_output(stat+path1+line, shell=True)
            bcheck= subprocess.check_output(stat+path2+line, shell=True)
            if check != bcheck:
                index = line.find('.sh')
                linedate = line[:index]+date+line[index:]
                os.system("mv "+path2+line+" "+path2+linedate)
                os.system("cp -R "+path1+line+" "+path2+line)
                i += 1
        else:
            os.system("cp -R "+path1+line+" "+path2+line)
            i += 1

    print("Backup procedure complete: "+str(i)+" files copied")
    os.system("unlink "+file1)

  except Exception as e:
    print("FAIL! "+ str(e))
    os.system("unlink "+file1)
    sys.exit()

if __name__ == "__main__":
    main()


