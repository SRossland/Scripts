#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./Masterupkeep.py

# Created by: Steven P. Rossland, 2018

# This script will attempt to update the master list of OBS in the NuSTAR
# catalog kept my University of Utah.  The criteria of selection is by the
# presence of an event_defcl dir, nu############A01_cl.evt (or B, 02) event
# file, or the presence of the obsid on the badobs list.

# Note:  As of Aug. 23, 2018, in each case of criteria, it will only be
# possible to check those obsids passed by the user.

import sys, os, numpy as np
def main():
  print("NuSTAR OBSID upkeep assistant")
  print("=============================")
  print("What can I help you with today?")
  print("  [1] Update processed OBS *")
  print("  [2] Update failed processed OBS")
  print("  [3] Create OBS list *")
  print("  [4] Exposure and count for a list of OBS (coming soon)")
  print("=============================")
  print("* You will be prompted by a second choice")
  print("  ** The screen will clear after this selection!")
  choice = eval(raw_input("  Choose:"))
  while choice not in [1,2,3]:
    print("Not a valid choice, Program exiting!")
    print("SUCK IT!!!!!!")
    sys.exit()
  
  dir_obs = '/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/xray/NuSTAR/OBS'
  dir_list = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs'

  pfile = dir_list+'/Master_list.txt'
  obsid, A_01, B_01, A_02, B_02, fail = getmaster(dir_list)

  if choice == 1:

    os.system("clear")
    print("It is assumed that you are only passing good obs")
    check_A = raw_input("Is this the case? [y/n]")
    if check_A == 'n':
      print("Fix your shit man!")
      sys.exit()
    
    os.system("clear")
    print("Update Processed Data")
    print("=====================")
    print("Which observation type?")
    print("  [1] Science (01 files)")
    print("  [2] Occulted (02 files)")
    print("=====================")
    choice_1 = eval(raw_input("  choice:"))
    while choice_1 not in [1,2]:
      print("Not a valid choice, Program exiting!")
      sys.exit()

    choice_1 = str(choice_1)
    os.system("clear")  
      
    print("What obsids would you like to update?")
    print("=====================================")
    li = raw_input(" File name with file extension (in logs dir):")  

    f = open(li,'r+')
    obs = []
    for line in f.readlines():
      ob = line[0:11]
      obs.append(ob)
      
    print("This will first check for a def_cl dir and 0"+choice_1+" file in event_cl dir")
    print("If both exits, the update will happen.  Otherwise, a printout of ")
    print("failed updates will be listed")

    for i in range(len(obs)):

	idx = obsid.index(obs[i])

        if os.path.isdir(dir_obs+'/'+obs[i]+'/event_defcl') == False:
          print(obs[i]+"   No defcl")
          continue
        if os.path.isfile(dir_obs+'/'+obs[i]+'/event_cl/nu'+obs[i]+'A0'+choice_1+'_cl.evt') == True:
	 if choice_1 == '1': 
	  A_01[idx] = '1'
	 if choice_1 == '2':
	  A_02[idx] = '1'
	else:
	  if choice_1 == '1':
            A_01[idx] = '0'
          if choice_1 == '2':
            A_02[idx] = '0'
	  print(obs[i]+ "  No A0"+choice_1)
  	if os.path.isfile(dir_obs+'/'+obs[i]+'/event_cl/nu'+obs[i]+'B0'+choice_1+'_cl.evt') == True:
	 if choice_1 == '1':
	  B_01[idx] = '1'
	 if choice_1 == '2':
	  B_02[idx] = '1'
	else:
	  if choice_1 == '1':
            B_01[idx] = '0'
          if choice_1 == '2':
            B_02[idx] = '0'
	  print(obs[i]+"   No B01"+choice_1)

    masterprint(pfile, obsid, A_01, B_01, A_02, B_02, fail)

  if choice == 2:
    os.system("clear")
    print("This will set all values except the failed to 0")
    blah = raw_input("Are you sure you want to proceed? [y/n]:")
    while blah not in ['y','Y','N','n']:
      print("Not a valid choice")
      sys.exit()
    if (blah == 'n') or (blah == 'N'):
      sys.exit()
    
    print("What obsids would you like to update?")
    print("=====================================")
    li = raw_input(" File name with file extension (in logs dir):")

    f = open(li,'r+')
    obs = []
    for line in f.readlines():
      ob = line[0:11]
      obs.append(ob)

    for i in range(len(obs)):
      idx = obsid.index(obs[i])
      A_01[idx] = '0'
      B_01[idx] = '0'
      A_02[idx] = '0'
      B_02[idx] = '0'
      fail[idx] = '1'
    
    masterprint(pfile,obsid,A_01,B_01,A_02,B_02,fail)

  if choice == 3:
    os.system("clear")

    print("This will create a list of OBS based on your input")
    print("==================================================")
    print("  [1] Unprocessed list")
    print("  [2] Science Processed")
    print("  [3] Occulted Processed")
    print("  [4] Failed Processed")
    print("  [5] Det A or B not Processed")
    print("==================================================")
    choice_3 = eval(raw_input("  choose:"))
    while choice_3 not in [1,2,3,4,5]:
      print("Not a valid choice")
      sys.exit()

    unfile = raw_input("Give a file name to save the list to:")

    if choice_3 == 1:
      os.system("clear")
      print("This will create a list of unprocessed OBS")
      blah = raw_input('')
      unobs = []

      for i in range(len(obsid)):
        if A_01[i] == '0' and B_01[i] == '0' and A_02[i] == '0' and B_02[i] == '0' and fail[i] == '0':
          unobs.append(obsid[i])

      for i in range(len(unobs)):
        with open(dir_list+'/'+unfile,'a+') as u:
          u.write(unobs[i]+'\n')
      print("There were "+str(i)+" unprocessed OBS")
      sys.exit()

    if choice_3 == 2:
      os.system("clear")

      print("This will give you a list of obs that both det A and B processed")
      print("See choice 5 for A and B choices")
      blah = raw_input('')
      unobs = []

      for i in range(len(obsid)):
        if A_01[i] == '1' and B_01[i] == '1':
	  unobs.append(obsid[i])

      for i in range(len(unobs)):
        with open(dir_list+'/'+unfile,'a+') as u:
          u.write(unobs[i]+'\n')
      sys.exit()

    if choice_3 == 3:
      os.system("clear")
      
      print("This will give you a list of obs that both det A and B processed")
      print("See choice 5 for A and B choices")
      blah = raw_input('')
      unobs = []

      for i in range(len(obsid)):
        if A_02[i] == '1' and B_02[i] == '1':
          unobs.append(obsid[i])
         
      for i in range(len(unobs)):
        with open(dir_list+'/'+unfile,'a+') as u:
          u.write(unobs[i]+'\n')
      sys.exit()         

    if choice_3 == 4:
      os.system("clear")

      print("This will give a list of failed OBS")
      blah = raw_input('')
      unobs = []

      for i in range(len(obsid)):
        if fail[i] == '1':
	  unobs.append(obsid[i])

      for i in range(len(unobs)):
        with open(dir_list+'/'+unfile,'a+') as u:
          u.write(unobs[i]+'\n')
      sys.exit()

    if choice_3 == 5:
      os.system('clear')
 
      print("Would you like a list of 01 or 02 differences?")
      print("==============================================")
      print("  [1] 01")
      print("  [2] 02")
      print("  [3] both")
      print("==============================================")
      print("NOTE: This will give you a list of OBS for A or B that ARE processed")
      ch_3 = eval(raw_input("   choice:"))

      while ch_3 not in [1,2,3]:
        print("Not a valid choice")
        sys.exit()
      
      unobs_A1, unobs_B1, unobs_A2, unobs_B2 = [],[],[],[]
      idx = unfile.index('.')

      for i in range(len(obsid)):
        if (ch_3 == 1) or (ch_3 == 3):
	  if A_01 == '1' and B_01 == '0':
	    unobs_A1.append(obsid[i])
          if A_01 == '0' and B_01 == '1':
    	    unobs_B1.append(obsid[i])

        if (ch_3 == 2) or (ch_3 == 3):
	  if A_02 == '1' and B_02 == '0':
            unobs_A2.append(obsid[i])
          if A_02 == '0' and B_02 == '1':
            unobs_B2.append(obsid[i])
       
      if ch_3 == 1:
	newfileA = unfile[:idx]+'_A01'+unfile[idx:]
        newfileB = unfile[:idx]+'_B01'+unfile[idx:]
	
	for i in range(len(unobs_A1)):
	  with open(dir_list+'/'+newfileA,'a+') as o:
	    o.write(unobs_A1+'\n')
	for i in range(len(unobs_B1)):
          with open(dir_list+'/'+newfileB,'a+') as o:
            o.write(unobs_B1+'\n')
	sys.exit()

      if ch_3 == 2:
	newfileA = unfile[:idx]+'_A02'+unfile[idx:]
        newfileB = unfile[:idx]+'_B02'+unfile[idx:]      
	
        for i in range(len(unobs_A2)):
          with open(dir_list+'/'+newfileA,'a+') as o:
            o.write(unobs_A2+'\n')
        for i in range(len(unobs_B2)):
          with open(dir_list+'/'+newfileB,'a+') as o:
            o.write(unobs_B2+'\n')
	sys.exit()

      if ch_3 == 3:
        newfileA1 = unfile[:idx]+'_A01'+unfile[idx:]
        newfileB1 = unfile[:idx]+'_B01'+unfile[idx:]
        newfileA2 = unfile[:idx]+'_A02'+unfile[idx:]
        newfileB2 = unfile[:idx]+'_B02'+unfile[idx:]


        for i in range(len(unobs_A1)):
          with open(dir_list+'/'+newfileA1,'a+') as o:
            o.write(unobs_A1+'\n')
        for i in range(len(unobs_B1)):
          with open(dir_list+'/'+newfileB1,'a+') as o:
            o.write(unobs_B1+'\n')
        for i in range(len(unobs_A2)):
          with open(dir_list+'/'+newfileA2,'a+') as o:
            o.write(unobs_A2+'\n')
        for i in range(len(unobs_B2)):
          with open(dir_list+'/'+newfileB2,'a+') as o:
            o.write(unobs_B2+'\n')
        sys.exit()
###############################################################################
def getmaster(dirl):
  fi = open(dirl+'/Master_list.txt','r+')
  o,a1,b1,a2,b2,f = [],[],[],[],[],[]
  lines = fi.readlines()[2:]
  fi.close()
  for line in lines:
    ob,As,Bs,Ao,Bo,fa = line.split('|')
    o.append(ob.strip())
    a1.append(As.strip())
    b1.append(Bs.strip())
    a2.append(Ao.strip())
    b2.append(Bo.strip())
    f.append(fa.strip())
  return(o,a1,b1,a2,b2,f)
###############################################################################
def masterprint(wfile,obsi, A01, B01, A02, B02, fai):
      with open(wfile, 'w+') as o:
	o.write('{:^11} | {:^15} | {:^15} | {:^6} '.format('OBSID','01 Processed','02 Processed','Failed'))
        o.write('\n')
        o.write('{:^11} | {:^6} | {:^6} | {:^6} | {:^6} | {:^6} '.format('','A','B','A','B',''))
        o.write('\n')
      
      for i in range(len(obsi)):
	with open(wfile,'a+') as o:
	  o.write('{:^11} | {:^6} | {:^6} | {:^6} | {:^6} | {:^15} '.format(obsi[i],A01[i],B01[i],A02[i],B02[i],fai[i]))
	  o.write('\n')
      print("Complete")    	  
      sys.exit()
##############################################################################
if __name__ == "__main__":
    main()
