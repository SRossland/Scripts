#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python
# syntax: ./caldb_corr.py obstype obsid

# This runs in the batch.slurm file to reprocess the files given
# obsid is the id of the observation

# obstype is what kind of observation (01 or 02)

import sys, os

obsid = sys.argv[2][0:11]
obstype = sys.argv[1]

# What I would like to do is move the file to the scratch drive, process it
# scratch = /scratch/local/u1019304

work_dir = '/scratch/kingspeak/serial/u1019304/NuSTAR'
if obstype == '01':
	runpipe = 'run_pipe_usrgti_notstrict.sh'
elif obstype == '02':
	runpipe = 'run_pipe_usrgti_notstrict_02.sh'
	
try:
	for ab in ['A','B']:
		os.system("/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/Scripts/"+runpipe+" "+work_dir+"/"+obsid+" "+ab)
		with open('/uufs/astro.utah.edu/common/home/u1019304/temp/batchcomplete.txt','a+') as O:
                	O.write(obsid+' '+ab+'\n')
	if obstype == '02':
		os.system("/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/Scripts/elev_corr.py "+work_dir+" "+obsid+" 0")
		
except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(str(exc_type) + '\n' + str(e) + '\n' + str(exc_tb.tb_lineno))
	with open('/uufs/astro.utah.edu/common/home/u1019304/temp/batcherrors.txt','a+') as O:
		O.write(obsid+' '+ab+' '+str(e)+'\n')
        pass
	
