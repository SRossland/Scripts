#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# This program is to simply createa series of observational lists from
# a given list based on parameters of glactic coordinates.

# NOTE: THESE LISTS USE THE UPPER VALUE AS THE SOFT CUTOFF IN LONGITUDE
#       IT IS ASSUMED YOU WILL ALWAYS HAVE AN EXCLUSION REGION AROUND THE
#       GALACTIC CENTER SO THIS IS HANDELED DIFFERENTLY 

import numpy as np, os, sys
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord as sc

input_list = sys.argv[1]  # your input list of obs
div_long = sys.argv[2]    # how many divisions you want in your longitude
div_lat = sys.argv[3]     # how many divisions you want in your latitude
root_file_name = sys.argv[4] # base file name used in the output files
# optional input if you want to change the latitude region being excluded
if len(sys.argv) == 6:
  excl_region_area = sys.argv[5]
else:
  excl_region_area = 30.0 # current point i'm working with

###########################
# definitions
def make_list(long_div, long_low, long_high, lat_div, lat_low, lat_high, excl_range):
# Passing are all parameters needed to create equal areas for the lists
  long_coords, lat_coords = [],[]
# longitude points:
  sets_long = long_div + 1
  step = long_high/sets_long
  for i in np.arange(long_low,long_high,step):
    long_coords.append((i,i+step))

  sets_lat = lat_div
  step = (lat_high-excl_range)/sets_lat
  low_temp = -lat_high
  for i in np.arange(low_temp,-excl_range,step):
    lat_coords.append((i,i+step))
  for i in np.arange(excl_range, lat_high,step):
    lat_coords.append((i,i+step))
  
  num_of_lists = sets_long*(sets_lat+1)
  return long_coords, lat_coords, num_of_lists

###########################


# Assuming you are using the archive established
obs_arch = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS'
output_dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs/'

# all information will be attempted to be retrieved from either the A or B
# cleaned event file (nu{OBS#}{A/B}01_cl.evt)

# create list of obs
obs_list = np.genfromtxt(output_dir+input_list,dtype="|U11")


# make output lists
# how many lists? division of long = 1+div
# division of lat = above*(1+div)

# longitude is 360 degrees, if the divisions do not make a good seperation, 
# an error will be thrown 
coords_long, coords_lat, listnum = make_list(float(div_long), 0., 360.,float(div_lat), -90., 90., float(excl_region_area))

#print(coords_long,coords_lat,listnum)
#sys.exit()
obs_and_pos = []

for i,ob in enumerate(obs_list):
  obs_file = obs_arch+'/'+ob+'/event_cl/nu'+ob+'A01_cl.evt'
  if not os.path.isfile(obs_file):
    obs_file = obs_file.replace('A01_cl.','B01_cl.')
  if not os.path.isfile(obs_file): 
    print('No event file found for '+ob)
    sys.exit()
  with fits.open(obs_file) as hdul:
    RA = hdul[0].header['RA_PNT']
    DEC = hdul[0].header['DEC_PNT']
  c = sc(ra = RA*u.degree, dec = DEC*u.degree, frame = 'fk5')
  obs_and_pos.append((ob, c.galactic.b.deg, c.galactic.l.deg))

length = np.zeros(int(listnum))
end_file_name = map(str,range(int(listnum)))

# failed first attempt, need to do list comprehensions

north = [obs_and_pos[i] for i in range(len(obs_and_pos)) if obs_and_pos[i][1] > 0]
south = [obs_and_pos[i] for i in range(len(obs_and_pos)) if obs_and_pos[i][1] < 0]

n1 = [north[i][0] for i in range(len(north)) if north[i][2] <= 180]
n2 = [north[i][0] for i in range(len(north)) if north[i][2] > 180]
s1 = [south[i][0] for i in range(len(south)) if south[i][2] <= 180]
s2 = [south[i][0] for i in range(len(south)) if south[i][2] > 180]

arrays = [n1,n2,s1,s2]

for i,li in enumerate(arrays):
  for ob in li:
    with open(output_dir+root_file_name+'_'+end_file_name[i], 'a+') as O:
      O.write(ob+'\n') 


