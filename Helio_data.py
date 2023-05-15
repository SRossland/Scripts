#!/uufs/astro.utah.edu/common/home/u1019304/VENV2.7.10/bin/python

# syntax: Helio_data.py obslist 


# Purpose: This script is to gather information from a list of NuSTAR 
# Observations about the satellite position (RA, DEC, Roll, etc.) and 
# Find the apperent position of the Sun in ECI and in orbital frames of 
# reference. It will then group those observations into a selection of
# pre-coded sky position segemnts and return a list of OBSIDs that can
# be used in conjunction with another script to make DET1 images.

import numpy as np, sys, os 
from astropy.io import fits
from astropy.coordinates import get_sun 
from astropy.time import Time
from math import sin,cos,atan2,acos
from astropy.coordinates import SkyCoord as sc
############
# Functions
def normalize(v, tolerance=0.00001):
    mag2 = sum(n * n for n in v)
    if abs(mag2 - 1.0) > tolerance:
        mag = sqrt(mag2)
        v = tuple(n / mag for n in v)
    return v

def q_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return w, x, y, z

def q_conjugate(q):
    w, x, y, z = q
    return (w, -x, -y, -z)

def qv_mult(q1, v1):
    q2 = (0.0,) + v1
    return q_mult(q_mult(q1, q2), q_conjugate(q1))[1:]

def axisangle_to_q(v, theta):
    v = normalize(v)
    x, y, z = v
    theta /= 2
    w = cos(theta)
    x = x * sin(theta)
    y = y * sin(theta)
    z = z * sin(theta)
    return w, x, y, z

def q_to_axisangle(q):
    w, v = q[0], q[1:]
    theta = acos(w) * 2.0
    return normalize(v), theta

def ra_dec_to_xyz(ra,dec):
    x = cos(ra)*cos(dec)
    y = sin(ra)*cos(dec)
    z = sin(dec)
    return x,y,z

def get_vector(q,v1):
    v2 = qv_mult(q,v1)
    return np.round(v2,decimals=4)
############


# initial parameters

obslist = sys.argv[1]

list_path = "/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/Scripts/logs"

if not "/uufs/astro.utah.edu/" in obslist:
    obslist = os.path.join(list_path,obslist)

if not os.path.isfile(obslist): 
    print("Can not find file {}, please check path!".format(obslist))
    sys.exit()

obs = np.genfromtxt(obslist,dtype="|U11")

# information I need to get:
# Date for sun position
# attitude quaternion for rotation of sun vector

# The assumptions that are being made are that the Sun has a nominal motion when compared to the motion of the satellite during a typical orbit. 
# 

obs_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS'
write_dir = '/uufs/astro.utah.edu/common/home/u1019304/NuSTAR/SolarData'

for i,ob in enumerate(obs):
    eventA = obs_dir+'/{}/event_cl/nu{}A.attorb'.format(ob,ob)
    eventB = eventA.replace('A.attorb','B.attorb')
    if os.path.isfile(eventA):
        event = eventA
    else:
        event = eventB
    with fits.open(event) as hdul:
        quat = hdul[1].data['QUATERNION']
        mjd = hdul[0].header['MJDREFI']
        mjdf = hdul[0].header['MJDREFF']
    tim = Time(mjd+mjdf, format='mjd')
    sun_coords = get_sun(tim)
    sun_vect = (sun_coords.ra.radian, sun_coords.dec.radian, sun_coords.distance.au)
    Quater = quat[len(quat)/2]
    sun_xyz = ra_dec_to_xyz(sun_vect[0],sun_vect[1])
    v2 = get_vector(Quater,sun_xyz)
    v2norm = v2/np.sqrt(np.sum(v2**2))
    c = sc(x=v2norm[0],y=v2norm[1],z=v2norm[2],representation='cartesian')
    c.representation='spherical'
    with open(write_dir+'/Solar_vectors.txt','a+') as O:
        O.write('{} {} {} {} {} {}\n'.format(ob,v2[0],v2[1],v2[2],c.ra.degree,c.dec.degree))



# The rotation will be qv_mult(q,v1), where q is the quaternion of the attitude
# of the satelite, and v1 is a vector that represents the Sun position
# The sun position will be achieved from get_sun by passing the start date
# of the obs from the observation information
# The Sun position will be given in ECI coords and output as such a vector but
# in the coordinate frame of the satellite




