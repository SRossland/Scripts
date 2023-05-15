import numpy as np
import copy, sys, os
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.modeling.rotations import rotation_matrix
from astropy.coordinates import get_sun
from astropy.time import Time
from astropy.io import fits
from astropy.coordinates import SkyCoord as sc
from scipy.optimize import curve_fit
obslist = sys.argv[1]



def attitude_matrix(nu2, nu3, ra, dec, pa):
    if isinstance(nu2, u.Quantity):
        nu2_value = nu2.to(u.deg).value
    else:
        nu2_value = nu2

    if isinstance(nu3, u.Quantity):
        nu3_value = nu3.to(u.deg).value
    else:
        nu3_value = nu3

    if isinstance(pa, u.Quantity):
        pa_value = pa.to(u.deg).value
    else:
        pa_value = pa
    if isinstance(ra, u.Quantity):
        ra_value = ra.to(u.deg).value
    else:
        ra_value = ra
    if isinstance(dec, u.Quantity):
        dec_value = dec.to(u.deg).value
    else:
        dec_value = dec

    # Get separate rotation matrices
    # astropy's rotation matrix takes inverse sign compared to rotations.rotate
    mv2 = rotation_matrix(-1*-nu2_value, axis='z')
    mv3 = rotation_matrix(-1*nu3_value, axis='y')
    mra = rotation_matrix(-1*ra_value, axis='z')
    mdec = rotation_matrix(-1*-dec_value, axis='y')
    mpa = rotation_matrix(-1*-1*pa_value, axis='x')

    m = np.dot(mv3, mv2)
    m = np.dot(mpa, m)
    m = np.dot(mdec, m)
    m = np.dot(mra, m)

    return m


def convert_quantity(x_in, to_unit, factor=1.):
    x = copy.deepcopy(x_in)
    if isinstance(x, u.Quantity):
        x_out = x.to(to_unit).value
    else:
        x_out = x*factor
    return x_out

def polar_angles(vector, positive_azimuth=False):
    if len(vector) != 3:
        raise ValueError('Input is not a vector')

    norm = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
    nu2 = np.arctan2(vector[1], vector[0])*u.rad
    nu3 = np.arcsin(vector[2]/norm)*u.rad

    if positive_azimuth:
        if np.isscalar(nu2.value) and nu2.value < 0.0:
            nu2 += 360.0 *u.deg
        if not np.isscalar(nu2.value) and np.any(nu2.value < 0.0):
            index = np.where(nu2.value < 0.0)[0]
            nu2[index] += 360.0 *u.deg

    return nu2, nu3

def unit_vector_sky(ra,dec):
    ra_rad = convert_quantity(ra, u.rad, factor=np.deg2rad(1.))
    dec_rad = convert_quantity(dec, u.rad, factor= np.deg2rad(1.))
    vector = np.array([np.cos(ra_rad)*np.cos(dec_rad), np.sin(ra_rad)*np.cos(dec_rad), np.sin(dec_rad)])
    return vector



def sky_to_tel(attitude, ra, dec):
    if attitude.shape != (3,3):
        raise ValueError('Attitude has to be a 3x3 matrix')

    unit_vector_sky_side = unit_vector_sky(ra, dec)
    inverse_attitude = np.transpose(attitude)

    unit_vector_tel = np.dot(inverse_attitude, unit_vector_sky_side)

    nu2, nu3 = polar_angles(unit_vector_tel)

    return nu2, nu3


def check_pos(attitude, nu2, nu3, positive_ra=True):
    nu2_deg = convert_quantity(nu2, u.deg, factor=u.arcsec.to(u.deg))
    nu3_deg = convert_quantity(nu3, u.deg, factor=u.arcsec.to(u.deg))
    unit_vector_tel = unit_vector_sky(nu2_deg, nu3_deg)
    unit_vector_sky_side = np.dot(attitude, unit_vector_tel)
    ra, dec = polar_angles(unit_vector_sky_side, positive_azimuth=positive_ra)
    return ra, dec


def parabola(x,a,b,c):
    return a+b*x+c*x**2

# need to feed  
obs = np.genfromtxt(obslist,'|U11')
obsdir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/OBS/'
ra, dec = [],[]
for i,ob in enumerate(obs):
    event = os.path.join(obsdir,'{}/event_cl/nu{}A02_cl.evt'.format(ob,ob))
    if not os.path.isfile:
        event = event.replace('A02_cl','B02_cl')
    with fits.open(event) as hdul:
        NuSTAR_RA = hdul[0].header['RA_PNT']
        NuSTAR_DEC = hdul[0].header['DEC_PNT']
        NuSTAR_PA = hdul[0].header['PA_PNT']
        NuSTAR_OBS_DATE = hdul[0].header['DATE-OBS']
    t = Time(NuSTAR_OBS_DATE, format='isot', scale='tt')
    Sun_coords = get_sun(t)
    Sun_RA = Sun_coords.ra.value
    Sun_DEC = Sun_coords.dec.value
    event_attorb = event.replace('02_cl.evt','.attorb')
    with fits.open(event_attorb) as hdul:
        Sun_Angle = hdul[1].data['SUN_ANGLE']
    SA = np.average(Sun_Angle)

    Rotation_Matrix = attitude_matrix(0., 0., NuSTAR_RA, NuSTAR_DEC, NuSTAR_PA)

    Relative_RA, Relative_Dec = sky_to_tel(Rotation_Matrix, Sun_RA, Sun_DEC)

    ra.append(Relative_RA.value)
    dec.append(Relative_Dec.value)

    Coord_Sun = sc(Relative_RA.value,Relative_Dec, unit='rad', frame='gcrs')
    Coord_NuSTAR = sc(0.0, 0.0, unit='deg', frame='gcrs')
    radial_dist_2points = Coord_NuSTAR.separation(Coord_Sun)
    #print(SA,radial_dist_2points.deg)

ra_deg = [i*180.0/np.pi for i in ra]
dec_deg = [i*180.0/np.pi for i in dec]

#start,pcov = curve_fit( parabola, ra_deg, dec_deg)

#Ylist = [parabola(x, *start) for x in ra_deg] 
ra_deg = [i*-1 for i in ra_deg]
#plt.scatter(dec_deg,ra_deg,s=4)
#plt.scatter(ra_deg, Ylist, color='red', s=4)
#plt.savefig('CartesianSolar2_reverse.png')
#plt.close()
#ra_deg2 = [i*-1 for i in ra_deg]
fig = plt.figure()
ax = fig.add_subplot(projection='polar')
#plt.axes(projection='polar')
#ax.rgrids((0,25,50,75,100,125,150,175),angle=330)
ax.scatter(dec,ra_deg, s=4)
ax.set_rorigin(-1)
ax.set_thetamin(-10)
ax.set_thetamax(75)
#plt.scatter(ra, Ylist, color='red', s=4)
plt.savefig('Polar3.png')





    #with open('Solar_positions_wrt_NuSTAR.txt','a+') as O:
    #    O.write('{} {} {}\n'.format(ob,Relative_RA.value*180/np.pi, Relative_Dec.value*180/np.pi))
    #Check_rotation = attitude_matrix(0., 0., 0., 0., -23.5)
    #Check_RA_Sol, Check_DEC_Sol = sky_to_tel(Check_rotation, 90.0,23.5)
    #Check_RA, Check_DEC = check_pos(Rotation_Matrix, Relative_RA, Relative_Dec)
    
#print(Check_RA_Sol.value*180/np.pi, Check_DEC_Sol.value*180/np.pi)
    #print(Relative_RA*180/np.pi, Relative_Dec*180/np.pi)
    #print(Check_RA*180/np.pi, Sun_RA, Check_DEC*180/np.pi, Sun_DEC)




