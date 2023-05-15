from __future__ import absolute_import, print_function, division

import numpy as np
import copy
import astropy.units as u
from astropy.modeling.rotations import rotation_matrix



def rot(ra, dec, rol):
  r11 = np.cos(ra)*np.cos(dec)
  r12 = np.cos(dec)*np.sin(ra)
  r13 = -np.sin(dec)
  r21 = np.sin(rol)*np.sin(dec)*np.cos(ra)-np.cos(rol)*np.sin(ra)
  r22 = np.sin(rol)*np.sin(dec)*np.sin(ra)+np.cos(rol)*np.cos(ra)
  r23 = np.sin(rol)*np.cos(dec)
  r31 = np.cos(rol)*np.sin(dec)*np.cos(rol)+np.sin(rol)*np.sin(ra)
  r32 = np.cos(rol)*np.sin(dec)*np.sin(rol)-np.sin(rol)*np.cos(ra)
  r33 = np.cos(dec)*np.cos(rol)
  return r11, r12, r13, r21, r22, r23, r31, r32, r33

def sun_to_xyz(RA,DEC):
  ra = RA*np.pi/180.
  de = DEC*np.pi/180.
  SA = np.array([[np.cos(de)*np.cos(ra)],[np.cos(de)*np.sin(ra)],[np.sin(de)]])
  return SA

def matrix(RA, DEC, ROLL):
  ra = RA*np.pi/180.
  de = np.abs(90.-DEC)*np.pi/180.
  ro = ROLL*np.pi/180.
  RM = np.array([[np.cos(ro),np.sin(ro),0],[-np.sin(ro),np.cos(ro),0],[0,0,1]])
  PM = np.array([[1,0,0],[0,np.cos(de),np.sin(de)],[0,-np.sin(de),np.cos(de)]])
  YM = np.array([[np.cos(ra),np.sin(ra),0],[-np.sin(ra),np.cos(ra),0],[0,0,1]])
  return RM,PM,YM

def OBF_sun(ra,dec,roll,sra,sdec):
  RollM,DecM,RaM = matrix(ra,dec,roll)
  SA = sun_to_xyz(sra,sdec)
  M32 = np.matmul(RaM,SA)
  M21 = np.matmul(DecM,M32)
  M   = np.matmul(RollM,M21)
  return M

def newsun(rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz,sra,sdec,srol):
  newra = rxx*sra+rxy*sdec+rxz*srol
  newdec = ryx*sra+ryy*sdec+ryz*srol
  newrol = rzx*sra+rzy*sdec+rzz*srol
  return newra, newdec, newrol

def euler_to_quat(yaw, pitch, roll):
  cy = np.cos(yaw*0.5)
  sy = np.sin(yaw*0.5)
  cp = np.cos(pitch*0.5)
  sp = np.sin(pitch*0.5)
  cr = np.cos(roll*0.5)
  sr = np.sin(roll*0.5)

  w = cr * cp * cy + sr * sp * sy
  x = sr * cp * cy - cr * sp * sy
  y = cr * sp * cy + sr * cp * sy
  z = cr * cp * sy - sr * sp * cy

  return w,x,y,z

def quat_to_euler(w,x,y,z):
    sinr_cosp = 2 * (w * x + y * z)
    cosr_cosp = 1 - 2 * (x**2 + y**2)
	# Roll values
    roll = np.arctan2(sinr_cosp, cosr_cosp)

    sinp = 2 * (w * y - z * x)

	# Pitch values
    pitch = np.where(np.abs(sinp) >= 1,
                     np.sign(sinp) * np.pi / 2,
                     np.arcsin(sinp))

    siny_cosp = 2 * (w * z + x * y)
    cosy_cosp = 1 - 2 * (y**2 + z**2)
	# Yaw values
    yaw = np.arctan2(siny_cosp, cosy_cosp)

    return roll, pitch, yaw

def quat_rot(ref_yaw,ref_pitch,ref_roll,p_yaw,p_pitch,p_roll):
    w,x,y,z = euler_to_quat(ref_yaw,ref_pitch,ref_roll)

    xx = x*x
    yy = y*y
    zz = z*z
    xy = x*y
    xz = x*z
    xw = x*w
    yz = y*z
    yw = y*w
    zw = z*w

    m00 = 1 - 2*(yy+zz) 
    m01 = 2*(xy-zw)
    m02 = 2*(xz+yw)
    m10 = 2*(xy+zw)
    m11 = 1 - 2*(xx + zz)
    m12 = 2*(yz-xw)
    m20 = 2*(xz-yw)
    m21 = 2*(yz+xw)
    m22 = 1 - 2*(xx+yy)
    m33 = 1
    m03,m13,m23,m30,m31,m32=0,0,0,0,0,0

    rot_yaw = m00*p_yaw+m01*p_pitch+m02*p_roll
    rot_pitch = m10*p_yaw+m11*p_pitch+m12*p_roll
    rot_roll = m20*p_yaw+m21*p_pitch+m22*p_roll

    sw, sx, sy, sz = euler_to_quat(p_yaw,p_pitch,p_roll)
    M = np.array([[m00,m01,m02,m03],[m10,m11,m12,m13],[m20,m21,m22,m23],[m30,m31,m32,m33]])
    sun = np.array([sw,sx,sy,sz])
    sun_rt = np.matmul(M,sun)  
    sun_rot_1, sun_rot_2, sun_rot_3 = quat_to_euler(sun_rt[0], sun_rt[1],sun_rt[2],sun_rt[3])
    return rot_yaw, rot_pitch, rot_roll, sun_rot_1, sun_rot_2, sun_rot_3

def q_rot(ya, pit, rol, sra, sdec):
    w,x,y,z = euler_to_quat(ra, pit, rol)
    vax = np.sin(sdec)*np.cos(sra)
    vay = np.sin(sdec)*np.sin(sra)
    vaz = np.cos(sdec)
    
    q = np.array([w,x,y,z])
    qp = np.array([w,-x,-y,-z])
    V = np.array([0,vax,vay,vaz]) 

    vb = np.matmul(q,V,qp)

def sky_to_tel(attitude, ra, dec):
  if attitude.shape != (3,3):
    raise ValueError('Attitude has to be a 3x3 matrix')

  unitVectorSky = unit_vector_sky(ra,dec)
  inverse_attitude = np.transpose(attitude)

  unitVectorTel = np.dot(inverse_attitude,unitVectorSky)
  nu2, nu3 = polar_angles(unitVectorTel)
  return nu2, nu3

def unit_vector_sky(ra,dec):
  ra_rad = convert_quantity(ra, u.rad, factor = np.deg2rad(1.))
  dec_rad = convert_quantity(dec, u.rad, factor = np.deg2rad(1.))
  vector = np.array([np.cos(ra_rad)*np.cos(dec_rad), np.sin(ra_rad)*np.cos(dec_rad), np.sin(dec_rad)])
  return vector

def convert_quantity(x_in, to_unit, factor=1.):
  x = copy.deepcopy(x_in)
  if isinstance(x, u.Quantity):
    x_out = x.to(to_unit).value
  else:
    x_out = x * factor
  return x_out

def polar_angles(vector, pos_az = False):
  if len(vector) != 3:
    raise ValueError('Input is not a vector or an array of vectors')
  norm = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
  nu2 = np.arctan2(vector[1], vector[0]) * u.rad
  nu3 = np.arcsin(vector[2]/norm) * u.rad

  if pos_az:
    if np.isscalar(nu2.value) and nu2.value < 0.0:
      nu2 += 360.0 * u.deg
    if not np.isscalar(nu2.value) and np.any(nu2.value < 0.0):
      index = np.where(nu2.value < 0.0)[0]
      nu2[index] += 360.0 * u.deg
  return nu2, nu3

def attitude_matrix(ra, dec, pa):

  nu2 = 0.0
  nu3 = 90.0  #This sets a position on the sky that is straight up from NuSTAR

  mv2 = rotation_matrix(-1*-nu2, axis='z')
  mv3 = rotation_matrix(-1*nu3, axis='y')
  mra = rotation_matrix(-1*ra, axis='z')
  mdec = rotation_matrix(-1*-dec, axis='y')
  mpa = rotation_matrix(-1*pa, axis='x')
  m = np.dot(mv3, mv2)
  m = np.dot(mpa, m)
  m = np.dot(mdec, m)
  m = np.dot(mra, m)
  return m

def posangle(attitude, v2, v3):
  v2r = np.radians(v2)
  v3r = np.radians(v3)

  x = -(attitude[2, 0] * np.cos(v2r) + attitude[2, 1] * np.sin(v2r)) * np.sin(v3r) \
        + attitude[2, 2] * np.cos(v3r)
  y = (attitude[0, 0] * attitude[1, 2] - attitude[1, 0] * attitude[0, 2]) * np.cos(v2r) \
        + (attitude[0, 1] * attitude[1, 2] - attitude[1, 1] * attitude[0, 2]) * np.sin(v2r)
  pa = np.degrees(np.arctan2(y, x))
  return pa
  
