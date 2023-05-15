import sys, inspect, os, datetime
import numpy as np
import scipy.spatial as spatial
from os.path import join
from astropy.io import fits



#def main():
##############################################
# Definitions
def Get_RA_DEC_5deg(RA,DEC,dis=5.0):
	coord_dir = dirs('swiftbat_coords')
	f = coord_dir+'/BAT_105m_catalog_16sep2017.txt'
	data = np.loadtxt(f,delimiter='|',skiprows=3,usecols=(2,3,9))
	RA_full = data[:,0]; DEC_full = data[:,1]; flux_full = data[:,2]
	ra_high = RA+dis; ra_low = RA-dis; dec_high = DEC+dis; dec_low = DEC-dis  # dis = 5 if it does't work
	idx_ra = np.where(np.logical_and(RA_full>=ra_low,RA_full<=ra_high))
	idx_dec = np.where(np.logical_and(DEC_full[idx_ra]>=dec_low,DEC_full[idx_ra]<=dec_high))
	RA_within = RA_full[idx_ra[0][idx_dec]]; DEC_within = DEC_full[idx_ra[0][idx_dec]]; flux_within = flux_full[idx_ra[0][idx_dec]]
	de1 = np.radians(DEC_within); de2 = np.radians(DEC); al1 = np.radians(RA_within); al2 = np.radians(RA)
	dist = np.degrees(np.arccos(np.sin(de1)*np.sin(de2)+np.cos(de1)*np.cos(de2)*np.cos(al1-al2)))
	RA_reg, DEC_ret, dist_ret, flux_ret = [],[],[],[]
	for i in range(len(RA_within)):
            if dist[i] <= dis:
		print str(i+1),' RA: ',RA_within[i], ' DEC: ',DEC_within[i], ' Dist: ',dist[i], 'Flux: ',flux_within[i]
		RA_reg.append(RA_within[i]); DEC_ret.append(DEC_within[i]); dist_ret.append(dist[i]); flux_ret.append(flux_within[i])
	RA_re = np.array(RA_reg); DEC_re = np.array(DEC_ret); dist_re = np.array(dist_ret); flux_re = np.array(flux_ret)
	return RA_re, DEC_re, dist_re, flux_ret

def dirs(obj):
	home = '/uufs/astro.utah.edu/common/home/u1019304'
	xray = '/uufs/chpc.utah.edu/common/home/astro/wik'
	if obj == 'swiftbat_coords':
	    return home+'/NuSTAR/Figures/'
	if obj == 'OBS':
	    return xray+'/NuSTAR/OBS'
	if obj == 'logs':
	    return home+'/NuSTAR/Scripts/logs'
	if obj == 'scripts':
	    return home+'/NuSTAR/Scripts/Scripts'
	if obj == 'home':
	    return home
	if obj == 'spec':
	    return home+'/NuSTAR/Scripts/spec'

def fits_header(obsid, det):
	obs_dir = dirs('OBS')+'/'+str(obsid)
	evt_dir = obs_dir + '/event_cl'
	clfile = evt_dir+'/nu'+str(obsid)+str(det)+'01_cl.evt'
	if os.path.isfile(clfile) == True:
	    with fits.open(clfile) as hdul:
		RA = hdul[0].header['RA_PNT']
		DEC = hdul[0].header['DEC_PNT']
		PA = hdul[0].header['PA_PNT']
		exp = hdul[0].header['EXPOSURE']
	    return RA, DEC, PA, exp
	else:
	    print('File '+str(clfile)+' does not exist, sucks to suck')
	    
def Exposure(obsid):
	expA = 0; expB = 0
	eventA = dirs('OBS')+'/'+obsid+'/event_cl/nu'+obsid+'A01_cl.evt'
	eventB = dirs('OBS')+'/'+obsid+'/event_cl/nu'+obsid+'B01_cl.evt'
	if os.path.isfile(eventA) == True:
		with fits.open(eventA) as hdul:
			expA += hdul[0].header['EXPOSURE']
	if os.path.isfile(eventB) == True:
		with fits.open(eventB) as hdul:
			expB += hdul[0].header['EXPOSURE']
	return expA, expB

def obs_dates(obsid):
	eventA = dirs('OBS')+'/'+obsid+'/event_cl/nu'+obsid+'A01_cl.evt'
        eventB = dirs('OBS')+'/'+obsid+'/event_cl/nu'+obsid+'B01_cl.evt'
	dateA, dateB = [],[]
        if os.path.isfile(eventA) == True:
                with fits.open(eventA) as hdul:
                        dateA.append(hdul[0].header['DATE-OBS'])
	else:
		dataA.append('####')
        if os.path.isfile(eventB) == True:
                with fits.open(eventB) as hdul:
                        dateB.append(hdul[0].header['DATE-OBS'])
	else:
		dateB.append('####')
        return dateA, dateB

def GET_OBS(filelist):
	if os.path.isfile(filelist) == True:
		obsid = []
		fl = open(filelist,'r+')
		for line in fl.readlines():
			obsid.append(line.rstrip())
		fl.close()
		return obsid
	elif os.path.isfile(dirs('logs')+'/'+filelist) == True:
		obsid = []
                fl = open(filelist,'r+')
                for line in fl.readlines():
                        obsid.append(line.rstrip())
                fl.close()
                return obsid
	elif os.path.isdir(dirs('OBS')+'/'+filelist) == True:
		return filelist
	else:
		print('No file or obsid found')
		sys.exit()

def list_obs_date(filelist):
 # Need to return the obsid's in the sorted array
	obsids = GET_OBS(filelist)
	dateA, dateB = [],[]
	for i in range(len(obsids)):
		dA, dB = obs_dates(obsids[i])
		if (dA == '####') and (dB == '####'):
			print('no date available',obids[i])
			continue
		if dA == ['####']:
			dA = dB
		if dB == ['####']:
			dB = dA
		dateA.append(dA); dateB.append(dB)
	blah = zip(dateA,obsids)
	blah.sort();dateB.sort()
	dateA_2, obs = zip(*blah)
	return dateA_2, dateB, obs	
		
def list_obs_location(filelist, radial_distance_from_center=5.0):
	obsids = GET_OBS(filelist)
	dis = radial_distance_from_center
	# Get info (RA, DEC)
	points = np.zeros((len(obsids),2))
	RA, DEC, crap, crap2 = [],[],[],[]
	for i in range(len(obsids)):
		if os.path.isfile(dirs('OBS')+'/'+obsids[i]+'/event_cl/nu'+obsids[i]+'A01_cl.evt') == True:
			R,D,P,E = fits_header(obsids[i],'A')
			RA.append(R); DEC.append(D);points[i,0] = R; points[i,1] = D
		else:
			R,D,P,E = fits_header(obsids[i],'B')
			RA.append(R); DEC.append(D);points[i,0] = R; points[i,1] = D
	point_tree = spatial.cKDTree(points)
	# find the centers 
	obs_to_group = obsids[:]
	group_num = 0
	while len(points) > 0.:
		points_near = point_tree.data[point_tree.query_ball_point(points[0],dis)]
		if len(points_near) != 1:
			pt_count = np.zeros(len(points_near))
			for i in range(len(points_near)):
				pt_count[i] = len(point_tree.data[point_tree.query_ball_point(points_near[i],dis)])
		else:
			pt_count = 0
		idx = np.argmax(pt_count)
		points_group = point_tree.data[point_tree.query_ball_point(points_near[idx],dis)]
		with open('OBS_GROUPS_'+dis+'.txt','a+') as O:
			O.write('#'*5+' GROUP: '+str(group_num)+' '+'#'*5+'\n')
		group_num += 1
		for i in range(len(points_group)):
			idx_pt = np.where(points == points_group[i])[0][0]
			obs_pt = obs_to_group[idx_pt]
			with open('OBS_GROUPS_'+dis+'.txt','a+') as O:
				O.write(obs_pt+'\n')
			obs_to_group.pop(idx_pt)
			points = np.reshape(np.delete(points,[idx_pt*2,idx_pt*2+1]),(len(points)-1,2))
		
	return 'DONE'	

def Sun_Roll_Ang(obsid):
	obs_dir = dirs('OBS')+'/'+obsid+'/event_cl/nu'+obsid+'A.attorb'
	obs_dir_B = dirs('OBS')+'/'+obsid+'/event_cl/nu'+obsid+'B.attorb'
	if os.path.isfile(obs_dir) == True:
		obs_dir_true = obs_dir
	elif os.path.isfile(obs_dir_B) == True:
		obs_dir_true = obs_dir_B
	else:
		return 'NO FILE FOUND'
	with fits.open(obs_dir_true) as hdul:
		sunang = hdul[1].data['SUN_ANGLE']
		roll = hdul[1].data['ROLL']
	return sunang, roll
		
def spec_files(det, mode):
	spec_dir = dirs('spec')
	for arg in ['full','nosun','sun']:
	    for arg_2 in ['_det0','_det1','_det2','_det3','']:
		file_root = 'spec'+det+'_'+mode+'_'+arg+arg_2
		file = spec_dir+'/'+file_root+'.pha'
		if os.path.isfile(file) == True:
			os.system('mv '+file+' '+spec_dir+'/spec_archive/'+file_root+datetime.datetime.today().strftime('%Y_%m_%d')+'.pha')
		os.system('cp '+spec_dir+'/spec'+det+'_orig.pha '+file)
		
'''def arguments():
	"""Returns tuple containing dictionary of calling function's
           named arguments and a list of calling function's unnamed
           positional arguments.
        """
	from inspect import getargvalues, stack
	posname, kwname, args = getargvalues(stack()[1][0])[-3:]
	posargs = args.pop(posname, [])
	args.update(args.pop(kwname, []))
	return args, posargs
'''

def getargspec(obj):
    """Get the names and default values of a callable's
       arguments

    A tuple of four things is returned: (args, varargs,
    varkw, defaults).
      - args is a list of the argument names (it may
        contain nested lists).
      - varargs and varkw are the names of the * and
        ** arguments or None.
      - defaults is a tuple of default argument values
        or None if there are no default arguments; if
        this tuple has n elements, they correspond to
        the last n elements listed in args.

    Unlike inspect.getargspec(), can return argument
    specification for functions, methods, callable
    objects, and classes.  Does not support builtin
    functions or methods.
    """
    if not callable(obj):
        raise TypeError, "%s is not callable" % type(obj)
    try:
        if inspect.isfunction(obj):
            return inspect.getargspec(obj)
        elif hasattr(obj, 'im_func'):
            # For methods or classmethods drop the first
            # argument from the returned list because
            # python supplies that automatically for us.
            # Note that this differs from what
            # inspect.getargspec() returns for methods.
            # NB: We use im_func so we work with
            #     instancemethod objects also.
            spec = list(inspect.getargspec(obj.im_func))
            spec[0] = spec[0][1:]
            return spec
        elif inspect.isclass(obj):
            return getargspec(obj.__init__)
        elif isinstance(obj, object) and \
             not isinstance(obj, type(arglist.__get__)):
            # We already know the instance is callable,
            # so it must have a __call__ method defined.
            # Return the arguments it expects.
            return getargspec(obj.__call__)
    except NotImplementedError:
        # If a nested call to our own getargspec()
        # raises NotImplementedError, re-raise the
        # exception with the real object type to make
        # the error message more meaningful (the caller
        # only knows what they passed us; they shouldn't
        # care what aspect(s) of that object we actually
        # examined).
        pass
    raise NotImplementedError, \
          "do not know how to get argument list for %s" % \
          type(obj)	

##############################################

if len(sys.argv) == 1:
      print(" The purpose of this script is to create a general repository of usable scripts that can be called to from anyscript without the need to 'recreate the wheel'.  For a list of availabe scripts pass the argument 'LIST'.")

if len(sys.argv) == 2:
      if sys.argv[1] == 'LIST':
	l = []
	for key, value in locals().items():
          if callable(value) and value.__module__ == __name__:
            l.append(key)
	print l

#if __name__ == "__main__":
#	main()
