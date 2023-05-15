#! /bin/python3

import os, sys, numpy as np


end_tag = sys.argv[1]

####### Definitions
def getValues(fi,det,root):
    path = os.path.join(root,fi)
    low_energy = fi.split('_')[2]
    high_energy = fi.split('_')[3]
    with open(path) as f:
        lines = f.readlines()
    check = 'det{} norm'.format(det)
    for line in lines:
        if check in line:
            li = line.rstrip()
            norm_val_raw = li.split()[3]
            norm_val = float(norm_val_raw)/((float(high_energy)-0.04) - float(low_energy))
        else: 
            continue
    return low_energy,high_energy,norm_val;

def get_exp(det,root):
    paramfi = det+'_params.txt'
    path = os.path.join(root,paramfi)
    with open(path) as f:
        lines = f.readlines()
    exp = lines[0].rstrip().split()[5]
    return exp;

##########




root_file = "/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/CXB/01/fixed_det_cstatFit/"
#root_file = "/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/CXB/02/nosun/3.0_41.28_36steps"
#root_file = root_file+"nosun_AB_3_40ish"
root_file = root_file+"full_AB_3_40ish"


write_file = '/uufs/astro.utah.edu/common/home/u1019304/detvals/'
list_of_files = os.listdir(root_file)

A_list,B_list = [],[]

for li in list_of_files:
    if 'A_nosun' in li:
        A_list.append(li)
    if 'B_nosun' in li:
        B_list.append(li)

A_idx, B_idx = [],[]
for id in A_list:
    A_idx.append(float(id.split('_')[2]))
for id in B_list:
    B_idx.append(float(id.split('_')[2]))

A_sort = [x for _,x in sorted(zip(A_idx,A_list))]
B_sort = [x for _,x in sorted(zip(B_idx,B_list))]


sorted_list = [A_sort,B_sort]
dets = ['A','B']

for i in range(4):
    for j in range(2):
        exp = get_exp(dets[j],root_file)
        writefile = "{}{}det{}normvals_{}.txt".format(write_file,dets[j],i,end_tag)
        for k,fi in enumerate(sorted_list[j]):
            low, high, norm = getValues(fi,i,root_file)
            fixednorm = float(norm)/float(exp)
            with open(writefile,'a+') as O:
                O.write("{} {} {}{}".format(low,high,fixednorm,'\n'))







