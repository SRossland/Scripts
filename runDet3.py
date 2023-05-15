import os


home_dir = '/uufs/astro.utah.edu/common/home/u1019304'
stat_program = '{}/NuSTAR/Scripts/Scripts/count_stat_fixeddetvals.py'.format(home_dir)

write_dir = '{}/test_FPMA_noDET3/fits_file/'.format(home_dir)
data_dir = '/uufs/chpc.utah.edu/common/home/astro/wik/NuSTAR/CXB/01/full/3.0_43.76_25steps_112GAP_for_det3_test/fits_file/'

E_low = [3.0, 3.36, 3.76, 4.2, 4.68, 5.24, 5.84, 6.52, 7.28, 8.12, 9.08, 10.16, 11.36, 12.72, 14.24, 15.92, 17.8, 19.92, 22.28, 24.92, 27.88, 31.2, 34.92, 39.08]
E_high = [3.36, 3.76, 4.2, 4.68, 5.24, 5.84, 6.52, 7.28, 8.12, 9.08, 10.16, 11.36, 12.72, 14.24, 15.92, 17.8, 19.92, 22.28, 24.92, 27.88, 31.2, 34.92, 39.08, 43.76]

for i, el in enumerate(E_low):
    os.system('python {} A {}DataA_01_sun_{}_{}keV.fits {} {} full {} 0 0'.format(stat_program, data_dir, el, E_high[i], el, E_high[i], write_dir))

###### This is where I will write the results to a final txt file. ########















