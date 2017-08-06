#!/bin/bash 

rm compute_covariances_fourier

gcc -I/usr/local/include -L/usr/local/lib -o ./compute_covariances_fourier compute_covariances_fourier.c -lfftw3 -lgsl -lgslcblas -lm

./compute_covariances_fourier 100000

# setenv LD_LIBRARY_PATH /home/teifler/lib
# gcc -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier compute_covariances_fourier.c -lfftw3 -lgsl -lgslcblas -lm 

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier compute_covariances_fourier.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass

#  ./covariances_powerspec ../../Tully_Fisher/paralist_Coyote_fid ../../Tully_Fisher/cov_parallel2/ NCoyote_fid_LSST_conti 1 0 1 20 30 5000 1 2 2 2 2 1 1 0 1 4 4 4 4 1 8
#./covariances_powerspec ../../Tully_Fisher/paralist_Coyote_fid ../../emulator/cov_playground2/ s-lhs.100_6_CosmicEmuOnly_range 1 0 1 20 30 5000 1 2 2 2 2 1 1 0 1 4 4 4 4 1 8

 #  ./covariances_powerspec ../../Tully_Fisher/paralist_Coyote_fid ../../Tully_Fisher/cov_parallel2/ NCoyote_fid_LSST_conti 1 0 1 20 30 5000 1 2 2 2 2 1 1 0 1 4 4 4 4 1 8

# 
# ./covariances_powerspec ../../Tully_Fisher/paralist_Coyote_fid ../../Tully_Fisher/cov/ Coyote_fid_DES_conti 1 0 1 20 30 5000 0 2 1 1 1 1 1 0 0 0 0 0 0 1
# 
# ./covariances_powerspec ../../Tully_Fisher/paralist_Coyote_fid ../../Tully_Fisher/cov/ Coyote_fid_LSST_conti 1 0 1 20 30 5000 0 2 3 3 3 1 1 0 0 0 0 0 0 1
	

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier compute_covariances_fourier.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier2 compute_covariances_fourier2.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier3 compute_covariances_fourier3.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier4 compute_covariances_fourier4.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier5 compute_covariances_fourier5.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier6 compute_covariances_fourier6.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier7 compute_covariances_fourier7.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier8 compute_covariances_fourier8.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier9 compute_covariances_fourier9.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier10 compute_covariances_fourier10.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier11 compute_covariances_fourier11.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier12 compute_covariances_fourier12.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier13 compute_covariances_fourier13.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier14 compute_covariances_fourier14.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass


gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier15 compute_covariances_fourier15.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier16 compute_covariances_fourier16.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier17 compute_covariances_fourier17.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier18 compute_covariances_fourier18.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier19 compute_covariances_fourier19.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier20 compute_covariances_fourier20.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier21 compute_covariances_fourier21.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier22 compute_covariances_fourier22.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier23 compute_covariances_fourier23.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../../class -lclass


