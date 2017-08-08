class:
	cd ../cosmolike_core/class; $(MAKE)

jpl:
	gcc -shared -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib like_fourier.c -o like_fourier.so -fPIC -lgsl -lfftw3 -lgslcblas -std=gnu99 -O3 -ffast-math -funroll-loops -L../cosmolike_core/class -lclass
	gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances compute_covariances.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
	gcc -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o like_fourier like_fourier.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass





