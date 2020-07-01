### Example build script to compile AWSIM with gcc, including FFTW library
gcc -O3 -ffast-math -c rk.c ab.c defs.c
gcc -O3 -ffast-math -lm AWSIM.c -o AWSIM.exe -Xlinker rk.o ab.o defs.o -l fftw3
rm *.o


