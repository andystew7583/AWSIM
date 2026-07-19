### Example build script to compile AWSIM with gcc, including FFTW library
gcc -O3 -ffast-math -c rk.c ab.c defs.c multigrid.c sor.c pressure.c
gcc -O3 -ffast-math -lm AWSIM.c -o AWSIM.exe -Xlinker rk.o ab.o defs.o multigrid.o sor.o pressure.o -l fftw3
rm *.o
