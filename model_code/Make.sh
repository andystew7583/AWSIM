### Example build script to compile AWSIM with gcc
gcc -O3 -ffast-math -c rk.c ab.c defs.c multigrid.c sor.c pressure.c
gcc -O3 -ffast-math -lm AWSIM.c -o AWSIM.exe -Xlinker rk.o ab.o defs.o multigrid.o sor.o pressure.o
rm *.o
