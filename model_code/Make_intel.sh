
### Example build script to compile AWSIM with intel compiler suite
source /opt/intel/bin/compilervars.sh intel64
icc -parallel -par-report3 -O3 -c rk.c ab.c defs.c multigrid.c sor.c pressure.c
icc -parallel -par-report3 -O3 AWSIM.c -o AWSIM.exe -Xlinker rk.o ab.o defs.o multigrid.o sor.o pressure.o
rm *.o
