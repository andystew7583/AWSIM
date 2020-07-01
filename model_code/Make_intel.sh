
### Example build script to compile AWSIM with intel compiler suite
source /opt/intel/bin/compilervars.sh intel64
icc -parallel -par-report3 -O3 -c rk.c ab.c defs.c
icc -parallel -par-report3 -O3 AWSIM.c -o AWSIM.exe -Xlinker rk.o ab.o defs.o
rm *.o
