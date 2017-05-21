gcc 4kinFL.c -g -o 4kinFL -I /opt/intel/compilers_and_libraries_2016.1.150/linux/mkl/include/ -I /usr/local/Wolfram/Mathematica/10.3/SystemFiles/IncludeFiles/C/ -L/opt/intel/compilers_and_libraries_2016.1.150/linux/mkl/lib/intel64_lin/ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lm -lpthread -L /opt/intel/compilers_and_libraries_2016.1.150/linux/compiler/lib/intel64_lin/ -liomp5 -O3 -march=native -fopenmp
module load mkl/2016.1.056

#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/compilers_and_libraries_2016.1.150/linux/compiler/lib/intel64_lin/ gdb ./4kin4_1 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/compilers_and_libraries_2016.1.150/linux/compiler/lib/intel64_lin/ 

#one run:
./4kinFL path/hamiltonian.dat Efermi k1 k2 k3

#reaction path

for ((i=n;i<N;i++)); do
./4kinFL path/hamiltonian${i}_hr.dat Efermi k1 k2 k3
done
#time
