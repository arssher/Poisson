The problem descripton is in doc/manual.pdf
My task number is 6, which means function set 2, uniform grid and max norm.

To compile and run it on BlueGene/P without OpenMP, issue something like
   make clean && MODE=debug COMPILER=ibm make && llsubmit poisson.jcf
Example of poisson.jcf is available in this repo.
Or, simpler, specify the params directly:
   make clean && MODE=release COMPILER=gcc make && mpisubmit.bg -n 128 -m smp --stderr poisson_gcc_2000_128.err poisson -args "0 0 2 2000 1"
Building with OpenMP support:
   make clean && MODE=release OMP=t make

To compile and run it on Lomonosov, first of all issue
   module add slurm/2.5.6

Then, if you want to compile the program with gcc+openmpi, run
   module add openmpi/1.8.4-gcc
Or, if you prefer icpc+impi,
   module add intel/15.0.090
   module add impi

After that mpicxx will point to the needed compiler.
Run the program like that:
   sbatch -p regular4 -n128  --time=0-00:2:00 --error=poisson_gcc_2000_128.err ompi ./poisson 0 0 2 2000 1
Replace 'ompi' with 'impi' if you use Intel stuff.

'debug' mode sets O0 optimizations and enables logging, see Makefile for the
details.
The program args are described in main.cpp

The graphs were drawn using the commands like
    cd Poisson
    rm -rf res/ && rsync -rav --exclude "core*.data" edu-cmc-stud16-628-06@bluegene.hpc.cs.msu.ru:/home/edu-cmc-stud16-628-06/Poisson/* . && python scripts/visualize.py
What to draw exactly have to be hardcoded inside visualize.py
