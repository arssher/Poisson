The problem descripton is in doc/manual.pdf
My task number is 6, which means function set 2, uniform grid and max norm.

To compile and run it on BlueGene/P without OpenMP, issue something like
   make clean && MODE=debug COMPILER=ibm make && llsubmit poisson.jcf
Example of poisson.jcf is available in this repo.
Or, simpler, specify the params directly:
   make clean && MODE=release COMPILER=gcc make && mpisubmit.bg -n 128 -m smp --stderr poisson_gcc_2000_128.err poisson -args "0 0 2 2000 1"

'debug' mode sets O0 optimizations and enables logging, see Makefile for the
details.
The program args are described in main.cpp

The graphs were drawn using the commands like
    cd Poisson
    rm -rf res/ && rsync -rav --exclude "core*.data" edu-cmc-stud16-628-06@bluegene.hpc.cs.msu.ru:/home/edu-cmc-stud16-628-06/Poisson/* . && python scripts/visualize.py
What to draw exactly have to be hardcoded inside visualize.py
