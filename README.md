# Mater_thesis
Work authored by Francesco Bosia, under the kind supervision of Carsten Magnus.
My master thesis project, tool to simulate HIV Within/between host evolution with mutation and
selection of viral variants.

Executive are created by compiling def.cpp, main.cpp with header.hpp as a header.

an example of compilation could be:

g++ -fopenmp -O3 -std=c++11 main.cpp def.cpp -o CEvo.o

or, in Windows using Cygwin

g++ -fopenmp -O3 -std=c++11 main.cpp def.cpp -o CEvo.exe

execution is regulated by the parameter file parameters.dat or parameters_cluster.dat,
depending on the machine being used. 
Generally, the parameters are tuned to simulate a possible infection of a population of individuals each having a blood volume of 1e-2 L. To scale up/down, fiddle with the parameters.
Most importantly for a standard simulations are the paths parameters and the chunk size.
The paths are important so that the program knows where to go to fish information and print output. 
The chunk size is important to tune the amount of parallelization one wants to attain in the
computation of the within host infection. Very heavy and nasty computation, optimal chunk size is
still not clear.

Output:
path/to/Output/dyn:
- /host_NR.dat contains detailed per-strain information of viral abundances at each time step
- /host_NR_healthy_cells.dat contains general within host dynamics data.
- /Infection_history.dat contains information about between host infections.

path/to/Output/seq:
- /host_NR_seq.dat contains very detailed information about each sequence present at each timestep
and its fitness.

execute with:
./CEvo.o parameters.dat
or
bsub _arguments_ ./CEvo.o parameters_cluster.dat 
on the cluster. For _arguments_ description please refer to the Euler wiki.

Have fun evolving!

Francesco Bosia
