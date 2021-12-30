# README: section 1
All the programs were compiled with the last openmpi library and the command used to load it is
```sh
module load openmpi-4.1.1+gnu-9.3.0
```
# Ring
Compile the ring with the command
```sh
mpic++ ring.cc -o ring.o
```
and tun it with the line
```sh
mpirun --map-by node -np num ./ring.o
```
where the *num* is the number of processes used by the ring (it will use all the processes available). The *--map-by node* flag were used to have comparable results when it were run on more than one node.
# Matrix
Compile the matrix with the command
```sh
mpic++ sum3Dmatrix.cc -o sum3Dmatrix.o
```
while 
```sh
mpirun -np 24 ./sum3Dmatrix.o Grid_one Grid_two Grid_three Matrix_one Matrix_two Matrix_three
```
where the *Grid* are the grid of processes used for the run and the *Matrix* arguments are the matrix dimensions. So for example to add two random matrices of 1200x200x100 with a 3D decomposition on the grid 3x2x4 the following line has to be used:
```sh
mpirun -np 24 ./sum3Dmatrix.o 3 2 4 1200 200 100
```