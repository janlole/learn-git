#include <stdio.h>
#include "mpi.h"
#include <iostream>

#define ITERATION 100'000

int main(int argc, char *argv[])
{
	int numproc, rank;
	double start_time, end_time, initial_time;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	initial_time = MPI_Wtime();

	const int first{0};
	const int last{numproc-1};
	
	int counter;
	MPI_Status status;

	// inizializing message destinations
	int Rdest, Ldest;
	Ldest = (rank != first) ? (rank-1) : last ;
	Rdest = (rank != last) ? (rank+1) : first ;
	
	int Rvalue, Lvalue;
	int Rtag, Ltag;
	
	start_time = MPI_Wtime();

	for (auto i{0}; i < ITERATION; ++i){
		// inizializing variables for the first sending
		Rtag = rank*10;
		Ltag = rank*10;
		Rvalue = 0;
		Lvalue = 0;
		counter = 0;
		do{
	
			// sending to right (--->)
			Rvalue += rank ;
			MPI_Send(&Rvalue, 1, MPI_INT, Rdest, Rtag,
				MPI_COMM_WORLD);
			// sending to left (<---)
			Lvalue -= rank;
			MPI_Send(&Lvalue, 1, MPI_INT, Ldest, Ltag,
				MPI_COMM_WORLD);
	
			// adjusting tag to receive, Rtag is the previous one
			// receiving from left (--->)
			Rtag = (Rtag > first) ? (Rtag-10) : last*10 ;
			MPI_Recv(&Rvalue, 1, MPI_INT, Ldest, Rtag,
				MPI_COMM_WORLD, &status);
	
	
			// adjusting tag to receive, Ltag is the next one
			// receiving from right (<---)
			Ltag = (Ltag < last*10) ? (Ltag+10) : first;
			MPI_Recv(&Lvalue, 1, MPI_INT, Rdest, Ltag,
				MPI_COMM_WORLD, &status);
	
			++counter;

			// to quit the loop each processor has to receive both from left 
			// and rigth its first sent message, with the tag proportional to
			// its id
	
		} while( (Rtag - rank*10 ) || (Ltag-rank*10));
	}
	end_time = MPI_Wtime();

	std::cout << "I am process " << rank
			  << " and i have received " << counter << " messages."
			  << "My final messages have tag " << Rtag 
			  << " and left-value " << Lvalue << " ,right-value " << Rvalue 
			  << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == first){
		std::cout << "---TIME---\n"
				  << "numproc\texec-time\tinit-time\n"
				  << numproc << "\t"
				  << (end_time - start_time) / ITERATION << "\t"
				  << (start_time - initial_time) 
				  << std::endl;
		}
	MPI_Finalize();
	return 0;
}
/*
 *	std::cout << "TIME : "
				  << "We were " << numproc << " and the ring execution time was "
				  << (end_time - start_time) / ITERATION 
				  << " and the initializaton time was "
				  << (start_time - initial_time) 
				  << std::endl;

 *
 *
 * */
