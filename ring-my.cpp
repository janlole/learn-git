/*
	This is an implementation of a double message passing ring.

	" Implement in c or C++ an MPI program using P processors on a ring (i.e. a simple 1D topology where each 
    processor has a left and right neighbour). The program should implement a stream of messages in both directions:

    - As first step P sends a message ( msgleft = rank ) to its left neighbour (P-1) and receives from its right 
    neighbour (P+1) and send aother message ( msgright = -rank) to P+1 and receive from P-1.
    
    - It then does enough iterations till all processors receive back the initial messages. At each iteration each 
    processor add its rank to the received message if it comes from left, substracting if it comes from right. 
    Both messages originating from a certain processor P should have a tag proportional to its rank (i.e. itag=P*10) "
*/

#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[])
{

	int rank, numproc; // Variables to memorize id processors and number of processors
	int Ldest, Rdest; // Left and Right destination of the messages
	// Ltag, Lmes and Lvalue are associated with the message from right to left (<---)
	// Rtag, Rmes and Rvalue are associated with the message from left to right (--->)
	int Ltag, Rtag;
	int Lmes, Rmes;
	int Lvalue, Rvalue;
	// variables for loop condition and counter
	int check{1};
	int counter{0};
	MPI_Status status;
	
	const int iteration{1000};
	double start_time, initial_time, end_time;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const int first{0};
	const int last{numproc-1};
	
	start_time = MPI_Wtime();

	// inizializing message destinations
	Ldest = (rank != first) ? (rank-1) : last ;
	Rdest = (rank != last) ? (rank+1) : first ;

	// inizializing variables for the first sending
	Rtag = rank*10;
	Ltag = rank*10;

	Rmes = rank;
	Lmes = - rank;

	// first sending
	MPI_Send(&Lmes, 1, MPI_INT, 
		Ldest, Ltag, MPI_COMM_WORLD);
	MPI_Send(&Rmes, 1, MPI_INT, 
		Rdest, Rtag, MPI_COMM_WORLD);

	initial_time = MPI_Wtime();

	for ( auto i{0}; i < iteration; ++i ){
		while(check){
			// tags have to be adjust for receiving and sending
			Rtag = (Rtag > 0) ? (Rtag-10) : (numproc-1)*10 ;
			Ltag = (Ltag < (numproc-1)*10) ? (Ltag+10) : 0 ;

			// from left to right (--->)
			MPI_Recv(&Rvalue, 1, MPI_INT, Ldest, Rtag,
				MPI_COMM_WORLD, &status);
			Rmes = Rvalue + rank ;
			MPI_Send(&Rmes, 1, MPI_INT, Rdest, Rtag,
				MPI_COMM_WORLD);

			// from right to left (<---)
			MPI_Recv(&Lvalue, 1, MPI_INT, Rdest, Ltag,
				MPI_COMM_WORLD, &status);
			Lmes = Lvalue - rank;
			MPI_Send(&Lmes, 1, MPI_INT, Ldest, Ltag,
				MPI_COMM_WORLD);

			// to quit the loop each processor has to receive both from left 
			// and rigth its first sent message, with the tag proportional to
			// its id
			check = (Rtag - rank*10) || (Ltag - rank*10);
			++counter;
		} 
	}

	end_time = MPI_Wtime();

	printf("I am process %d and i have received %d messages. "
		"My final messages have tag %d and left-value %d ,right-value %d\n",
		rank, counter, Rtag, Lvalue, Rvalue);
	printf("I am process %d: initial_time %1.8e total_time %1.8e\n",
		rank, initial_time - start_time,  initial_time - start_time + (end_time - initial_time)/iteration );

	MPI_Finalize();
	return 0;
}