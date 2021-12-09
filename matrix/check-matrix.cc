	// NOTE: Starting 1D decomposition
/*
	my 3D matrix is rows*cols*mats matrix
	so the 1D decomposition has to divide the number of mats per processor
		the 2D decompotition has to divede the mats by rows (the element per column are contiguos)

	THIS IMPLEMENTATION ALLOWS TO COMPUTE AN ADDITION OF Nr x Nc x Nm MATRICES
	ONLY WHEN Nr x Nc x Nm IS MULTIPLE OF THE TOTAL NUMBER OF PROCESSES

	THE AIM OF THIS PROGRAM IS JUST TO OBTAIN A BASELINE FOR matrix.cc
*/

#include <iostream>
#include <vector>
#include <random>
#include "mpi.h"

#define MAX_DIM_DEC 3 													// maximum number of dimensions on which make decomposition
#define LIMIT 100'000

int main(int argc, char *argv[])
{
	const int root{0};
	int numproc,rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numproc);

	// __________________________________________________________________
	// inizializing matrices dimensions
	
	int dim_vect[MAX_DIM_DEC];
	dim_vect[0] = std::stoi(argv[1]);
	dim_vect[1] = std::stoi(argv[2]);
	dim_vect[2] = std::stoi(argv[3]);
	int dimension{dim_vect[0]*dim_vect[1]*dim_vect[2]};
	//__________________________________________________________________
	// Inizializing the two matrices

	std::vector<double> matA (dimension);						
	std::vector<double> matB (dimension);						
	std::vector<double> matResult (dimension);						
	std::vector<double> matChek (dimension);						
	std::uniform_real_distribution<double> unif(-LIMIT,LIMIT);
	std::default_random_engine re;
	if (rank == root){	
		// fastest way
		for (auto i{0}; i < dimension; ++i){
			matA[i] =  unif(re);
			matB[i] =  unif(re);
			matChek[i] = matA[i] + matB[i];
		}
	}	
	//__________________________________________________________________
	// timing the pre-sending operations
	double start_pre_send;
	start_pre_send = MPI_Wtime();
	//__________________________________________________________________	
	// Declaring two vectors for the submatrices				


	int elems{dimension/numproc};
	
	std::vector<double> subA (elems);
	std::vector<double> subB (elems);
	//__________________________________________________________________
	// timing the pre-sending operations
	double pre_send;
	pre_send = MPI_Wtime() - start_pre_send;
	//__________________________________________________________________
	// Sending to each processor the main decomposition

	MPI_Scatter(&matA[0], elems, MPI_DOUBLE, 
		&subA[0], elems, MPI_DOUBLE,
		root, MPI_COMM_WORLD);
	MPI_Scatter(&matB[0], elems, MPI_DOUBLE, 
		&subB[0], elems, MPI_DOUBLE,
		root, MPI_COMM_WORLD);
	//__________________________________________________________________
	// each processor execute the computation on its data

	for (auto i{0}; i < elems; ++i){
		subA[i] += subB[i];
	}
	//__________________________________________________________________
	// collecting the results from all the processors, main decomposition

	MPI_Gather(&subA[0], elems, MPI_DOUBLE, 
		&matResult[0] , elems, MPI_DOUBLE,
			root, MPI_COMM_WORLD );
	//__________________________________________________________________
	// timing the computation of the main part
	double main_matrix;
	main_matrix = MPI_Wtime() - pre_send;
	//__________________________________________________________________

	// checking the results
	int check_result{0};
	if (rank == root){
		for (auto i{0}; i < dimension; ++i){
			if (matResult[i] - matChek[i]){
				check_result = 1;
				break;
			}
		}
		if (check_result){
			std::cout << "----- WARNING ----- \n"
					  << "I am process " << rank
					  << " and something went WRONG\n" 
					  << "----- WARNING ----- \n"
					  << std::endl;
		}
		else{		
			std::cout << "I am process " << rank 
					  << " and the computation went fine" << std::endl;
			std::cout << "\n--- TIME ---\n"
					  << "pre_send\tmain_matrix\ttotal\n"
					  << pre_send <<'\t'<< main_matrix
					  << '\t' <<pre_send+main_matrix
					  << std::endl;
		}
	}
	//__________________________________________________________________


	MPI_Finalize();
	return 0;
}
