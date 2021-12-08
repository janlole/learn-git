// NOTE: Starting 1D decomposition
/*
	my 3D matrix is rows*cols*mats matrix
	so the 1D decomposition has to divide the number of mats per processor
		the 2D decompotition has to divede the mats by rows (the element per column are contiguos)
	to do:
		- find a way to partitionate the 24 cores in a friendly-matrix-dimensions way

*/

#include <iostream>
#include <vector>
#include <random>
#include "mpi.h"

#define MAX_DIM_DEC 3 													// maximum number of dimensions on which make decomposition
#define LIMIT 100'000

void round_robin_distribution(int remainder, int cut, int square, int skip, int distance,
	int start_index, 
	const std::vector<double>& matA, const std::vector<double>& matB, std::vector<double>& matResult,
	MPI_Comm communicator);

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
	dim_vect[0] = std::stoi(argv[4]);
	dim_vect[1] = std::stoi(argv[5]);
	dim_vect[2] = std::stoi(argv[6]);
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
	// constructing virtual topology dimensions

	int dims [MAX_DIM_DEC];
	dims[0] = std::stoi(argv[1]);
	dims[1] = std::stoi(argv[2]);
	dims[2] = std::stoi(argv[3]);
	
	int block[MAX_DIM_DEC];
	int remainder[MAX_DIM_DEC];
	for (auto i{0}; i < MAX_DIM_DEC; ++i){
		block[i] = dim_vect[i]/dims[i];
		remainder[i] = dim_vect[i]%dims[i];
	}
	//__________________________________________________________________
	// Constructing blocks with Type_vector

	MPI_Datatype tmpblock;												// temporary datatype for block construction
	MPI_Datatype block_square;											// a square-shape portion of a 2D matrix
	MPI_Type_vector(block[1], block[0], dim_vect[0], MPI_DOUBLE, &tmpblock);
	MPI_Type_create_resized( tmpblock, 0, sizeof(double), &block_square);
	MPI_Type_commit(&block_square);

	MPI_Datatype block_cube;											// a cube-shape portion of a 3D matrix
	MPI_Type_vector( block[2], 1 ,dim_vect[0]*dim_vect[1], block_square, &tmpblock);
	MPI_Type_create_resized( tmpblock, 0, sizeof(double), &block_cube);
	MPI_Type_commit(&block_cube);
	//__________________________________________________________________
	// Placement of the blocks in which the matrix should be split

	int disps[numproc];
	int count[numproc];

	for (auto k{0}; k < dims[2]; ++k){
		for (auto j{0}; j < dims[1]; ++j){
			for (auto i{0}; i < dims[0]; ++i){
				disps[ k * dims[0]*dims[1]+ j*dims[0] + i] = 
				k * dim_vect[0] * dim_vect[1] * block[2]
				+ j * dim_vect[0] * block[1] 
				+ i * block[0] ;
				count[ k * dims[0]*dims[1]+ j*dims[0] + i] = 1;
				
			}
		}
	}
	//__________________________________________________________________
	// Declaring two vectors for the submatrices				

	int elems{block[0]*block[1]*block[2]};
	
	std::vector<double> subA (elems);
	std::vector<double> subB (elems);
	//__________________________________________________________________
	// timing the pre-sending operations
	double pre_send;
	pre_send = MPI_Wtime() - start_pre_send;
	//__________________________________________________________________
	// Sending to each processor the main decomposition

	MPI_Scatterv(&matA[0], count, disps, block_cube, 
		&subA[0], elems, MPI_DOUBLE,
		root, MPI_COMM_WORLD);
	MPI_Scatterv(&matB[0], count, disps, block_cube, 
		&subB[0], elems, MPI_DOUBLE,
		root, MPI_COMM_WORLD);
	//__________________________________________________________________
	// each processor execute the computation on its data

	for (auto i{0}; i < elems; ++i){
		subA[i] += subB[i];
	}
	//__________________________________________________________________
	// collecting the results from all the processors, main decomposition

	MPI_Gatherv(&subA[0], elems, MPI_DOUBLE, 
		&matResult[0] , count, disps, block_cube,
			root, MPI_COMM_WORLD );
	//__________________________________________________________________
	// timing the computation of the main part
	double main_matrix;
	main_matrix = MPI_Wtime() - pre_send;
	//__________________________________________________________________
	// Starting a round robin distribution of the remainders

	if (remainder[0]+remainder[1]+remainder[2]){

		if  (remainder[2]  ){
			round_robin_distribution(remainder[2], dim_vect[1], dim_vect[0], 
				dim_vect[0] , dim_vect[0]*dim_vect[1],
				dim_vect[0]*dim_vect[1]*(dim_vect[2]-remainder[2]), 
				matA,matB, matResult, MPI_COMM_WORLD);
		}

		if (remainder[1] ){
			round_robin_distribution(remainder[1], dim_vect[2]-remainder[2] , dim_vect[0], 
				dim_vect[0]*dim_vect[1] , dim_vect[0],
				dim_vect[0]*(dim_vect[1]-remainder[1]), 
				matA,matB, matResult,  MPI_COMM_WORLD);
		}
		
		if (remainder[0]  ){
			round_robin_distribution( dim_vect[1]-remainder[1] , dim_vect[2]-remainder[2], remainder[0], 
				dim_vect[0]*dim_vect[1], dim_vect[0],
				dim_vect[0]-remainder[0], 
				matA, matB, matResult, MPI_COMM_WORLD);
		}
	}
	//__________________________________________________________________
	// timing the computation of the remaining part
	double remain_matrix;
	remain_matrix = MPI_Wtime() - main_matrix;
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
					  << "pre_send\tmain_matrix\tremain_matrix\ttotal.comput\ttotal\n"
					  << pre_send <<'\t'<< main_matrix<<'\t'<<remain_matrix
					  << '\t' << main_matrix + remain_matrix
					  << '\t' <<pre_send+main_matrix+remain_matrix
					  << std::endl;
		}
	}
	//__________________________________________________________________


	MPI_Finalize();
	return 0;
}

void round_robin_distribution(int remainder, int cut, int square, int skip, int distance,
	int start_index, const std::vector<double>& matA, const std::vector<double>& matB, std::vector<double>& matResult,
	MPI_Comm communicator){

	MPI_Datatype tmpblock;
	MPI_Datatype block_remainder;										// block type for round robin distribution
	MPI_Type_vector( remainder, square, distance, MPI_DOUBLE, &tmpblock);
	MPI_Type_create_resized( tmpblock, 0, sizeof(double), &block_remainder);
	MPI_Type_commit(&block_remainder);

	int rank;
	int how_many_proc;
	MPI_Status status;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &how_many_proc);
	--how_many_proc;

	int elems {remainder*square};
	std::vector<double> subA(elems);
	std::vector<double> subB(elems);

	if (rank == 0){
		int start{1}, end;
		int check{cut/how_many_proc + (0 < cut%how_many_proc)};
		int max_row{(cut > how_many_proc)? how_many_proc : cut};
		int rem_row{( cut%how_many_proc )? cut%how_many_proc : max_row };

		for (auto j{0}; j < check ; ++j){
			int trasl{start_index + j * skip * how_many_proc};
			end = ( j < check - 1) ? max_row : rem_row ;
			start = 1;
			for (auto i{0}; i < end; ++i){
				MPI_Send(&matA[ trasl + i*skip], 1, block_remainder, 
					start+i, 0, communicator );
				MPI_Send(&matB[ trasl + i*skip], 1, block_remainder, 
					start+i, 0, communicator );
			}
			for (auto i{0}; i < end; ++i){
				MPI_Recv(&matResult[ trasl + i*skip], 1, block_remainder, 
					start+i, 0, communicator, &status );
			}
		}
	}
	else {
		for (auto j{0}; j < cut/how_many_proc + ( rank - 1 < cut%how_many_proc); ++j){
			MPI_Recv(&subA[0], elems, MPI_DOUBLE,
				0, 0, communicator, &status);
			MPI_Recv(&subB[0], elems, MPI_DOUBLE,
				0, 0, communicator, &status);

			for (auto i{0}; i < elems; ++i)
				subA[i] += subB[i] ;
			MPI_Send(&subA[0], elems, MPI_DOUBLE,
				0, 0, communicator);
		}			
	}
}
	/*
	// printing the overall result

	if (rank == root){
		std::cout << "final result\n";
		for (auto z{0}; z <dim_vect[2]; ++z){
			for (auto y{0}; y < dim_vect[1]; ++y){
				for (auto x{0}; x < dim_vect[0]; ++x){
					std::cout << matResult[ x + y * dim_vect[0] 
						+ z *dim_vect[0]*dim_vect[1]] 
					<< '\t';
				}
				std::cout << '\n';
			}
			std::cout << "\n\n";
		}
		std::cout  <<std::endl;
	}
	//__________________________________________________________________
	/**/