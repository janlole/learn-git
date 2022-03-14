#include <iostream>
#include <vector>
#include <random>
#include <memory>
#include <cmath>
#include <chrono>
#include "mpi.h"


#if !defined(NDIM)
#define NDIM 2
#endif

#define LIMIT 1

#if !defined(NUMPOINTS)
#define NUMPOINTS 4095
#endif

#if !defined(DOUBLE_PRECISION_KPOINT)
#define float_t float
#define MPI_float_t MPI_FLOAT
#else
#define float_t double
#define MPI_float_t MPI_DOUBLE
#endif


struct kpoint {
	float_t coord[NDIM];

	kpoint() : coord{} {}
	~kpoint() noexcept {}

	const float_t* begin() const { return &coord[0]; }
	float_t* begin() { return &coord[0]; }
	
	const float_t* end() const { return &coord[NDIM]; }
	float_t* end() { return &coord[NDIM]; }

	float_t& operator[](const unsigned int i) noexcept { return coord[i]; }
	const float_t& operator[](const unsigned int i) const noexcept { return coord[i]; }
	
	// Move ctor and assignment
	kpoint(kpoint&& x) noexcept = default;
	kpoint& operator=(kpoint&& x) noexcept = default;
	//_____________________________________________
	
	// Copy ctor and assignment:  the ORIGINAL copy ctor is NOT exception safe
	kpoint(const kpoint& v) : coord{} {
		std::copy(v.begin(), v.end(), this->begin());
	}
	kpoint& operator=(const kpoint& x){ // if acquiring resources DO NOT mark "noexcept"
		// coord.reset();
		auto tmp = x;
		(*this) = std::move(tmp);
		return *this;
	}
	//_____________________________________________

	bool operator==( kpoint& l){
		for ( auto i{0} ; i < NDIM ; ++i ){
			if ( l.coord[i] != this -> coord[i] ){
				return false;}
		}
		return true;
	}
	bool operator!=( kpoint& l){
		return !( *this == l);
	}
	//_____________________________________________
	// swap method
	void swap(kpoint& p){
		auto tmp {std::move(p)};
		p = std::move(*this);
		*this = std::move(tmp);
	}
	//_____________________________________________

};

struct knode
{
	int axis;						// splitting dimension
	kpoint split;					// splittting element
	struct knode *left, *right;		// left and right subtrees

	knode(): axis{}, split{}, left{nullptr}, right{nullptr} {}

	bool operator==( knode& l){
		if ((l.axis != this -> axis) | (l.split != this->split))
			return false;
		if ( (l.left == nullptr) ^ (this -> left == nullptr) )
			return false;
		if ( (l.right == nullptr) ^ (this -> right == nullptr) )
			return false;
		if ( (l.left != nullptr) && (this -> left != nullptr) && (*(l.left) != *(this->left)) )
			return false;
		if ( (l.right != nullptr) && (this -> right != nullptr) && (*(l.right) != *(this->right)) )
			return false;
		return true;

	}

	bool operator!=(knode& l){
		return !(*this == l);
	}
};


kpoint* partition( const float_t median, const int axis, kpoint* low, kpoint* high);

kpoint* core_algorithm(kpoint* low, kpoint* high, knode* node, int axis);

typedef kpoint* (COMP)(kpoint*,kpoint*, knode*, int);

template<COMP choosen_algorithm>
knode* build_kdtree_recursive(kpoint* low, kpoint* high, knode* node, int depth);

template<COMP choosen_algorithm>
kpoint* build_one_knode_mpi(kpoint* first_kpoint, kpoint* last_kpoint, knode* node, int axis);

void sorting(const int axis, kpoint* low, kpoint* high);

bool check_sorting( const int axis, kpoint* low, kpoint* high);

kpoint* select(kpoint* start, kpoint* end, const int position, const int axis);

// function usefull to control the quality of the kdtrees

int height(knode* node);
int count_subtree(knode* node, int height);
bool check_kdtree(knode* node);
knode* find_kpoint_fast(knode* node, kpoint target);


int main(int argc, char *argv[])
{
	const int root{0};
	int numproc,rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numproc);
	
	// MPI_KPOINT AND MPI_KNODE DATATYPE
	struct kpoint dummy_kpoint;
	MPI_Datatype MPI_KPOINT;
	int blocklens[NDIM];
	MPI_Datatype old_type[NDIM];
	MPI_Aint indices[NDIM];
	MPI_Aint base_address_kpoint;
	MPI_Get_address( &dummy_kpoint, &base_address_kpoint);
	for ( auto i{0}; i < NDIM;++i){
		blocklens[i] = 1;
		old_type[i] = MPI_float_t;
		MPI_Get_address( &dummy_kpoint.coord[i], &indices[i]);
		indices[i] = MPI_Aint_diff(indices[i], base_address_kpoint);
	}
	MPI_Type_create_struct(NDIM, blocklens, indices, old_type, &MPI_KPOINT );
	MPI_Type_commit(&MPI_KPOINT);

	struct knode dummy_knode;
	MPI_Datatype MPI_KNODE;
	int blocklens_knode[4] = {1,1,1,1};
	MPI_Aint indices_knode[4];
	MPI_Aint base_address;
	MPI_Get_address(&dummy_knode, &base_address);
	MPI_Get_address(&dummy_knode.axis, &indices_knode[0]);
	MPI_Get_address(&dummy_knode.split, &indices_knode[1]);
	MPI_Get_address(&dummy_knode.left, &indices_knode[2]);
	MPI_Get_address(&dummy_knode.right, &indices_knode[3]);
	indices_knode[0] = MPI_Aint_diff(indices_knode[0], base_address);
	indices_knode[1] = MPI_Aint_diff(indices_knode[1], base_address);
	indices_knode[2] = MPI_Aint_diff(indices_knode[2], base_address);
	indices_knode[3] = MPI_Aint_diff(indices_knode[3], base_address);
	MPI_Datatype old_type_knode[4] = {MPI_INT, MPI_KPOINT, MPI_AINT, MPI_AINT};
	MPI_Type_create_struct(4, blocklens_knode, indices_knode, old_type_knode, &MPI_KNODE);
	MPI_Type_commit(&MPI_KNODE);
	
	//_____________________________________________

	// to be more precise: the child cores need less space than the parent ones (exactly (parent's kpoint - 1)/2)
	// and it is decided from the start, so here must be already set
	// BUT due to the fact that actually there could be some problem when NUMPOINTS is not a power of two, 
	// I decided to keep the size one ore two kpoints larger
	int n_step{static_cast<int>(std::floor(log2(numproc)))};
	int sub_tree_size{0};
	for ( auto step{0}; step < n_step; ++step){
		if (rank < pow(2,step+1) && rank >= pow(2,step)){
			sub_tree_size = static_cast<int>(std::ceil(NUMPOINTS/(pow(2,step+1))));
			// std::cout << "I'm processor\t" << rank << "\tand my space is\t" << sub_tree_size << std::endl;
		}
	}
	if ( rank == root )
		sub_tree_size = NUMPOINTS;

	// std::unique_ptr<kpoint[]> Grid{new kpoint[sub_tree_size]};
	// std::unique_ptr<knode[]> Nodes{new knode[sub_tree_size]};

	std::vector<kpoint> Grid(sub_tree_size);
	std::vector<knode> Nodes(sub_tree_size);


	// Random assignment and tree of control
	MPI_Status status;

	if ( rank == root ){
		std::uniform_real_distribution<float_t> unif(-LIMIT,LIMIT);
		std::default_random_engine re{static_cast<long unsigned int>(time(0))};
		for ( auto x{0}; x < sub_tree_size; ++x ){
			for ( auto& y : Grid[x]){
				y = unif(re);
			}
		}
	}

	//_____________________________________________

	double start_time, end_time, start_building, end_building;
	start_time = MPI_Wtime();

	kpoint* first_kpoint{&Grid[0]};
	kpoint* mide;
	kpoint* last_kpoint{&Grid[NUMPOINTS-1]}; 
	knode* working_node{&Nodes[0]};
	knode* relink_node; knode* relink_node_received; knode* tmp_knode_ptr_1; 
	int axis{0}, process_message, size_message{NUMPOINTS};

	// step is the number of sending messages and 2^(step+1) is total number of process
	for (auto step{0}; step < n_step; ++step){
		if (rank < pow(2,step) ){
			// construct one node with saved information
			axis = (axis + 1)%NDIM; 
			mide = build_one_knode_mpi<core_algorithm>(first_kpoint, last_kpoint, working_node, axis);
			// send instruction for build next node to the (rank + 2^step) core
			process_message = rank + pow(2,step);
			MPI_Send((mide+1), (last_kpoint - mide) , MPI_KPOINT, process_message,0, MPI_COMM_WORLD );
			// axis is needed only if there is no WIDDE_AXIS activate (note WIDE_AXIS must be set yet)
			MPI_Send(&axis, 1, MPI_INT, process_message, 0, MPI_COMM_WORLD); 
			// pointer to the node that will be construct by the other process 
			// (that will be usefull to "re-link" togheter the substrees)
			MPI_Send(&working_node, 1, MPI_AINT, process_message, 0, MPI_COMM_WORLD); 
			// save information to reconstruct the kdtree
			working_node = working_node -> left;
			last_kpoint = mide - 1;
		}
		else if (rank < pow(2,step+1) && rank >= pow(2,step)){
			process_message = rank - pow(2,step);
			// receive information from (rank - 2^step) core
			MPI_Probe(process_message, 0, MPI_COMM_WORLD, &status);
		    MPI_Get_count(&status, MPI_KPOINT, &size_message);
			MPI_Recv(&Grid[0], size_message , MPI_KPOINT, process_message, 0 , MPI_COMM_WORLD, &status);
			MPI_Recv(&axis, 1, MPI_INT, process_message, 0, MPI_COMM_WORLD, &status); 
			// SAVE THE POINTER SOMEWHERE TO RELINK THE SUBTREES
			MPI_Recv(&relink_node, 1, MPI_AINT, process_message, 0, MPI_COMM_WORLD, &status); 
			// save information for the job( arguments of core algorithm function)
			first_kpoint = &Grid[0];
			last_kpoint = &Grid[0] + size_message - 1;
		}
	}

	start_building = MPI_Wtime();

	build_kdtree_recursive<core_algorithm>(first_kpoint, last_kpoint, working_node, axis+1);
	working_node = working_node + ( last_kpoint - first_kpoint + 1) ;

	end_building = MPI_Wtime();

	// step is the number of sending messages and 2^(step+1) is total number of process
	for (auto step{n_step-1}; step >= 0; --step){
		if (rank < pow(2,step+1) && rank >= pow(2,step)){
			process_message = rank - pow(2,step);
			// send subtree and the link address to (rank - 2^step) core
			MPI_Send(&Nodes[0], working_node - &Nodes[0], MPI_KNODE, process_message, 0 , MPI_COMM_WORLD);
			MPI_Send(&relink_node, 1, MPI_AINT, process_message, 0, MPI_COMM_WORLD); 
			relink_node = &Nodes[0];
			MPI_Send(&relink_node, 1, MPI_AINT, process_message, 0, MPI_COMM_WORLD); 

		}
		else if (rank < pow(2,step) ){
			process_message = rank + pow(2,step);
			// send instruction for build next node to the (rank + 2^step) core
			MPI_Probe(process_message, 0, MPI_COMM_WORLD, &status);
		    MPI_Get_count(&status, MPI_KNODE, &size_message);
			MPI_Recv(working_node, size_message, MPI_KNODE, process_message,0, MPI_COMM_WORLD , &status);
			MPI_Recv(&relink_node_received, 1, MPI_AINT, process_message, 0, MPI_COMM_WORLD, &status); 

			// link process
			relink_node_received -> right = working_node;
			// re-arrange addresses of the committed/received subtree

			MPI_Recv(&relink_node_received, 1, MPI_AINT, process_message, 0, MPI_COMM_WORLD, &status); 

			knode* tmp_knode_ptr{working_node};
			while ( tmp_knode_ptr < working_node + size_message ){
				if ( tmp_knode_ptr -> left != nullptr )
					tmp_knode_ptr -> left = tmp_knode_ptr->left - relink_node_received + working_node;
				if ( tmp_knode_ptr -> right != nullptr )
					tmp_knode_ptr -> right = tmp_knode_ptr->right - relink_node_received + working_node;
				++tmp_knode_ptr;
			}
			working_node = tmp_knode_ptr ;
		}
	}

	end_time = MPI_Wtime();

	if (rank == root){
		std::cout << "NDIM,NUMPOINTS,numproc,MPI TIME,sending_time,building_time,receiving_time\n" 
				  << NDIM << ","									// ndim
				  << NUMPOINTS << ","								// numpoints
				  << numproc << ","									// numproc
				  << end_time - start_time << ","					// total-time= MPI_TIME
				  << start_building - start_time << ","				// sending_time
				  << end_building - start_building << ","			// building_time
				  << end_time - end_building 						// receiving_time
				  << std::endl;
	}




	#if defined(CHECK_CORRECTNESS)
		if(rank == root){
					std::vector<knode> Nodes_copy(NUMPOINTS);
					build_kdtree_recursive<core_algorithm>(&Grid[0],(&Grid[NUMPOINTS])-1, &Nodes_copy[0], 1);
					knode* research;
					int count_not_found{0};
					int num_point{(std::min(NUMPOINTS, int(pow(2,11) -1) ))};
					for ( auto kpoint_target {0}; kpoint_target < num_point; ++kpoint_target ){
						research = find_kpoint_fast(&Nodes[0], Grid[kpoint_target]);
						if ( research == nullptr ){
							++count_not_found;
						}
						research = find_kpoint_fast(&Nodes[0], Grid[NUMPOINTS - 1 - kpoint_target]);
						if ( research == nullptr ){
							++count_not_found;
						}
					}
					
					std::cout << "\n\n############\n\n";
					std::cout << "sizeof(float_t):\t" << sizeof(float_t) << std::endl;
					std::cout << "NUMBER OF POINTS WHICH WILL BE CHECKED\t" << num_point*2  << "\ton\t" << NUMPOINTS << "\n" <<std::endl;
					if (count_not_found){
						std::cout << "EXISTENCE TEST\t FAILED\n"
								  << "number of points not found\t" << count_not_found << std::endl;
					}
					else
						std::cout << "EXISTENCE TEST\t PASSED\n";
		
					std::cout << "\n\n###########\n\n";
		
					int h{height(&Nodes[0])};
					int h_sort{height(&Nodes_copy[0])};
					int teoretic_height{static_cast<int>(std::ceil(log2(NUMPOINTS+1)-1))};
					if ( (h - teoretic_height) | (h_sort - teoretic_height)){
						std::cout << "HEIGHT TEST\t FAILED\n"
								  << "teoretic height:\t" << teoretic_height << "\theight:\t" << h << "\theight_copy:\t" << h_sort << std::endl;
					}
					else
						std::cout << "HEIGHT TEST\t PASSED\n";
		
					std::cout << "\n\n###########\n\n";
					
					bool mammamia(Nodes[0]  == Nodes_copy[0]);
					std::cout << "COMPARISON\t" << mammamia << std::endl;
		
					std::cout << "\n\n###########\n\n";
					
					bool equili(check_kdtree(&Nodes[0]));
					bool equili_sort(check_kdtree(&Nodes_copy[0]));
					if ( !equili | !equili_sort){
						std::cout << "EQUILIBRIUM TEST\t FAILED\n"
								  << "EQUILIBRIUM Nodes\t" << equili << "\tEQUILIBRIUM_copy\t" << equili_sort << std::endl;
					}
					else
						std::cout << "EQUILIBRIUM TEST\t PASSED\n";

					std::cout << "\n\n###########\n\n";
					int tot_nod{ count_subtree(&Nodes[0], h)};
					int tot_nod_co{count_subtree(&Nodes_copy[0], h_sort)};
					
					if ( (NUMPOINTS - tot_nod) | (NUMPOINTS - tot_nod_co)){
						std::cout << "NUMBER KPOINT TEST\t FAILED\n"
								  << "TOTAL KPOINTS Nodes\t" << tot_nod << "\tTOTAL KPOINTS_copy\t" << tot_nod_co << "\tTOTAL original KPOINTS \t" << NUMPOINTS << std::endl;
					}
					else
						std::cout << "NUMBER KPOINT TEST\t PASSED\n";
		
					std::cout << "\n\n###########\n\n";
		}
	#endif

	return 0;
}
template<COMP choosen_algorithm>
kpoint* build_one_knode_mpi(kpoint* first_kpoint, kpoint* last_kpoint, knode* node, int axis){
	knode* last_node{node};
	kpoint* mide; 
	if ( first_kpoint == last_kpoint){
		node -> split = *first_kpoint;			
		node -> axis = 0;
	}
	else{
		mide = choosen_algorithm(first_kpoint,last_kpoint,node,axis);
		if ( first_kpoint != mide ){
			++last_node;
			last_node -> axis = NDIM ;			
			node -> left = last_node;
		}
	}
	return mide;
}

template<COMP choosen_algorithm>
knode* build_kdtree_recursive(kpoint* low, kpoint* high, knode* node, int depth){
	if ( low == high ){
		// if one point is available, return a leaf without any computation
		node -> split = *low;
		return node;		
	}
	int axis{depth%NDIM};	
	
	knode* node_working{node};
	kpoint* mide{choosen_algorithm(low,high,node,axis)};
	
	++depth;
	kpoint* tmp;
	if ( low != mide ){
		++node_working;
		tmp = mide-1;
		node -> left = build_kdtree_recursive<choosen_algorithm>(low, tmp, node_working, depth);
	}
	if ( high != mide ){
		node_working = node + (mide - low) + 1;
		tmp = mide+1;
		node -> right = build_kdtree_recursive<choosen_algorithm>(tmp, high, node_working, depth);
	}
	return node;
}


kpoint* core_algorithm(kpoint* low, kpoint* high, knode* node, int axis){
	kpoint* mide{select(low, high, int((high-low)/2), axis ) };

	node -> split = *(mide);
	node -> axis = axis;
	return mide;
}


kpoint* select(kpoint* start, kpoint* end, const int position, const int axis){
	if ( start == end )
		return start;
	if ( end - start < 5){
		sorting(axis, start, end);
		return start + position; 
	}

	int size{(int(end-start)/5) + (int(end-start)%5>0)};
	std::vector<kpoint> tmp (size) ;
	kpoint* start_tmp{start};
	kpoint* end_tmp{end};
	int start_index{0};
	while ( start_index < (int(end-start)/5) ){
		sorting(axis, start_tmp, start_tmp + 4);
		tmp[start_index] = *(start_tmp+2);
		++start_index;
		start_tmp += 5;
	}
	if ( start_index < (int(end-start)/5) +  (int(end-start)%5>0)){
		sorting(axis, start_tmp, end_tmp);
		tmp[start_index] = *((int(end - start_tmp))/2 + start_tmp);
	}

	kpoint* mide{select(&tmp[0], &tmp[size-1], int(size/2), axis)};

	mide = partition( mide->coord[axis], axis, start, end );
	int rank{int(mide-start)+1};

	if ( rank == (position + 1))
		return mide;
	else if ( rank > (position + 1))
		return select(start, --mide, position, axis);
	else
		return select(++mide, end, position - rank, axis);
}

bool check_sorting( const int axis, kpoint* low, kpoint* high){
	if ( low == high )
		return true;
	kpoint* tmp{low};
	++low;
	while ( low <= high - 1){
		if ( tmp -> coord[axis] > low -> coord[axis]){
			std::cout << "THE PROBLEM IS WITH\t" 
					  << tmp -> coord[axis] << "\t" 
					  << low -> coord[axis] << std::endl;
			return false;
		}
		++low;
		++tmp;
	}
	return true;
}

void sorting(const int axis, kpoint* low, kpoint* high){
	if ( low >= high )
		return ;
	kpoint* median{(int(high - low))/2 + low};
	median = partition(median->coord[axis], axis, low, high);
	if ( low != median )
		sorting( axis , low , median - 1);
	if ( high != median )
		sorting( axis , median + 1, high);
	return;
}

kpoint* partition( const float_t median, const int axis, kpoint* low, kpoint* high){
	// partion the data and return the pointer to the median
	// that's necessary for the continuation of the main algorithms
	kpoint* tmplow; kpoint* tmpmid; kpoint* tmphigh;
	tmplow = low;
	tmpmid = low;
	tmphigh = high;

	while ( tmpmid <= tmphigh ){
		if ( tmpmid -> coord[axis] < median ){
			tmpmid -> swap(*tmplow);
			++tmplow;
			++tmpmid;
		}
		else if ( tmpmid -> coord[axis] > median ){
			tmpmid -> swap(*tmphigh);
			--tmphigh;
		}
		else
			++tmpmid;
	}
	tmplow = low;
	while (tmplow -> coord[axis] != median)
		++tmplow;

	return tmplow;
}


// function usefull to control the quality of the kdtrees

int height(knode* node){
	if ( (node -> left == nullptr ) && ( node -> right == nullptr ) )
		return 0;
	if ( (node -> left == nullptr ) && ( node -> right != nullptr ) )
		return 1 + height(node -> right );
	if ( (node -> left != nullptr ) && ( node -> right == nullptr ) )
		return 1 + height(node -> left );
	return 1 + std::max( height(node -> left), height(node -> right) );
}

int count_subtree(knode* node, int height){
	if ( height == 0 )
		return 1;
	--height;
	if ( (node -> left == nullptr ) && ( node -> right == nullptr ) )
		return 1;
	if ( ((node -> left) == nullptr ) && ( (node -> right) != nullptr ) )
		return 1 + count_subtree(node -> right, height);
	if ( ((node -> left) != nullptr ) && ( (node -> right) == nullptr ) )
		return 1 + count_subtree(node -> left , height);
	return 1 +  count_subtree(node -> left, height) + count_subtree(node -> right, height) ;
}
bool check_kdtree(knode* node){
	if ( (node -> left == nullptr ) && ( node -> right == nullptr ) )
		return true;
	if ( (node -> left == nullptr ) && ( node -> right != nullptr ) )
		return (node->split.coord[node->axis] <= node->right->split.coord[node->axis]) && check_kdtree(node->right);

	if ( (node -> left != nullptr ) && ( node -> right == nullptr ) )
		return (node->split.coord[node->axis] >= node->left->split.coord[node->axis]) && check_kdtree(node->left);

	return (node->split.coord[node->axis] <= node->right->split.coord[node->axis]) && 
	(node->split.coord[node->axis] >= node->left->split.coord[node->axis]) && 
	check_kdtree(node->right) && check_kdtree(node->left);
}

knode* find_kpoint_fast(knode* node, kpoint target){
	if ( node->split == target )
		return node;
	if ( node->left != nullptr && (target[node->axis] <= node->split[node->axis]) ){
		knode* tmp {find_kpoint_fast( node->left, target)};
		if ( tmp != nullptr)
			return tmp;
	}
	if ( node->right != nullptr && (target[node->axis] >= node->split[node->axis])){
		knode* tmp {find_kpoint_fast( node->right, target)};
		if ( tmp != nullptr)
			return tmp;
	}
	return nullptr;
}
