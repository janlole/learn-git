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

#if defined (SORTING)
#define CORE_ALGORITHM core_algorithm_sorting
#else
#define CORE_ALGORITHM core_algorithm
#endif

#if !defined(RANDOM)
#define SEED 1234
#else
#define SEED time(0)
#endif

#if !defined(NUMPOINTS)
#define NUMPOINTS 32767
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
	// swap method
	void swap(kpoint& p){
		auto tmp = p;
		p = *this;
		*this = tmp;
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


std::ostream& operator<<(std::ostream& os, const kpoint& p) {
	for (const auto& x : p)
		os << x << ",";
	return os;
}
std::ostream& operator<<(std::ostream& os, const knode& n) {
	os << "address: " << &n << "\n"
	   << "axis:" << n.axis << "\n"
	   << "kpoint:" << n.split << "\n";
	if ( n.left != nullptr ){
		os << "-L-" << &n << "-> " << n.left << "\n";
		// os << *n.left << "\n";
	}
	if ( n.right != nullptr ){
		os << "-R-" << &n << "-> " << n.right << "\n";
		// os << *n.right << "\n";
	}
	return os;
}

kpoint* partition( const float_t median, const int axis, kpoint* low, kpoint* high);

kpoint* core_algorithm(kpoint* low, kpoint* high, knode* node, int axis);

kpoint* core_algorithm_sorting(kpoint* low, kpoint* high, knode* node, int axis);

typedef kpoint* (COMP)(kpoint*,kpoint*, knode*, int);

template<COMP choosen_algorithm>
knode* build_kdtree_recursive(kpoint* low, kpoint* high, knode* node, int depth);

template<COMP choosen_algorithm>
void build_kdtree_iterative(kpoint* first_kpoint, knode* node);

template<COMP choosen_algorithm>
kpoint* build_one_knode(kpoint* first_kpoint, kpoint* last_kpoint, knode* node, int axis);

void sorting(const int axis, kpoint* low, kpoint* high);

bool check_sorting( const int axis, kpoint* low, kpoint* high);

kpoint* select(kpoint* start, kpoint* end, const int position, const int axis);

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
	int blocklens[1] = {1};
	MPI_Aint indices[1];
	MPI_Datatype old_type[1];
	old_type[0] = MPI_float_t;
	MPI_Get_address( &dummy_kpoint.coord, &indices[0]);
	indices[0] = 0;
	MPI_Type_create_struct(1, blocklens, indices, old_type, &MPI_KPOINT );
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
	int n_step{3};
	int sub_tree_size;
	for ( auto step{0}; step < n_step; ++step){
		if (rank < pow(2,step+1) && rank >= pow(2,step)){
			sub_tree_size = std::ceil(NUMPOINTS/(pow(2,step)));
		}
	}
	std::vector<kpoint> Grid(sub_tree_size);
	std::vector<knode> Nodes(sub_tree_size);
	
	std::vector<kpoint> Grid_copy(NUMPOINTS);
	std::vector<knode> Nodes_copy(NUMPOINTS);

	auto t1 = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();

	// Random assignment and tree of control
	std::uniform_real_distribution<float_t> unif(-LIMIT,LIMIT);
	std::default_random_engine re{static_cast<long unsigned int>(SEED)};

	MPI_Status status;

	if ( rank == root ){
		for ( auto& x : Grid ){
			for ( auto& y : x){
				y = unif(re);
				for (auto warm{0}; warm < 10; ++warm)
					unif(re);
			}
		}
		// grid copy
		for ( auto i{0}; i < NUMPOINTS; ++i){
			Grid_copy[i] = Grid[i];
		}
		t1 = std::chrono::high_resolution_clock::now();
		build_kdtree_recursive<core_algorithm>(&Grid_copy[0],(&Grid_copy[NUMPOINTS])-1, &Nodes_copy[0], 1);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "TEST KDtree\n"
				  << "TIME RECURSIVE ALGORITHM-core_algorithm\t"
				  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
				  << "\t milliseconds" << std::endl;

		std::cout << "i am sending\t" << Nodes_copy[0] << std::endl;
		MPI_Send(&Nodes_copy[0], 4, MPI_KNODE, 1, 0, MPI_COMM_WORLD);
	}

	//_____________________________________________

	if (rank == 1){
		std::cout << "before recv\t" << Nodes[0] << std::endl;
		MPI_Recv(&Nodes[0], 4, MPI_KNODE, 0, 0, MPI_COMM_WORLD, &status);
		std::cout << "after recv\t" << Nodes[0] << std::endl;
	}

	kpoint* first_kpoint{&Grid[0]};
	kpoint* mide;
	kpoint* last_kpoint{&Grid[NUMPOINTS-1]}; 
	knode* working_node{&Nodes[0]};
	knode* relink_node; knode* relink_node_received; knode* tmp_knode_ptr_1; 
	int axis{0}, process_message, size_message;

	// step is the number of sending messages and 2^(step+1) is total number of process
	for (auto step{0}; step < n_step; ++step){
		if (rank < pow(2,step) ){
			// construct one node with saved information
			axis = (axis + 1)%NDIM; 
			std::cout << "I'm processor\t" << rank << "\t and my working axis is\t" << axis << std::endl;
			std::cout << "I'm processor\t" << rank << "\t and step is\t" << step << "\t" 
					  << ""
					  << std::endl;

			mide = build_one_knode<core_algorithm>(first_kpoint, last_kpoint, working_node, axis);
			// send instruction for build next node to the (rank + 2^step) core
			process_message = rank + pow(2,step);
			MPI_Send(mide+1, (last_kpoint - mide) + 1, MPI_KPOINT, process_message,0, MPI_COMM_WORLD );
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
		    MPI_Get_count(&status, MPI_INT, &size_message);
			MPI_Recv(&Grid[0], size_message , MPI_KPOINT, process_message, 0 , MPI_COMM_WORLD, &status);
			MPI_Recv(&axis, 1, MPI_INT, process_message, 0, MPI_COMM_WORLD, &status); // axis
			MPI_Recv(&relink_node, 1, MPI_AINT, process_message, 0, MPI_COMM_WORLD, &status); // SAVE THE POINTER SOMEWHERE TO RELINK THE SUBTREES
			// save information for the job( arguments of core algorithm function)
			first_kpoint = &Grid[0];
			last_kpoint = &Grid[0] + size_message - 1;
		}
	}

	build_kdtree_recursive<core_algorithm>(first_kpoint, last_kpoint, working_node, axis+1);
	std::cout<< "I'm processor\t" << rank << "\tand i complete my sub-tree\n";

	// step is the number of sending messages and 2^(step+1) is total number of process
	for (auto step{n_step-1}; step >= 0; --step){
		if (rank < pow(2,step+1) && rank >= pow(2,step)){
			process_message = rank - pow(2,step);
			// send subtree and the link address to (rank - 2^step) core
			MPI_Send(&Nodes[0], working_node - &Nodes[0] + 1 , MPI_KNODE, process_message, 0 , MPI_COMM_WORLD);
			MPI_Send(&relink_node, 1, MPI_AINT, process_message, 0, MPI_COMM_WORLD); 
			tmp_knode_ptr_1 = &Nodes[0];
			MPI_Send(&tmp_knode_ptr_1, 1, MPI_AINT, process_message, 0, MPI_COMM_WORLD); 
		}
		if (rank < pow(2,step) ){
			process_message = rank + pow(2,step);
			// send instruction for build next node to the (rank + 2^step) core
			MPI_Probe(process_message, 0, MPI_COMM_WORLD, &status);
		    MPI_Get_count(&status, MPI_INT, &size_message);
			MPI_Recv(1 + working_node, size_message, MPI_KNODE, process_message,0, MPI_COMM_WORLD , &status);
			MPI_Recv(&relink_node_received, 1, MPI_AINT, process_message, 0, MPI_COMM_WORLD, &status); 

			// link process
			relink_node_received -> right = 1 + working_node;
			// re-arrange addresses of the committed/received subtree
			
			MPI_Recv(&relink_node_received, 1, MPI_AINT, process_message, 0, MPI_COMM_WORLD, &status); 
			knode* tmp_knode_ptr{1+working_node};
			while ( tmp_knode_ptr < 1+working_node + size_message ){
				tmp_knode_ptr -> left = tmp_knode_ptr->left - relink_node_received + working_node + 1;
				tmp_knode_ptr -> right = tmp_knode_ptr->right - relink_node_received + working_node +1;
				++tmp_knode_ptr;
			}
			working_node = tmp_knode_ptr;
		}
	}







	// t1 = std::chrono::high_resolution_clock::now();
	// build_kdtree_recursive<core_algorithm_sorting>(&Grid_copy[0],(&Grid_copy[NUMPOINTS])-1, &Nodes_copy[0], 0);
	// t2 = std::chrono::high_resolution_clock::now();
	// std::cout << "TIME RECURSIVE ALGORITHM-core_algorithm_sorting\t"
	// 		  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
	// 		  << "\t milliseconds" << std::endl;

	// OUTSIDE FUNCTION: first node must be set before the function call 
	// the idea is that it can be set and construct by different processes

	// t1 = std::chrono::high_resolution_clock::now();
	// Nodes_copy[0].axis = NDIM ;
	// Nodes_copy[0].split.coord[0] = 0;
	// Nodes_copy[0].split.coord[1] = NUMPOINTS - 1 ;
	// build_kdtree_iterative<core_algorithm>(&Grid_copy[0], &Nodes_copy[0]);
	// t2 = std::chrono::high_resolution_clock::now();
	// std::cout << "\n\n############\n\n";
	// std::cout << "TIME ITERATIVE ALGORITHM\t"
	// 		  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
	// 		  << "\t milliseconds" << std::endl;



	#if defined(CHECK_CORRECTNESS)
			knode* research;
			int count_not_found{0};
			int num_point{(std::min(NUMPOINTS, int(pow(2,10) -1) ))};
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
			count_not_found = 0;
			for ( auto kpoint_target {0}; kpoint_target < num_point; ++kpoint_target ){
				research = find_kpoint_fast(&Nodes_copy[0], Grid_copy[kpoint_target]);
				if ( research == nullptr ){
					++count_not_found;
				}
				research = find_kpoint_fast(&Nodes_copy[0], Grid_copy[NUMPOINTS - 1 - kpoint_target]);
				if ( research == nullptr ){
					++count_not_found;
				}
			}
			
			std::cout << "\n\n############\n\n";
			std::cout << "sizeof(float_t):\t" << sizeof(float_t) << std::endl;
			std::cout << "NUMBER OF POINTS WHICH WILL BE CHECKED\t" << num_point*2 << "\n" <<std::endl;
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
				std::cout << "BALANCE TEST\t FAILED\n"
						  << "EQUILIBRIUM Nodes\t" << equili << "\tEQUILIBRIUM_copy\t" << equili_sort << std::endl;
			}
			else
				std::cout << "BALANCE TEST\t PASSED\n";

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

	#endif





	MPI_Finalize();
	return 0;
}
template<COMP choosen_algorithm>
kpoint* build_one_knode(kpoint* first_kpoint, kpoint* last_kpoint, knode* node, int axis){
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
		if ( last_kpoint != mide ){
			// if we do not allocate memory for the right node the "re-link" will be simplier
			// just send the address of the node that do not have the right child yet
			// then re-send that address to the parent core so it can set on that address the correct right-child-subtree

			// ++last_node;
			// last_node -> axis = NDIM ;			
			// node -> right = last_node;
		}
	}
	return mide;
}
template<COMP choosen_algorithm>
void build_kdtree_iterative(kpoint* first_kpoint, knode* node){
	knode* working_node{node};	knode* last_node{node};
	int axis{0}; kpoint* mide; 
	while ( (working_node - node) >= 0 ){
		if ( working_node->axis >= NDIM  ){
			kpoint* start{first_kpoint + static_cast<int>(working_node->split[0]) };
			kpoint* end{first_kpoint + static_cast<int>(working_node->split[1]) };

			if ( start == end){
				working_node -> split = *start;			
				working_node -> axis = 0;
			}
			else{
				axis = working_node->axis - NDIM ;
				mide = choosen_algorithm(start,end,working_node,axis);
				axis = (axis + 1)%NDIM + NDIM; 

				if ( start != mide ){
					++last_node;
					last_node -> axis = axis ;
					last_node -> split[0] = (start - first_kpoint);
					last_node -> split[1] = mide - first_kpoint - 1 ;
					
					working_node -> left = last_node;
				}
				if ( end != mide ){
					++last_node;
					last_node -> axis = axis ;
					last_node -> split[0] = mide - first_kpoint + 1;
					last_node -> split[1] = (end - first_kpoint);
					
					working_node -> right = last_node;
				}
				working_node = last_node;
			}
			
		}
		else{
			--working_node;
		}
		// std::cout << working_node << std::endl;

	}
}

template<COMP choosen_algorithm>
knode* build_kdtree_recursive(kpoint* low, kpoint* high, knode* node, int depth){
	if ( low == high ){
		// if one point is available, return a leaf without any computation
		node -> split = *low;
		return node;		
	}

	#if defined(WIDE_AXIS)
		// the following code allow the algorithm to search for the axis with the most
		// wide extension, instead of a symple 
		std::vector<float_t> tmp(NDIM*2);
		int axis{0};
		kpoint* low_func{low};
		while ( low_func <= high ){
			for ( int ax{0}; ax < NDIM; ++ax){
				if ( low_func -> coord[ax] < tmp[ax*NDIM] )
					tmp[ax*NDIM] = low_func ->coord[ax];
				if ( low_func ->coord[ax] > tmp[ax*NDIM+1] )
					tmp[ax*NDIM+1] = low_func ->coord[ax];
			}
			++low_func;
		}
		for ( int ax{0}; ax < NDIM; ++ax){
				if ( (tmp[ax*NDIM + 1] - tmp[ax*NDIM]) > (tmp[axis*NDIM + 1] - tmp[axis*NDIM]))
					axis = ax;
		}
	#else
		int axis{depth%NDIM};
	#endif
	
	static knode* node_working;
	node_working = node;
	kpoint* mide{choosen_algorithm(low,high,node,axis)};
	++depth;
	kpoint* tmp;
	if ( low != mide ){
		++node_working;
		tmp = mide-1;
		node -> left = build_kdtree_recursive<choosen_algorithm>(low, tmp, node_working, depth);
	}
	if ( high != mide ){
		++node_working;
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

kpoint* core_algorithm_sorting(kpoint* low, kpoint* high, knode* node, int axis){
	kpoint* clow{low};
	kpoint* chigh{high};
	sorting(axis, low, high);
	kpoint* mide{(int(chigh - clow))/2 + clow };

	node -> split = *mide;
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

/*
	#if defined(CHECK_CORRECTNESS)
			knode* research;
			int count_not_found{0};
			int num_point{(std::min(NUMPOINTS, int(pow(2,10) -1) ))};
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
			count_not_found = 0;
			for ( auto kpoint_target {0}; kpoint_target < num_point; ++kpoint_target ){
				research = find_kpoint_fast(&Nodes_copy[0], Grid_copy[kpoint_target]);
				if ( research == nullptr ){
					++count_not_found;
				}
				research = find_kpoint_fast(&Nodes_copy[0], Grid_copy[NUMPOINTS - 1 - kpoint_target]);
				if ( research == nullptr ){
					++count_not_found;
				}
			}
			
			std::cout << "\n\n############\n\n";
			std::cout << "sizeof(float_t):\t" << sizeof(float_t) << std::endl;
			std::cout << "NUMBER OF POINTS WHICH WILL BE CHECKED\t" << num_point*2 << "\n" <<std::endl;
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
				std::cout << "BALANCE TEST\t FAILED\n"
						  << "EQUILIBRIUM Nodes\t" << equili << "\tEQUILIBRIUM_copy\t" << equili_sort << std::endl;
			}
			else
				std::cout << "BALANCE TEST\t PASSED\n";

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
	#endif

*/
