#include <iostream>
#include <vector>
#include <random>
#include <memory>
#include <cmath>
#include <chrono>
#include <omp.h>


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

#if !defined(DOUBLE_PRECISION)
#define float_t float
#else
#define float_t double
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

struct ksplit
{
	int axis;
	kpoint* median_kpoint;
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
kpoint* build_one_knode(kpoint* first_kpoint, kpoint* last_kpoint, knode* node, int depth);

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

knode* wrong_node( knode* node ){
	if ( (node -> left != nullptr ) && ( node -> right == nullptr ) ){
		if ( (node->split.coord[node->axis] >= node->left->split.coord[node->axis]) )
			return wrong_node(node->left);
		else
			return node;
	}
	if ( (node -> left == nullptr ) && ( node -> right != nullptr ) ){
		if ( (node->split.coord[node->axis] <= node->right->split.coord[node->axis]) )
			return wrong_node(node->right);
		else
			return node;
	}
	if ( ((node -> left) != nullptr ) && ( (node -> right) != nullptr ) ){
		if ( (node->split.coord[node->axis] <= node->right->split.coord[node->axis])  && 
			(node->split.coord[node->axis] >= node->left->split.coord[node->axis]) ){
			if ( ! check_kdtree(node->left))
				return wrong_node(node-> left);
			else if ( ! check_kdtree(node->right))
				return wrong_node(node-> right);
			else nullptr;
		}
		else{
			return node;
		}
	}
	return nullptr;

}

knode* find_kpoint(knode* node, kpoint target){
	if ( node->split == target )
		return node;
	if ( node->left != nullptr ){
		knode* tmp {find_kpoint( node->left, target)};
		if ( tmp != nullptr)
			return tmp;
	}
	if ( node->right != nullptr ){
		knode* tmp {find_kpoint( node->right, target)};
		if ( tmp != nullptr)
			return tmp;
	}
	return nullptr;
}
knode* find_kpoint_fast(knode* node, kpoint target){
	if ( node->split == target )
		return node;
	if ( node->left != nullptr && (target[node->axis] <= node->split[node->axis]) ){
		knode* tmp {find_kpoint( node->left, target)};
		if ( tmp != nullptr)
			return tmp;
	}
	if ( node->right != nullptr && (target[node->axis] >= node->split[node->axis])){
		knode* tmp {find_kpoint( node->right, target)};
		if ( tmp != nullptr)
			return tmp;
	}
	return nullptr;
}




int main(int argc, char const *argv[])
{
	std::cout << "sizeof(float_t):\t" << sizeof(float_t) << std::endl;
	std::cout << "NUMPOINTS:\t" << NUMPOINTS << std::endl;

	std::vector<kpoint> Grid(NUMPOINTS);
	std::vector<kpoint> Grid_copy(NUMPOINTS);
	std::vector<knode> Nodes(NUMPOINTS);
	std::vector<knode> Nodes_copy(NUMPOINTS);


	// Random assignment
	std::uniform_real_distribution<float_t> unif(-19,-1);
	std::default_random_engine re{static_cast<long unsigned int>(SEED)};

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

	//_____________________________________________

	auto t1 = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();
	omp_set_dynamic(true);
	t1 = std::chrono::high_resolution_clock::now();
	#pragma omp parallel  num_threads(4) 
	{
		#pragma omp single
		{
			int threads {omp_get_num_threads()};
			std::cout << "THREADS\t" << threads << std::endl;
			#pragma omp task 
			{
				build_one_knode<core_algorithm>(&Grid[0], &Grid[NUMPOINTS-1], &Nodes[0], 0);
			}
		}
	}
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "TIME OMP ALGORITHM-core_algorithm_sorting\t"
			  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
			  << "\t milliseconds" << std::endl;



	#if defined(CHECK_CORRECTNESS)
			t1 = std::chrono::high_resolution_clock::now();
			Nodes_copy[0].axis = NDIM ;
			Nodes_copy[0].split.coord[0] = 0;
			Nodes_copy[0].split.coord[1] = NUMPOINTS - 1 ;
			build_kdtree_recursive<core_algorithm>(&Grid_copy[0], &Grid_copy[NUMPOINTS-1] ,&Nodes_copy[0], 0);
			t2 = std::chrono::high_resolution_clock::now();
			std::cout << "\n\n############\n\n";
			std::cout << "TIME RECURSIVE ALGORITHM\t"
					  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
					  << "\t milliseconds" << std::endl;
			
			knode* research;
			int count_not_found{0};
			int num_point{(std::min(NUMPOINTS/2, int(pow(2,13) -1) ))};
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

	return 0;
}

template<COMP choosen_algorithm>
kpoint* build_one_knode(kpoint* first_kpoint, kpoint* last_kpoint, knode* node, int depth){
	knode* last_node{node}; 
	int axis;
	kpoint* mide; kpoint* tmp;
	if ( first_kpoint == last_kpoint){
		node -> split = *first_kpoint;			
	}
	else{
		axis = (depth%NDIM);
		mide = choosen_algorithm(first_kpoint,last_kpoint,node,axis);		
		++depth;		

		if ( first_kpoint != mide ){
			last_node = node + 1;
			node -> left = last_node;
			tmp = mide - 1;

			if ( depth < 3){
				#pragma omp task firstprivate(first_kpoint, tmp, last_node, depth)
				{
					build_one_knode<choosen_algorithm>(first_kpoint, tmp, last_node, depth);
				}
			}
			else{
				 #pragma omp task firstprivate(first_kpoint, tmp, last_node, depth)
				 {
				 	 build_kdtree_recursive<choosen_algorithm>(first_kpoint, tmp, last_node, depth);
				 }
			}
		}
		
		if ( last_kpoint != mide ){
			last_node = node + ( mide - first_kpoint ) + 1;
			node -> right = last_node;
			tmp = mide + 1;

			if ( depth < 3){
				#pragma omp task firstprivate(last_kpoint, tmp, last_node, depth)
				{
					build_one_knode<choosen_algorithm>(tmp, last_kpoint, last_node, depth);
				}
			}
			else{
				#pragma omp task firstprivate(last_kpoint, tmp, last_node, depth)
				{
					build_kdtree_recursive<choosen_algorithm>(tmp, last_kpoint, last_node, depth);
				}
			}
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
































/*	int zerosR{0};
	int zerosL{0};
	std::vector<int> bugged_knode_L(NUMPOINTS/2);
	std::vector<int> bugged_knode_R(NUMPOINTS/2);
	for ( auto& x : Nodes ){
		if ( (&Nodes[NUMPOINTS-1] - x.left < 0) ){
			bugged_knode_L[zerosL] = &x - &Nodes[0];
			++zerosL;
			std::cout << "L" << std::endl;
			std::cout << x ;
			std::cout << &x - &Nodes[0] << std::endl;
			std::cout << x.left - &Nodes[NUMPOINTS-1] << "\n" << std::endl;
		}
		if ( (&Nodes[NUMPOINTS-1] - x.right < 0)){
			bugged_knode_R[zerosR] = &x - &Nodes[0];
			++zerosR;
			std::cout << "R" << std::endl;
			std::cout << x ;
			std::cout << &x - &Nodes[0] << std::endl;
			std::cout << x.right - &Nodes[NUMPOINTS-1] << "\n" << std::endl;
		}
	}
	if ( zerosL | zerosR ){
		std::cout << "found L\t" << zerosL << "\t"
				  << "found R\t" << zerosR << std::endl;
	}

	int existence{0}; int check_existence{0};
	for ( auto y : Grid){
		for ( auto x : Nodes ){
			if (x.split == y ){
				++check_existence;
				break;
			}
		}
	}

	std::cout << "FOUNDED KPOINTS\t" << check_existence << std::endl;

	for (auto& x : Nodes){
		x.left = nullptr;
		x.right = nullptr;
	}
	for ( auto i{0}; i < NUMPOINTS; ++i){
		Grid[i] = Grid_copy[i];
	}

	t1 = std::chrono::high_resolution_clock::now();
	#pragma omp parallel num_threads(1) 
	{
		#pragma omp single
		{
			int threads {omp_get_num_threads()};
			std::cout << "THREADS\t" << threads << std::endl;
			#pragma omp task // shared(Grid, Nodes)
			{
				build_one_knode<core_algorithm>(&Grid[0], &Grid[NUMPOINTS-1], &Nodes[0], 0);
			}
		}
	}
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "TIME OMP ALGORITHM-core_algorithm_sorting\t"
			  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
			  << "\t milliseconds" << std::endl;
	std::cout << "\n\n LEFT \n\n";
	for ( auto i: bugged_knode_L ){
		if (i != 0)
			std::cout << Nodes[i] << std::endl;
	}
	std::cout << "\n\n RIGHT \n\n";
	for ( auto i: bugged_knode_R ){
		if (i != 0)
			std::cout << Nodes[i] << std::endl;
	}

	// OUTSIDE FUNCTION: first node must be set before the function call 
	// the idea is that it can be set and construct by different processes
	zerosL = 0; zerosR = 0;
	for ( auto& x : Nodes ){
		if ( (&Nodes[NUMPOINTS-1] - x.left < 0) ){
			bugged_knode_L[zerosL] = &x - &Nodes[0];
			++zerosL;
			std::cout << "L" << std::endl;
			std::cout << x ;
			std::cout << &x - &Nodes[0] << std::endl;
			std::cout << x.left - &Nodes[NUMPOINTS-1] << "\n" << std::endl;
		}
		if ( (&Nodes[NUMPOINTS-1] - x.right < 0)){
			bugged_knode_R[zerosR] = &x - &Nodes[0];
			++zerosR;
			std::cout << "R" << std::endl;
			std::cout << x ;
			std::cout << &x - &Nodes[0] << std::endl;
			std::cout << x.right - &Nodes[NUMPOINTS-1] << "\n" << std::endl;
		}
	}
	if ( zerosL | zerosR ){
		std::cout << "found L\t" << zerosL << "\t"
				  << "found R\t" << zerosR << std::endl;
	}
*/


















