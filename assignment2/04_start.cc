#include <iostream>
#include <vector>
#include <random>
#include <memory>
#include <cmath>
#include <chrono>



#define NDIM 2
#define LIMIT 10
#define MIN_DIST 1

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

	bool operator==(const kpoint& l){
		for ( auto i{0} ; i < NDIM ; ++i ){
			if ( l.coord[i] != this -> coord[i] ){
				return false;}
		}
		return true;
	}
	bool operator!=(const kpoint& l){
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
		os << *n.left << "\n";
	}
	if ( n.right != nullptr ){
		os << "-R-" << &n << "-> " << n.right << "\n";
		os << *n.right << "\n";
	}
	return os;
}

ksplit choosing_split(kpoint* low, const kpoint* high);

kpoint* partition( const float_t median, const int axis, kpoint* low, kpoint* high);

knode* build_kdtree_recursive(kpoint* low, kpoint* high, knode* node);

knode* build_kdtree_iterative(kpoint* low, kpoint* high, std::vector<knode>& Nodes);

kpoint* core_algorithm(kpoint* low, kpoint* high, knode* node);

kpoint* core_algorithm_sorting(kpoint* low, kpoint* high, knode* node);

kpoint* core_algorithm_coparison(kpoint* low, kpoint* high, knode* node);

void sorting(const int axis, kpoint* low, kpoint* high);

bool check_sorting( const int axis, kpoint* low, kpoint* high);

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
		return (node->split.coord[node->axis] < node->right->split.coord[node->axis]) && check_kdtree(node->right);

	if ( (node -> left != nullptr ) && ( node -> right == nullptr ) )
		return (node->split.coord[node->axis] > node->left->split.coord[node->axis]) && check_kdtree(node->left);

	return (node->split.coord[node->axis] < node->right->split.coord[node->axis]) && check_kdtree(node->right) && 
	(node->split.coord[node->axis] > node->left->split.coord[node->axis]) && check_kdtree(node->left);
}


int main(int argc, char const *argv[])
{
	std::cout << "sizeof(float_t):\t" << sizeof(float_t) << std::endl;
	std::cout << "NUMPOINTS:\t" << NUMPOINTS << std::endl;



	std::vector<kpoint> Grid(NUMPOINTS);
	// std::vector<kpoint> Grid_copy(NUMPOINTS);
	std::vector<knode> Nodes(NUMPOINTS);
	// std::vector<knode> Nodes_copy(NUMPOINTS);


	// Random assignment
	std::uniform_real_distribution<float_t> unif(-LIMIT,LIMIT);
	std::default_random_engine re{static_cast<long unsigned int>(SEED)};

	for ( auto& x : Grid ){
		for ( auto& y : x){
			y = unif(re);
			for (auto warm{0}; warm < 10; ++warm)
				unif(re);
		}
	}
	//_____________________________________________
	
	// grid copy
	// for ( auto i{0}; i < NUMPOINTS; ++i){
	// 	Grid_copy[i] = Grid[i];
	// }
	auto t1 = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();

	t1 = std::chrono::high_resolution_clock::now();
	knode* root{build_kdtree_recursive(&Grid[0],(&Grid[NUMPOINTS])-1, &Nodes[0])};
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "TIME RECURSIVE ALGORITHM\t"
			  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
			  << "\t milliseconds" << std::endl;


	// std::cout << "\n\n############\n";
	// for ( auto &x : Grid_copy )
	// 	std::cout << x << "\t" << &x <<"\n";

	std::cout << "\n\n############\n";
	for (auto i{0}; i < 10; ++i){
		std::cout << count_subtree(&Nodes[0],i) << "\t" << pow(2,i+1) - 1 << "\n";
	}
	std::cout << "\n\n############\n";

	t1 = std::chrono::high_resolution_clock::now();
	// bool mammamia(Nodes[0]  == Nodes_copy[0]);
	bool equili(check_kdtree(&Nodes[0]));
	t2 = std::chrono::high_resolution_clock::now();
	// std::cout << "COMPARISON\t" << mammamia << std::endl;
	std::cout << "EQUILIBRIUM\t" << equili << std::endl;
	std::cout << "TIME COMPARISON ALGORITHM\t"
			  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
			  << "\t milliseconds" << std::endl;


	t1 = std::chrono::high_resolution_clock::now();
	int h{height(&Nodes[0])};
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "TIME HEIGHT ALGORITHM\t"
			  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
			  << "\t milliseconds" << std::endl;

	// t1 = std::chrono::high_resolution_clock::now();
	// int h_copy{height(&Nodes_copy[0])};
	// t2 = std::chrono::high_resolution_clock::now();
	// std::cout << "TIME HEIGHT ALGORITHM\t"
	// 		  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
	// 		  << "\t milliseconds" << std::endl;

	double teoretic_height{log2(NUMPOINTS+1)-1};

	std::cout << "height:\t" << h << std::endl;
	std::cout << "teoretic height:\t" << teoretic_height << std::endl;
	


	// std::cout << Nodes[0] << std::endl;
	

	return 0;
}

knode* build_kdtree_recursive(kpoint* low, kpoint* high, knode* node){
	if ( low == high ){
		// if one point is available, return a leaf without any computation
		node -> split = *low;
		return node;		
	}
	static knode* node_working{node};

	kpoint* mide{CORE_ALGORITHM(low,high,node)};

	if ( low != mide ){
		++node_working;
		node -> left = build_kdtree_recursive(low, mide-1, node_working);
	}
	if ( high != mide ){
		++node_working;
		node -> right = build_kdtree_recursive(mide+1, high, node_working);
	}
	return node;
}

kpoint* core_algorithm_coparison(kpoint* low, kpoint* high, knode* node){
	kpoint* mide_wrong{core_algorithm(low, high,  node)};
	knode fake_node{*node};
	kpoint* mide{core_algorithm_sorting(low,  high, node)};
	std::cout << "AXIS WRONG\t" << fake_node.axis << "\t"
			  << "AXIS\t" << node -> axis << "\n";
	std::cout << (node -> split.coord[node->axis] - fake_node.split.coord[fake_node.axis])/LIMIT << "\n" ;
	return mide;
}

kpoint* core_algorithm(kpoint* low, kpoint* high, knode* node){
	ksplit newsplit{choosing_split(low,high)};
	kpoint* mide{partition(newsplit.median_kpoint -> coord[newsplit.axis], newsplit.axis, low, high)};

	node -> split = *mide;
	node -> axis = newsplit.axis;
	return mide;
}

kpoint* core_algorithm_sorting(kpoint* low, kpoint* high, knode* node){
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

	sorting(axis, low, high);
	kpoint* mide{(int(high - low))/2 + low};

	node -> split = *mide;
	node -> axis = axis;
	return mide;

}

ksplit choosing_split(kpoint* low, const kpoint* high){
	// return a structure containing the most extense axis and its median
	std::vector<float_t> tmp(NDIM*2);
	ksplit result;
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
	float_t median{ (tmp[axis*NDIM + 1] + tmp[axis*NDIM])/2};

	low_func = low;
	kpoint* closest{low};
	float_t dist{std::abs( low -> coord[axis] - median)};
	while ( low_func <= high ){
		if ( std::abs( low_func -> coord[axis] - median) < dist ){
			dist = std::abs( low_func -> coord[axis] - median);
			closest = low_func;
		}
		++low_func;
	}
	result.axis = axis;
	result.median_kpoint = closest;
	return result;
}

bool check_sorting( const int axis, kpoint* low, kpoint* high){
	if ( low == high )
		return true;
	kpoint* tmp{low};
	++low;
	while ( low <= high - 1){
		if ( tmp -> coord[axis] > low -> coord[axis])
			return false;
		++low;
		++tmp;
	}
	return true;
}

void sorting(const int axis, kpoint* low, kpoint* high){
	if ( low == high )
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
	// that's necessary for the continuation of the main algorithm
	kpoint* mid{low};
	while ( mid <= high ){
		if ( mid -> coord[axis] < median ){
			mid -> swap(*low);
			++low;
			++mid;
		}
		else if ( mid -> coord[axis] > median ){
			mid -> swap(*high);
			--high;
		}
		else{
			++mid;
		}
	}
	return --mid;
}
/*
	t1 = std::chrono::high_resolution_clock::now();
	knode* base{build_kdtree_iterative(&Grid[0],(&Grid[NUMPOINTS])-1, Nodes)};
	t2 = std::chrono::high_resolution_clock::now();
	std::cout << "\n\n############\n";
	std::cout << "TIME ITERATIVE ALGORITHM\t"
			  << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
			  << "\t milliseconds" << std::endl;
*/
/*
knode* build_kdtree_iterative(kpoint* low, kpoint* high, std::vector<knode>& Nodes){
	Nodes[0].axis = NDIM + 1;
	Nodes[0].split.coord[0] = 0;
	Nodes[0].split.coord[1] = NUMPOINTS - 1 ;
	int woorking_node{0};
	int last_node{0};
	while ( woorking_node >= 0 ){
		if ( Nodes[woorking_node].axis == (NDIM + 1) ){
			kpoint* start{&( Grid[static_cast<int>(Nodes[woorking_node].split.coord[0])] )};
			kpoint* end{&( Grid[static_cast<int>(Nodes[woorking_node].split.coord[1])] )};
			if ( start == end){
				Nodes[woorking_node].split = *start;			
				Nodes[woorking_node].axis = 0;
			}
			else{
				kpoint* mide{core_algorithm(low,high,node)};	

				if ( start != mide ){
					++last_node;
					Nodes[last_node].axis = NDIM + 1;
					Nodes[last_node].split.coord[0] = Nodes[woorking_node].split.coord[0];
					Nodes[last_node].split.coord[1] = mide - start - 1 + Nodes[woorking_node].split.coord[0];
					Nodes[woorking_node].left = &Nodes[last_node];
				}
				if ( end != mide ){
					++last_node;
					Nodes[last_node].axis = NDIM + 1;
					Nodes[last_node].split.coord[0] = mide - start + 1 + Nodes[woorking_node].split.coord[0];
					Nodes[last_node].split.coord[1] = Nodes[woorking_node].split.coord[1];
					Nodes[woorking_node].right = &Nodes[last_node];
				}
				woorking_node = last_node;
			}
			
		}
		else{
			--woorking_node;
		}
	}
}
*/