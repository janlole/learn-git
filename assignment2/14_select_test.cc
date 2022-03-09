#include <iostream>
#include <vector>
#include <random>
#include <memory>
#include <cmath>
#include <chrono>
#include <algorithm>

#if !defined(NDIM)
#define NDIM 2
#endif

#define LIMIT 1

#if !defined(NUMPOINTS)
#define NUMPOINTS 4095
#endif

#if !defined(DOUBLE_PRECISION_KPOINT)
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
	kpoint& operator=(const kpoint& x){
		auto tmp = x;
		(*this) = std::move(tmp);
		return *this;
	}

	//_____________________________________________

	bool operator==( kpoint& l){
		for ( auto i{0} ; i < NDIM ; ++i )
			if ( l.coord[i] != this -> coord[i] )
				return false;
		return true;
	}

	bool operator!=( kpoint& l){
		return !( *this == l);
	}

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

typedef kpoint* (COMP)(kpoint*,kpoint*, knode*, int);

template<COMP choosen_algorithm>
knode* build_kdtree_recursive(kpoint* low, kpoint* high, knode* node, int depth);

kpoint* core_algorithm(kpoint* low, kpoint* high, knode* node, int axis);

// sort algorithm: naive implementation of quick sort
void sorting(const int axis, kpoint* low, kpoint* high);

kpoint* partition( const float_t median, const int axis, kpoint* low, kpoint* high);

kpoint* select(kpoint* start, kpoint* end, const int position, const int axis);

// function usefull to control the quality of the kdtrees
int main(int argc, char *argv[])
{
	
	std::cout << "NDIM,tree_dimension,power,axis,time_select,time_kdtree,theoretic_time_kdtree,ratio_real_theoretic"  << std::endl;
	int sub_tree_size{0};
	std::vector<float> times(50);
	int start_pow{2};
	std::uniform_real_distribution<float_t> unif(-LIMIT,LIMIT);
	std::default_random_engine re{static_cast<long unsigned int>(time(0))};

	for (auto power{20}; power < 31; ++power ){
		sub_tree_size = pow(2,power+start_pow)-1;
		std::vector<kpoint> Grid(sub_tree_size);
		std::vector<knode> Nodes(sub_tree_size);
		for ( auto& x : Grid ){
			for ( auto& y : x){
				y = unif(re);
			}
		}
		for ( auto axis{0}; axis < NDIM; ++axis){
			float sum_per_cycle{0};
			int how_many_times{2};
			for (auto i{0}; i < how_many_times; ++i){
				std::random_shuffle ( Grid.begin(), Grid.end() );
				auto t0 = std::chrono::high_resolution_clock::now();
				kpoint* tmp = select(&Grid[0], &Grid[sub_tree_size-1], sub_tree_size/2, axis);
				auto t1 = std::chrono::high_resolution_clock::now();
				auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
				sum_per_cycle += static_cast<float>(elapsed.count());
			}
			
			times[power] += (sum_per_cycle/how_many_times)/NDIM;

			float theoretic_time{0};
			for (auto k{0}; k <= power ; ++k){
				theoretic_time += ((pow(2,power-k)) * times[k] );
			}

			std::random_shuffle ( Grid.begin(), Grid.end() );
			auto t0 = std::chrono::high_resolution_clock::now();
			if ( axis == NDIM - 1 ){
				knode* tmp = build_kdtree_recursive<core_algorithm>(&Grid[0], &Grid[sub_tree_size-1], &Nodes[0], axis);
			}
			auto t1 = std::chrono::high_resolution_clock::now();
			auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);

			float time{static_cast<float>(elapsed.count())};

			std::cout << NDIM << ","										// NDIM
					  << sub_tree_size << ","								// tree_dimension
					  << power + start_pow << ","							// power
					  << axis << ","										// axis
					  << (sum_per_cycle/how_many_times) / 1e6 << ","		// time_select
					  << time / 1e6 << ","									// time_kdtree
					  << theoretic_time / 1e6 << ","						// theoretic_time_kdtree
					  << time / theoretic_time //<< ","						// ratio_real_theoretic
					  << std::endl;
		
		}
	}
	return 0;
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

