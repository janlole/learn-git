#include <iostream>
#include <vector>
#include <random>
#include <memory>
#include <cmath>

#define SEED 1234
#define NDIM 2
#define NUMPOINTS 10
#define LIMIT 1
#define MIN_DIST 1
// #define DOUBLE_PRECISION double

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
};

struct ksplit
{
	int axis;
	float_t median;
};

std::ostream& operator<<(std::ostream& os, const kpoint& p) {
	for (const auto& x : p)
		os << x << ",";
	return os;
}
std::ostream& operator<<(std::ostream& os, const knode& n) {
	os << "axis:\t" << n.axis << " "
	   << "kpoint:\t" << n.split << " ";
	if ( n.left != nullptr ){
		os << "\n\t-L-" << *n.left ;
		// return os;
	}
	if ( n.right != nullptr ){
		os << "\n\t-R-" << *n.right ;
		// return os;
	}
	return os;
}

ksplit axis_median(kpoint* low, const kpoint* high);

kpoint* closest_to_median(kpoint* low, const kpoint* high, const float_t median, const int axis);

kpoint* partition( const float_t median, const int axis, kpoint* low, kpoint* high);

int main(int argc, char const *argv[])
{
	
	std::cout << "sizeof(float_t):\t" << sizeof(float_t) << std::endl;
	std::cout << "sizeof(kpoint):\t" << sizeof(kpoint) << std::endl;
	std::cout << "sizeof(knode):\t" << sizeof(knode) << std::endl;
	std::cout << "sizeof(ksplit):\t" << sizeof(ksplit) << std::endl;
	std::cout << "typeid(float_t).name():\t" << typeid(float_t).name() << std::endl;
	std::cout << "typeid(double).name():\t" << typeid(double).name() << std::endl;


	std::vector<kpoint> Grid(NUMPOINTS);
	std::vector<knode> Nodes(NUMPOINTS);

	// Random assignment
	std::uniform_real_distribution<float_t> unif(-LIMIT,LIMIT);
	std::default_random_engine re{static_cast<long unsigned int>(SEED)};

	for ( auto& x : Grid ){
		for ( auto& y : x){
			y = unif(re);
		}
	}
	//_____________________________________________
	
	for ( auto& x : Grid ){
		std::cout << x << "\t" << &x << "\n";
	}

	Nodes[0].axis = NDIM + 1;
	Nodes[0].split.coord[0] = 0;
	Nodes[0].split.coord[1] = NUMPOINTS - 1 ;
	int check{50};
	int woorking_node{0};
	int last_node{0};
	while ( woorking_node >= 0 && check){
		if ( Nodes[woorking_node].axis == (NDIM + 1) ){
			kpoint* start{&( Grid[static_cast<int>(Nodes[woorking_node].split.coord[0])] )};
			kpoint* end{&( Grid[static_cast<int>(Nodes[woorking_node].split.coord[1])] )};


			ksplit newsplit{axis_median(start, end)};
			kpoint* newclosest{closest_to_median(start, end, newsplit.median, newsplit.axis)};
			newsplit.median = newclosest -> coord[newsplit.axis];
			// std::cout << "\n closest \n" << newclosest << "\t" << *newclosest << "\t" << "\n-----\n";
			kpoint* mide{partition(newsplit.median, newsplit.axis, start, end)};	
			// std::cout << "\n mide \n" << mide << "\n-----\n";

			if ( start != mide ){
				++last_node;
				Nodes[last_node].axis = NDIM + 1;
				Nodes[last_node].split.coord[0] = Nodes[woorking_node].split.coord[0];
				Nodes[last_node].split.coord[1] = mide - start - 1 + Nodes[woorking_node].split.coord[0];
				Nodes[woorking_node].left = &Nodes[last_node];
				// std::cout << "\nHHHH\n" << Nodes[last_node].split.coord[1] << "\nHHHH\n" << std::endl;
			}
			if ( end != mide ){
				++last_node;
				Nodes[last_node].axis = NDIM + 1;
				Nodes[last_node].split.coord[0] = mide - start + 1 + Nodes[woorking_node].split.coord[0];
				Nodes[last_node].split.coord[1] = Nodes[woorking_node].split.coord[1];
				Nodes[woorking_node].right = &Nodes[last_node];
				// std::cout << "\nHHHH\n" << Nodes[last_node].split.coord[0] << "\nHHHH\n" << std::endl;
			}

			Nodes[woorking_node].split = *mide;
			Nodes[woorking_node].axis = newsplit.axis;
			std::cout << Nodes[woorking_node]<< "\n" << std::endl;
			woorking_node = last_node;
		}
		else{
			--woorking_node;
			// std::cout << "\nworking_node:\t" << woorking_node << "\n";
		}
		--check;
	}

	std::cout << "\n\n############\n";
	
	for ( auto &x : Grid )
		std::cout << x << "\t" << &x <<"\n";
	std::cout << "\n\n############\n";
	std::cout << Nodes[0] << std::endl;
	return 0;
}


ksplit axis_median(kpoint* low, const kpoint* high){
	// return a structure containing the most extense axis and its median
	std::vector<float_t> tmp(NDIM*2);
	ksplit result;
	int axis{0};
	while ( low <= high ){
		for ( int ax{0}; ax < NDIM; ++ax){
			if ( low -> coord[ax] < tmp[ax*NDIM] )
				tmp[ax*NDIM] = low ->coord[ax];
			if ( low ->coord[ax] > tmp[ax*NDIM+1] )
				tmp[ax*NDIM+1] = low ->coord[ax];
		}
		++low;
	}
	for ( int ax{0}; ax < NDIM; ++ax){
			if ( (tmp[ax*NDIM + 1] - tmp[ax*NDIM]) > (tmp[axis*NDIM + 1] - tmp[axis*NDIM]))
				axis = ax;
	}
	result.median =  (tmp[axis*NDIM + 1] + tmp[axis*NDIM])/2;
	result.axis = axis;
	return result;
}


kpoint* closest_to_median(kpoint* low, const kpoint* high, const float_t median, const int axis){	
	// return the pointer to the closest kpoint to the median, given a dimension axis
	kpoint* closest{low};
	float_t dist{MIN_DIST};
	while ( low <= high ){
		if ( std::abs( low -> coord[axis] - median) < dist ){
			dist = std::abs( low -> coord[axis] - median);
			closest = low;
		}
		++low;
	}
	return closest;
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
	return -- mid;
}
/*

	ksplit newsplit{axis_median(&Grid[0], &Grid[NUMPOINTS]-1)};

	kpoint* newclosest{closest_to_median(&Grid[0], &Grid[NUMPOINTS] -1 , newsplit.median, newsplit.axis)};

	newsplit.median = newclosest -> coord[newsplit.axis];

	kpoint* mide{partition(newsplit.median, newsplit.axis, &Grid[0], &Grid[NUMPOINTS] - 1)};	
*/
/*
	std::cout << "\n\n############\n";
	
	for ( auto &x : Grid )
		std::cout << x << "\t" << &x <<"\n";
*/