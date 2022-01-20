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

ksplit axis_median(const std::vector<kpoint>& Grid);

std::size_t closest_to_median(const std::vector<kpoint>& Grid, const float_t median, const int axis);

void partition( const float_t median, const int axis, kpoint* low, kpoint* high);

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
	std::uniform_real_distribution<double> unif(-LIMIT,LIMIT);
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

	ksplit newsplit{axis_median(Grid)};

	std::size_t newclosest{closest_to_median(Grid, newsplit.median, newsplit.axis)};
	
	newsplit.median = Grid[newclosest].coord[newsplit.axis];

	partition(newsplit.median, newsplit.axis, &Grid[0], &Grid[NUMPOINTS] - 1);

	
	std::cout << "\n\n############\n";
	
	for ( auto &x : Grid )
		std::cout << x << "\t" << &x - &Grid[0] <<"\n";

	return 0;
}


ksplit axis_median(const std::vector<kpoint>& Grid){
	std::vector<float_t> tmp(NDIM*2);
	ksplit result;
	for ( auto x : Grid ){
		for ( int ax{0}; ax < NDIM; ++ax){
			if ( x.coord[ax] < tmp[ax*NDIM] )
				tmp[ax*NDIM] = x.coord[ax];
			if ( x.coord[ax] > tmp[ax*NDIM+1] )
				tmp[ax*NDIM+1] = x.coord[ax];
		}
	}
	for ( auto x : tmp )
		std::cout << x << "\n";
	int axis{0};
	for ( int ax{0}; ax < NDIM; ++ax){
			if ( (tmp[ax*NDIM + 1] - tmp[ax*NDIM]) > (tmp[axis*NDIM + 1] - tmp[axis*NDIM]))
				axis = ax;
	}
	result.median =  (tmp[axis*NDIM + 1] + tmp[axis*NDIM])/2;
	result.axis = axis;
	return result;
}


std::size_t closest_to_median(const std::vector<kpoint>& Grid, const float_t median, const int axis){	
	std::size_t closest{0};
	float_t dist{MIN_DIST};
	for ( auto& x : Grid ){
		if ( std::abs(x.coord[axis] - median) < dist ){
			dist = std::abs(x.coord[axis] - median);
			closest = &x - &Grid[0];
		}
	}
	return closest;
}

void partition( const float_t median, const int axis, kpoint* low, kpoint* high){
	kpoint* mid{low};
	while ( mid <= high ){
		std::cout << low << "\t" << mid << "\t" << high << "\n";
		std::cout << *low << "\t" << *mid << "\t" << *high << "\n";

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
}
/*
	// TEST : median
	float_t med{median(Grid, 1)};
	std::cout << med << "\n";

	for ( auto x : Grid ){
		std::cout << x << "\n";
	}
*/

/*
	\\ TEST : closest_to_median
	std::size_t close;
	close = closest_to_median(Grid, 0, 0);

	for ( auto x : Grid ){
		std::cout << x << "\n";
	}
	std::cout << "closest to zero\t" << Grid[close] << "\t" << close << std::endl;
*/
/*
	std::cout << "Grid address: " << &Grid << std::endl;
	std::cout << "sizeof(std::vector):\t" << sizeof(std::vector<float_t>) << std::endl;
	std::cout << "sizeof(std::unique_ptr):\t" << sizeof(std::unique_ptr<float_t>) << std::endl;
*/	

/*
	std::cout << "\n\n\n\n" << std::endl;
	std::cout << "###################" << std::endl;
	std::cout << "STD::VECTOR\n";
	std::cout << "###################\n" << std::endl;
	std::vector<float> vettore1{2,4,6};
	std::vector<float> vettore2{1,3,5};
	std::vector<std::vector<float>> vettorevettore{vettore1, vettore2};
	std::vector<std::vector<float>*> vettorevettore1{&vettore1, &vettore2};
	
	std::cout << "&vettore1:\t" << &vettore1 << "\n"
			  << "&vettore1[0]:\t" << &vettore1[0] << "\n"
			  << "vettore1[0]:\t" << vettore1[0] << "\n";
	
	std::cout << "###################\n" << std::endl;

	std::cout << "&vettore2:\t" << &vettore2 << "\n"
			  << "&vettore2[0]:\t" << &vettore2[0] << "\n"
			  << "vettore2[0]:\t" << vettore2[0] << "\n";

	std::cout << "###################\n" << std::endl;
	
	std::cout << "&vettorevettore:\t" << &vettorevettore << "\n"
			  << "&vettorevettore[0]:\t" << &vettorevettore[0] << "\n"
			  << "vettorevettore[0]:\t" << vettorevettore[0][0] << "\n"
			  << "&vettorevettore[1]:\t" << &vettorevettore[1] << "\n"
			  << "vettorevettore[1]:\t" << vettorevettore[1][0] << "\n";

	std::cout << "###################\n" << std::endl;
	
	std::cout << "&vettorevettore1:\t" << &vettorevettore1 << "\n"
			  << "&vettorevettore1[0]:\t" << &vettorevettore1[0] << "\n"
			  << "vettorevettore1[0]:\t" << vettorevettore1[0][0] << "\n"
			  << "&vettorevettore1[1]:\t" << &vettorevettore1[1] << "\n"
			  << "vettorevettore1[1]:\t" << vettorevettore1[1][0] << "\n";

*/
/*	for ( auto x : Grid ){
		std::cout << x << "\t" << &Grid[i] << "\t" << typeid(Grid[i]).name() 
				  << "\t" << &(Grid[i].coord[0])
				  // << "\t" << sizeof(Grid[i].coord[0])
				  << "\t" << &(Grid[i].coord[1]) << std::endl;
				  // << "\t" << sizeof(Grid[i].coord[1]) <<std::endl;
		++i;
		
	}

	kpoint tmp0{Grid[0]};
	Grid[0] = Grid[1];
	Grid[1] = tmp0;

	i =0 ;
	std::cout << "Grid address: " << &Grid << std::endl;
	for ( auto x : Grid ){
		std::cout << x << "\t" << &Grid[i] << "\t" << typeid(Grid[i]).name() 
				  << "\t" << &(Grid[i].coord[0])
				  << "\t" << &(Grid[i].coord[1]) 
				  << std::endl;
		++i;
	}
	
*/
	// kpoint(kpoint&& x) : coord{std::move(x.coord)} {
	// }
	// kpoint& operator=(kpoint&& x){
	// 	coord = std::move(x.coord);
	// 	return *this;
	// }