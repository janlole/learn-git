#include <iostream>
#include <memory>

#define DIMENSION 2

class Point {
	std::unique_ptr<double[]> elem;
public:
	Point() noexcept = default; // mark "noexcept" if not acquiring resources, makes code faster
  	Point(const double x, const double y) : elem{new double[DIMENSION]} { elem[0] = x; elem[1] = y;}

	double& operator[](const unsigned int i) noexcept { return elem[i]; }
	const double& operator[](const unsigned int i) const noexcept { return elem[i]; }

	// value of x and y
	double x() const {return elem[0]; }
	double y() const {return elem[1]; }
	//_____________________________________________

	// pointer to x and y
	double* xp() const {return &elem[0]; }
	double* yp() const {return &elem[1]; }
	//_____________________________________________

	// deconstructor
	~Point() noexcept { } // dtor should not rise exception
	//_____________________________________________

	// Move ctor and assignment
	Point(Point&& x) noexcept = default;
	Point& operator=(Point&& x) noexcept = default;
	//_____________________________________________
	
	// Copy ctor and assignment:  the ORIGINAL copy ctor is NOT exception safe
	Point(const Point& v) : elem{new double[DIMENSION]} {
		std::copy(v.begin(), v.end(), begin());
	}
	Point& operator=(const Point& x){ // if acquiring resources DO NOT mark "noexcept"
	  // clean resources with auxiliaty function
		auto tmp = x;
		(*this) = std::move(x);
		return *this;
	}
	//_____________________________________________

	// for range
	const double* begin() const { return &elem[0]; }
	double* begin() { return &elem[0]; }
	
	const double* end() const { return &elem[DIMENSION]; }
	double* end() { return &elem[DIMENSION]; }
	//_____________________________________________

};

/*
class Grid {
	std::size_t _size;
	std::unique_ptr<Point[]> elem;
public:
	Grid() noexcept = default; // mark "noexcept" if not acquiring resources, makes code faster
  	Grid(const unsigned int l) : _size{l}, elem{new double[l]} {}

	double& operator[](const unsigned int i) noexcept { return elem[i]; }
	const double& operator[](const unsigned int i) const noexcept { return elem[i]; }
	
	std::size_t size() const { return _size; }
	
	~Grid() noexcept { } // dtor should not rise exception
	

	// for range
	const double* begin() const { return &elem[0]; }
	double* begin() { return &elem[0]; }
	
	const double* end() const { return &elem[_size]; }
	double* end() { return &elem[_size]; }
	//_____________________________________________

};

*/


std::ostream& operator<<(std::ostream& os, const Point& p) {
	for (const auto& x : p)
		os << x << ",";
	os << std::endl;
	return os;
}


int main() {
	Point vuoto();
	Point pieno(3,8);

	std::cout << "pieno\t" << pieno << std::endl;
	std::cout << "pieno\t" << pieno.xp() << "\t" << pieno.yp() << std::endl;
	std::cout << pieno << std::endl;

	return 0;
}
