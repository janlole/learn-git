#include <iostream>
#include <random>
#include <fstream>
#include <string>

#define LIMIT 1
#define POINTS 10

class Point
{
private:
	double xco,yco;
public:
	Point(double xcoord, double ycoord): xco{xcoord}, yco{ycoord} {}
	Point(): xco{0}, yco{0} {}

	double x() const {return xco;}
	double y() const {return yco;}
	void x(double d){xco = d;}
	void y(double d){yco = d;}

	~Point(){}
};

class Grid
{
private:
	Point* elem;
	std::size_t sz;
public:
	Grid(std::size_t s) : elem{new Point[s]}, sz{s} {} 
	Grid() : elem{new Point[0]}, sz{0} {} 

	Point* elements() const {return elem;}

	Point& operator[](std::size_t i) { return elem[i]; }

	~Grid() {delete[] elem; }
};

void savegrid(const std::size_t npoint, char const name_file[]);

Grid loadgrid(const std::size_t npoint, char const name_file[]);

int main(int argc, char const *argv[])
{
	char const name_file[]{"example.txt"};
	savegrid(POINTS, name_file);
	Grid gridnew(POINTS);
	gridnew = loadgrid(POINTS, name_file);

	return 0;
}

void savegrid(const std::size_t npoint, char const name_file[]){
	// random generator inizialitation
	std::uniform_real_distribution<double> unif(-LIMIT,LIMIT);
	std::default_random_engine re{static_cast<long unsigned int>(time(0))};
	//__________________________________________________________

	Grid newgrid(npoint);
	
	for (std::size_t i{0}; i < npoint; ++i){
		newgrid[i].x(unif(re));
		newgrid[i].y(unif(re));
		
	}
	std::ofstream myfile;
	myfile.open(name_file);

	for (std::size_t i{0}; i < npoint; ++i){
		myfile << newgrid[i].x() << "\t" << newgrid[i].y() << "\n";
		
	}
  	myfile.close();
}


Grid loadgrid(const std::size_t npoint, char const name_file[]){
	double a,b;
	Grid newgrid(npoint);
	std::size_t counter{0};

	std::ifstream myfile;
	myfile.open(name_file);

	while (myfile >> a){
		myfile >> b;
		newgrid[counter].x(a);
		newgrid[counter].y(b);
		++counter;
	}
	myfile.close();

	return newgrid;
}
