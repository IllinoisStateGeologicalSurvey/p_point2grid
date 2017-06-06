/************************************************************
 * Grid.hpp
 * A basic grid implementation for uniform meshing
 * Used/for binning points for outcore map-reduce
 ***********************************************************/

#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include <iostream>
#include <iomanip>
#include "Point.hpp"
#include "DTypes.h"
#include "GridPoint.hpp"
using namespace std;

struct Grid
{
	int cols;
	int rows;
	struct point* origin;
	double resX;
	double resY;
	DType datatype; // Currently unused,need to parameterize this
	GridPoint* data;
	Grid();
	Grid(const struct point *origin, int cols, int rows, DType datatype, double resX, double resY);	
	~Grid();
	void * alloc(); // Allocate memory for grid pixels
	void dealloc(); // Free memory taken by grid's pixels
	size_t getSize(); // Get size of grid in bytes
	int getCol(double coord);
	int getRow(double coord);
	int write(char* outPath, char* projSRS);
	int within(const struct point *pt);
	double getMaxX();
	double getMinY();
	int set(const struct point *pt);
	int getCell(const struct point *pt, int* idx);
};


#endif
