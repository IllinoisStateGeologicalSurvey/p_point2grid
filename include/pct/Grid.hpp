/************************************************************
 * Grid.hpp
 * A basic grid implementation for uniform meshing
 * Used/for binning Points for outcore map-reduce
 ***********************************************************/

#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include <iostream>
#include <iomanip>
#include "pct/Point.hpp"
#include "pct/DTypes.hpp"
#include "pct/Pixel.hpp"
using namespace std;

typedef struct Grid
{
	int cols;
	int rows;
	struct Point* origin;
	double resX;
	double resY;
	DType datatype; // Currently unused,need to parameterize this
	Pixel* data;
	Grid();
	Grid(const struct Point *origin, int cols, int rows, DType datatype, double resX, double resY);	
	~Grid();
	void * alloc(); // Allocate memory for grid pixels
	void dealloc(); // Free memory taken by grid's pixels
	int cellCount(); // Return the number of cells (cols * rows)
	size_t getSize(); // Get size of grid in bytes
	int getCol(double coord);
	int getRow(double coord);
	int write(char* outPath, int epsg);
	int within(const struct Point *pt);
	double getMaxX();
	double getMinY();
	int set(const struct Point *pt);
	Pixel* get(int col, int row);
	float* getSumArr(); // Return array representing the sum for each pixel
	
	int getCell(const struct Point *pt, int* idx);
} Grid;


#endif
