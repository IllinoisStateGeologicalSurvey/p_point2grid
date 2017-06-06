/******************************************
 * BlockScheme.hpp
 * Basic functionality for block-cyclic 
 * decomposition of a raster which will not fit 
 * in memory. Should be useful for partitioning rasters 
 * across processors.
 ********************************************/

#ifndef BLOCK_SCHEME_HPP
#define BLOCK_SCHEME_HPP

#include <stdlib.h>
#include <iomanip>
#include "Point.hpp"
#include "Grid.hpp"
#include "DTypes.h"

typedef struct blockScheme
{
	int cols; // Number of blocks in row of grid
	int rows; // Number of blocks in column of grid
	int b_cols; // Number of columns within block
	int b_rows; // Number of rows within block
	int proc_count;
	DType datatype; // Type of data within block
	blockScheme();
	blockScheme(Grid* grid, int block_limit, DType datatype, int proc_count);
	void getBlockPosition(int rank, int *position); // Return the col,row index from  block ID
	int getBlockGrid(int blockId, Grid* grid, int res); // Return a grid object representing block
	int getBlockId(int x, int y); // Get block ID based on col,row indexing
	int getBlock(long x, long y); // Get the block ID based on pixel coordinates
	int getBlockRank(int blockId); // Return the rank specified block is assigned to
	int getBlockCount(int rank); // Return the number of blocks to be mapped to a particular processor
	int* getBlocks(int rank); // Return the indexes of the processors blocks
	~blockScheme();
	//getBlock(int p_x, int p_y, int *block_x, int *block_y);
} blockScheme;

#endif


