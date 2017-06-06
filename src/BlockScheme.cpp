#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include "Point.hpp"
#include "Grid.hpp"
#include "DTypes.h"
#include "BlockScheme.hpp"

blockScheme::blockScheme() {
	cols = 0;
	rows = 0;
	b_cols = 0;
	b_rows = 0;
	datatype = DT_Unknown;
	proc_count = 0;
};

blockScheme::~blockScheme() {
	cols = 0;
	rows = 0;
	b_cols = 0;
	b_rows = 0;
	datatype = DT_Unknown;
	proc_count = 0;
};

void blockScheme::getBlockPosition(int rank, int* position) {
	position[0] = rank / b_cols;
	position[1] = rank % b_cols;
}
// Get the ID of a block based on its col/row location
int blockScheme::getBlockId(int x, int y) {
	int idx = y * cols + x;
	return idx;
}

// Return the blockID based on the pixel coordinates
int blockScheme::getBlock(long x, long y) {
	int i = 0;
	int b_x = floor((float)x / (float)b_cols);
	int b_y = floor((float)y / (float)b_rows);
	return getBlockId(b_x, b_y);
}
//TODO: Function to return a grid object from a block index
int blockScheme::getBlockGrid(int blockId, Grid* grid, int res) {
	int pos[2] = {0, 0};
	// Make sure the x and y indexes are assigned properly by this function
	getBlockPosition(blockId, &pos[0]);
	struct point* _origin = new point();
	_origin->x = pos[0] * b_cols;
	_origin->y = pos[1] * b_rows;
	_origin->z = 0.0;
	grid = new Grid(_origin, b_cols, b_rows, datatype, res, res);
	return 0;
}
// Get the process rank based on the blockId
int blockScheme::getBlockRank(int blockId) {
	return blockId % proc_count;
}

int blockScheme::getBlockCount(int rank) {
	int gridSize = cols * rows;
	int count = (int)gridSize / proc_count;
	int rem = fmod(gridSize, proc_count);
	if (rem > 0 && rank <= rem) {
		count++;
	}
	return count;
}

int* blockScheme::getBlocks(int rank) {
	int gridSize = cols * rows;
	int counter = rank;
	int i = 0;
	int n_blocks = getBlockCount(rank);
	int* blockIds = (int*)malloc(sizeof(int) * n_blocks);
	for (i = 0; i < n_blocks; i++) {	
		blockIds[i] = counter;
		counter = counter + proc_count;
	}
	return blockIds;
};

/** 
 * Alternate constructor: based on grid dimensions and block
 * pixel limit.
 * @param grid: The grid to subdivide
 * @param block_limit: number of cells per block
 */
blockScheme::blockScheme(Grid* grid, int block_limit, DType _datatype, int _proc_count) {
	int gridSize = grid->cols * grid->rows;
	//Start with block that is basically square from square root
	b_cols = (int)sqrt(block_limit);
	b_rows = block_limit / b_cols;
	
	int c_rem = fmod(grid->cols, b_cols);
	int r_rem = fmod(grid->rows, b_rows);
	cols = ceil((float)grid->cols / (float)b_cols);
	rows = ceil((float)grid->rows / (float)b_rows);
	printf("Leftovers: rows %i, cols %i\n", c_rem, r_rem);
	printf("Block matrix: dims(%ix%i) cols %i, rows %i\n", b_cols, b_rows, cols, rows);
	proc_count = _proc_count;
	datatype = _datatype;
};


