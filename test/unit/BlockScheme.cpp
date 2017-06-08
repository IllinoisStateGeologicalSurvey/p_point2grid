#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
#include "Point.hpp"
#include "Grid.hpp"
#include "DTypes.h"
#include "BlockScheme.hpp"
#include <math.h>
int main(int argc, char* argv[]) {
	struct point origin;
	DType d_type = DT_Float32;
	int mem_limit = 8000000; // Memory limit for blocks
	origin.update(-142.0, 42.0, 1000);
	//Test with current NED times 10
	int g_cols = 63721200;
	int g_rows = 28081200;
	//int g_cols = 2040000;
	//int g_rows = 273443;
	int procs = 200;
	double res[2] = {1.0, 1.0};
	struct grid *grid1 = new grid(&origin, g_cols, g_rows, d_type, res[0], res[1]);
	int block_limit = mem_limit / GetDataTypeSize(d_type);
	blockScheme *bScheme = new blockScheme(grid1, block_limit, d_type, procs);
	int b_count = bScheme->cols * bScheme->rows;
	int i = 0;
	int b_rank;
	int bucket_count = 0;
	int position[2] = {0,0};
	for (i = 0; i < b_count; i++) {
		bScheme->getBlockPosition(i, &position[0]);
		b_rank = bScheme->getBlockRank(i);
		bucket_count = ceil((float)b_count / float(procs)); 
		printf("Block: %i has pos: {%i, %i} and rank %i, %i bpp\n", i, position[0], position[1], b_rank, bucket_count);
	}
	delete(bScheme);
	delete(grid1);

	return 0;
}


