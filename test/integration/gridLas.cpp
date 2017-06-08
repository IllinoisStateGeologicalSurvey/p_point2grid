#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <unistd.h>
#include "points2grid/lasfile.hpp"
#include "pct/Grid.hpp"
#include "pct/BBox.hpp"
#include "pct/DTypes.hpp"
#include "pct/Point.hpp"
#include "float.h"
#include "math.h"
#include <time.h>

int main(int argc, char* argv[]) {
	int i,j,k = 0;
	char filePath[1024] = "/projects/isgs/lidar/champaign/las/10421216.las";
	double res = 3.0f;
	time_t begin,end;
	struct Point* min = new Point(DBL_MAX, DBL_MAX, DBL_MAX);
	struct Point* max = new Point(-DBL_MAX, -DBL_MAX, -DBL_MAX);
	struct Point* p = new Point();
	struct Grid* grid = NULL;
	las_file las;
	double *mins = NULL;
	double *maxs = NULL;
	long long n_pts = 0;
	int cols, rows = 0;
	begin = time(NULL);
	las.open(filePath);
	n_pts = (long long)las.points_count();
	mins = las.minimums();
	maxs = las.maximums();
	DType datatype = DT_Float32;
	min->update(mins[0],maxs[1],mins[2]);
	cols = (int)ceil((maxs[0] - mins[0]) / res);
	rows = (int)ceil((maxs[1] - mins[1]) / res);
	grid = new Grid(min, cols, rows, datatype, res, res);
	printf("Grid dims: [%i, %i]\n", cols, rows);
	grid->alloc();
	int col, row = 0;
	printf("Grid size is %zu GB\n", grid->getSize() / 1000000000);
	for (i =0; i < n_pts; i++) {
		// Check for last returns
		if (las.getReturnNumber(i) == las.getReturnCount(i)) {	
			p->update(las.getX(i), las.getY(i), las.getZ(i));
	//	printf("Reading point %i with coords (%f,%f,%f)\n", i, p->x, p->y, p->z);
			grid->set(p);
		}
	}
	las.close();
	char outPath[1024] = "/gpfs_scratch/ncasler/data/tmp/gridTest.tif";
	printf("Writing grid out to %s\n", &outPath[0]);
	grid->write(&outPath[0], 3443);
	grid->dealloc();
	end = time(NULL);
	printf("Finished in %f seconds\n", difftime(end, begin));
	delete(grid);
	delete(p);
	delete(min);
	delete(max);



	return 0;
}
