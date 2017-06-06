#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include "Point.hpp"
#include "PointBucket.hpp"
#include "Grid.hpp"
#include "GridPoint.hpp"
#include <math.h>
#include <float.h>
#include <algorithm>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>

#include <mpi.h>


using namespace std;
PointBucket::PointBucket() {
	id = 0;
	limit = 10000;
	count = 0;
	pts = (struct point *)malloc(sizeof(struct point) * limit);
	//pts = new vector<point>();
};

PointBucket::PointBucket(int _id, int _limit, int _count, char* _baseDir) {
	id = _id;
	limit = _limit;
	count = _count;
	pts = (struct point *)malloc(sizeof(struct point) * limit);
	strcpy(&baseDir[0], _baseDir);
	//pts = new vector<point>();
};

PointBucket::~PointBucket() {
	id = 0;
	limit = 0;
	count = 0;
	delete pts;
};


bool PointBucket::insert(const struct point *pt) {
	if (count < limit) {
		pts[count].update(pt->x, pt->y, pt->z);
		//pts[count] = *pt;
		//printf("Point Added: (%lf, %lf %lf)\n", pts[count].x, pts[count].y, pts[count].z);

		count++;
		return true;

	} else {
		return false;
	}
};

void PointBucket::clear() {
	int i = 0;
	for (i = count-1; i >= 0; --i){
		//pts[i] = NULL;
		pts[i].update(0.0,0.0,0.0);
		count--;
	}

}

/** TODO: Fix this functionality, it really seems to work better 
with vectors
**/
/*const struct point PointBucket::remove(int idx) {
	if (idx >= 0 && idx < count) {
	struct point pt2 = pts[idx];
	pts[idx].clear();
	count--;
	return pt2;
	}
}
*/
int PointBucket::maxDim() {
	double minX = DBL_MAX;
	double maxX = DBL_MIN;
	double minY = DBL_MAX;
	double maxY = DBL_MIN;
	double minZ = DBL_MAX;
	double maxZ = DBL_MIN;
	double delts[3] = {0.0, 0.0, 0.0};
	int i = 0;
	struct point pt;
	for (i = 0; i < limit; i++) {
		pt = pts[i];
		if (pt.x < minX) {
			minX = pt.x;
		} else if (pt.x > maxX) {
			maxX = pt.x;
		}
		if (pt.y < minY) {
			minY = pt.y;
		} else if (pt.y > maxY) {
			maxY = pt.y;
		}
		if (pt.z < minZ) {
			minZ = pt.z;
		} else if (pt.z > maxZ) {
			maxZ = pt.z;
		}
	}
	delts[0] = maxX - minX;
	delts[1] = maxY - minY;
	delts[2] = maxZ - minZ;
	int dim = 0;
	double maxDim = 0.0;
	for (i = 0; i < 3; i++) {
		printf("Dim %i, size: %f\n", i, delts[i]);
		if (delts[i] > maxDim) {
			dim = i;
			maxDim = delts[i];
		}
	}
	return dim;
}
		
// This function will split the bucket along a given dimension at specified key, points, lower points will
// remain in this bucket, upper will be transferred to second bucket
//int PointBucket::split(struct PointBucket *bucket2, int dim) {
	
/**void PointBucket::sortByDim(int dim) {
	switch(dim) {
		case 0: std::sort(pts->begin(), pts->end(), compareX);
				break;
		case 1: std::sort(pts->begin(), pts->end(), compareY);
				break;
		case 2: std::sort(pts->begin(), pts->end(), compareZ);
				break;
	}
}
**/	
// Dumb splitting operation, doesnt care about dimensions or sorting, those will be 
// placed in wrapper
/*void PointBucket::split(struct PointBucket* bucket2) {
	int n = count / 2;
	std::vector<point>::reverse_iterator rit;
	int i = count;
	struct point pt;
	for (i = count-1; i >= n; --i) {
		pt = remove(i);
		bucket2->insert(&pt);
		count--;
	}
}
*/
int PointBucket::makeDir() {
	char dirPath[2056];
	sprintf(&dirPath[0], "%s/%i", baseDir, id);
	struct stat st = {0};
	if (stat(dirPath, &st) == -1) {
		mkdir(dirPath, 0700);
	}
	return 0;
}


int PointBucket::dump(const char* bucketFileName) {
	FILE *fp;
	char bucketPath[2056];
	sprintf(&bucketPath[0], "%s/%i/%s.bkt", baseDir, id,bucketFileName);
	if (access(bucketPath, F_OK) != -1) {
		//File exists, remove it
		remove(bucketPath);
	}
	fp = fopen(bucketPath, "wb");
	if (fp == NULL) {
		fprintf(stderr, "Cannot open file at %s\n", bucketPath);
		perror("WriteError ");
		return errno;
	}
	fwrite(&id, sizeof(int), 1, fp);
	printf("Writing %i points...\n", count);
	fwrite(&count, sizeof(int), 1, fp);
	int i = 0;
	fwrite(&pts[0], sizeof(struct point), count, fp);
	int last = count - 1;
	int mid = last / 2;
	printf("Point %i, has values(%lf, %lf, %lf)\n", 0, pts[0].x, pts[0].y, pts[0].z);
	//printf("Point %i, has values(%x, %x, %x)\n", 1, pts[1].x, pts[1].y, pts[1].z);
	//printf("Point %i, has values(%x, %x, %x)\n", mid, pts[mid].x, pts[mid].y, pts[mid].z);
	printf("Point %i, has values(%lf, %lf, %lf)\n", last, pts[last].x, pts[last].y, pts[last].z);
	//for ( i = 0; i < count; i++) {
	//	fwrite(&pts->at(i), sizeof(struct point), 1, fp);
	//}
	fclose(fp);
	return 1;
	
};

int PointBucket::read(const char* bucketFileName) {
	FILE *fp;
	char bucketPath[2056];
	sprintf(&bucketPath[0], "%s/%i/%s.bkt", baseDir, id, bucketFileName);
	fp = fopen(bucketPath, "rb");
	if (fp == NULL) {
		fprintf(stderr, "Cannot open file at %s\n", bucketPath);
		perror("Error ");
		return errno;
	}
	int _count = 0;
	int _id = 0;
	fread(&_id, sizeof(int), 1, fp);
	fread(&_count, sizeof(int), 1, fp);
	printf("Bucket[%i] Reading %i points...\n", _id, _count);
	count = _count;
	id = _id;
	int i = 0;
	//struct point pt;
	//struct point *points = (struct point*)malloc(sizeof(struct point)*count);
	
	fread(&pts[0], sizeof(struct point), count, fp);
	int last = count - 1;
	int mid = last/2;
	//printf("Point %i, has values(%lf, %lf, %lf)\n", 0, pts[0].x, pts[0].y, pts[0].z);
	//printf("Point %i, has values(%x, %x, %x)\n", 1, pts[1].x, pts[1].y, pts[1].z);
	//printf("Point %i, has values(%x, %x, %x)\n", mid, pts[mid].x, pts[mid].y, pts[mid].z);
	
	//printf("Point %i, has values (%lf, %lf, %lf)\n", last, pts[last].x, pts[last].y, pts[last].z);
	//for (i = 0; i < count; i++) {
	//	fread(&pt, sizeof(struct point), 1, fp);	
	//	pts->push_back(pt);
		//}
	fclose(fp);
	//free(points);
	return 0;
}

// Return the number of buckets covered by this bucked
// This assumes the coordinates xy refer to grid coordinates
// This could also be replaced by a comparator operation to sort points based
// on block
std::vector<int> PointBucket::getBlockList(blockScheme* bScheme) {
	int cntr = 0;
	int i,j = 0;
	int block = 0;
	//int blockArr[count] = (int*)malloc(sizeof(int) * count);
	std::vector<int> blocks;
	for (i = 0; i < count; i++) {
		block = bScheme->getBlock(pts[i].x, pts[i].y);
		//blockArr[i] = block;
		for (j = 0; j < cntr+1; j++) {
			// Check this first to prevent out-of-range errors
			if (j == cntr) {
				// No match found
				blocks.push_back(block);
				cntr++;
				break;
			} 

			if (block == blocks.at(j)) {
				break;
			}
		}
	}
	return blocks;
};

/* Method to apply the points to grid */
/* Note assumes points are in raster coords */
int PointBucket::grid(Grid* destGrid) {
	int i;
	int count;
	for (i = 0; i < count; i++) {
		if (destGrid->set(&pts[i])) {
			count++;
		}
	}
	return count;
};



int PointBucket::send(int rank, MPI_Comm* comm, MPI_Request* req, MPI_Status *status) {
	size_t headerSize = sizeof(int);
	int err = 0;
	size_t bufSize = (sizeof(point) * count)+ headerSize;
	unsigned char buffer[bufSize];
	memcpy(&buffer[0], &count, headerSize);
	memcpy(&buffer[headerSize], &pts[0], sizeof(point)*count);
	fprintf(stdout, "sending Buffer to %i with %i pts.\n", rank, count);
	err = MPI_Isend(&buffer[0], bufSize, MPI_BYTE, rank, 1, *comm, req); 
	MPI_Wait(req, status);
	if (err) {
		fprintf(stderr, "ERROR: Failed to send buffer, error %i\n", err);
	}
	//clear();
	return err;
};
void PointBucket::print() {
	printf("Bucket size: %i\n", count);
	
	//for (vector<point>::iterator it = pts->begin(); it != pts->end(); ++it) {
	for (int i = 0; i < count; i++) {
		
		pts[i].print();
		cout << ", ";

		if (i % 10 == 0){
			cout << endl;
		}

	}
	cout << endl;
}
