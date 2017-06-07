/*******************************************************************
 * PointBucket.hpp
 * A Basic PointBucket implementation
 * Used for storing/binning Points
 *******************************************************************/

#ifndef Point_BUCKET_HPP
#define Point_BUCKET_HPP

#include <vector>
#include <iostream>
#include <iomanip>
#include "pct/Point.hpp"
#include "BlockScheme.hpp"
#include <mpi.h>
#include "Grid.hpp"
using namespace std;

struct PointBucket
{
	int id;
	int limit;
	int count;
	char baseDir[1024];
	struct Point* pts;
	PointBucket();
	PointBucket(int _id, int _limit, int _count, char* _baseDir);
	~PointBucket();
	bool insert(const struct Point *pt);
	void clear();
	//int Point remove(int idx);
	int maxDim();
	void sortByDim(int dim);
	vector<int> getBlockList(BlockScheme *bScheme);
	int blockSort(BlockScheme *bScheme);
	int makeDir();
	int grid(Grid* destGrid);
	//void split(struct PointBucket* bucket2);
	int dump(const char* bucketFilePath);
	int read(const char* bucketFilePath);
	int send(int rank, MPI_Comm* comm, MPI_Request* req, MPI_Status* status);
	void print();
	// read();
	// write();
};


#endif
