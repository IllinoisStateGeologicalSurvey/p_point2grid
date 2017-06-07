#include <stdio.h>
#include <iomanip>
#include <vector>
#include <stdlib.h>
#include "Point.hpp"
#include "PointBucket.hpp"
#include "util.hpp"
#include "CRS.hpp"

using namespace std;

/*double randomVal(int _range, int _min) {
	int newMin = _min;
	bool isNeg = false;
	double randDbl = 0.0;
	if (_min < 0) {
		newMin = abs(_min);
		isNeg = true;
	} 
	if (isNeg) {
		randDbl = (double) (rand() % _range - newMin);
	} else {
		randDbl = (double) (rand() % _range + newMin);
	}
	return randDbl;
};
*/
int main(int argc, char* argv[])
{
	struct point test;
	int limit = 1000000;
	CRS *crs = new CRS();
	int bucket_id = 0;//double x = (double) (rand() % 360 - 180);
	char tmp_dir[1024] = "/gpfs_scratch/ncasler/tmp/";
	//double y = (double) (rand() % 180 - 90);
	struct PointBucket bucket(1, limit, 0, &tmp_dir[0]);
	struct PointBucket bucket2(2,limit, 0, &tmp_dir[0]);
	struct PointBucket bucket3(3, limit,0, &tmp_dir[0]);
	for (int i = 0; i < limit; i++) {
		//x = randomVal(360, -180);
		//y = randomVal(180, -90);
		test.randomize();
		crs->transform(&test);
		//test.print();
		//test.update(x, y);
		bucket.insert(&test);
	}
	//bucket.print();
	int dim = bucket.maxDim();
	printf("Max Dim: %i\n", dim);
	//bucket.sortByDim(dim);
	char bucketPath[31] = "test1";
	//bucket.print();
	
	//bucket.print();
	// Try to add one more point to see if we can segfault
	if (bucket.insert(&test)) {
		fprintf(stderr, "Insert successful\n");
	} else {
		fprintf(stderr, "Failed to insert, emptying\n");
	
	}
	int i = 0;
	bucket.dump(bucketPath);
	//bucket.split(&bucket2);
	//bucket.print();
	printf("Bucket1 has %i points\n", bucket.count);
	//bucket2.print();
	//for (i = limit -1; i > -1; --i) {
	//	bucket.remove(i);
	//}
	bucket3.read(bucketPath);
	printf("New Bucket count: %i\n", bucket3.count);
	//printf("%i\n",bucket.count);
	//bucket.print();
	delete crs;
	return 0;
}
