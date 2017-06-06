#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdexcept>
#include <mpi.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <float.h>
#include <algorithm>
#include "Point.hpp"
#include "Grid.hpp"
#include "GridPoint.hpp"
#include "DTypes.h"

#include "gdal.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"
using namespace std;

Grid::Grid() {
	cols = 0;
	rows = 0;
	datatype = DT_Unknown;
	data = NULL;
	origin = new point();
}
/* Note Currently origin is assuming lower-left corner for index function,
  This is not standard for raster data nor is it recommended. Need to change for
  interoperability
*/
Grid::Grid(const struct point *_origin, int _cols, int _rows, DType _datatype, double _resX, double _resY) 
{
	origin = new point(*_origin);
	cols = _cols;
	rows = _rows;
	datatype = _datatype;
	data = NULL;
	resX = _resX;
	resY = _resY;
}

Grid::~Grid() {
	cols = 0;
	rows = 0;
	datatype = DT_Unknown;
	dealloc();
	resX = 0.0;
	resY = 0.0;
	delete(origin);
}
/** Allocate the memory necessary for the given grid **/
void * Grid::alloc() {
	int d_size = GetDataTypeSize(datatype);
	if (!d_size > 0) {
		fprintf(stderr, "GridError: Invalid datatype prevented allocation\n");
		return NULL;
	}
	if (data != NULL) {
		fprintf(stderr, "GridError: Memory already allocated for grid data\n");
		return NULL;
	}
	fprintf(stdout, "Allocating %ix%i %s cells", cols, rows, GetDataTypeName(datatype));
	data = (GridPoint*)malloc(sizeof(GridPoint) * cols * rows);
	//data = malloc(sizeof(cols * rows * d_size));
	return data;
}

void Grid::dealloc() {
	if (data != NULL) {
		free(data);
		data = NULL;
	}
}

// Create a geotiff from the raster
int Grid::write(char* outPath, char* projSRS) {
	GDALDatasetH hDataset;
	
	GDALAllRegister();
	const char *pszFormat = "GTiff";
	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	char **papszMetadata;
	char **papszOptions = NULL;
	OGRSpatialReference hSRS;
	char *pszSRS_WKT = NULL;
	GDALRasterBandH hBand;

	if ( hDriver == NULL )
		exit (1 );

	papszMetadata = GDALGetMetadata( hDriver, NULL );
	if ( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ) )
		printf( "Driver %s supports Create() method.\n", pszFormat );
	if ( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATECOPY, FALSE ) )
		printf( "Driver %s supports CreateCopy() method.\n", pszFormat );
	hDataset = GDALCreate( hDriver, outPath, cols, rows, 1, GDT_Float32, papszOptions);
	

	double adfGeoTransform[6] = { origin->x, resX, 0, origin->y, 0, resY * -1};
	// Geotransform for the raster
	GDALSetGeoTransform(hDataset, adfGeoTransform);
	hSRS.importFromProj4(projSRS);
	// Create WKT definition of projection
	OSRExportToWkt(hSRS, &pszSRS_WKT);
	OSRDestroySpatialReference(hSRS);
	GDALSetProjection(hDataset, pszSRS_WKT);
	CPLFree( pszSRS_WKT);

	hBand = GDALGetRasterBand(hDataset, 1);
	GDALRasterIO(hBand, GF_Write, 0, 0, cols, rows, data, cols, rows, GDT_Float32, 0,0);
	/*Close the dataset*/
	GDALClose( hDataset);
	return 1;
}




size_t Grid::getSize() {
	size_t gridSize = (size_t)cols * rows * GetDataTypeSize(datatype);
	return gridSize;
}


double Grid::getMaxX() {
	double max = origin->x + (resX * cols);
	return max;
}

double Grid::getMinY() {
	double max = origin->y - (resY * rows);
	return max;
}

int Grid::getCol(double coord) {
	if (coord > getMaxX() || coord < origin->x) {
		return -1;
	}
	int idx = ceil((coord - origin->x) / resX);
	return idx;
}

int Grid::getRow(double coord) {
	if (coord < getMinY() || coord > origin->y) {
		return -1;
	}
	int idx = ceil((origin->y - coord) / resY);
	return idx;
}
/** Test if a point intersects with the grid **/
int Grid::within(const struct point *pt) {
	int flag = 0;
	if (pt->x >= origin->x && pt->x <= getMaxX()) {
		if (pt->y >= getMinY() && pt->y <= origin->y) {
			flag = 1;
		} else {
			cout << "Point out of Row range" << endl;
		} 
	} else {
			cout << "Point out of Col range" << endl;
	}
	return flag;
}
/* Assumes that the point is in raster coordinates */
int Grid::set(const struct point* pt) {
	int i, j;
	int count = cols * rows;
	i = (pt->y * rows) + pt->x;
	if (i < 0 || i > count) {
		return 0;
	}
	GridPoint* cell = &data[i];
	if (cell->empty) {
		cell->sum = pt->z;
		cell->count = 1;
		return 1;
	} else {
		float tmpSum = cell->sum + pt->z;
		cell->sum = tmpSum;
		cell->count++;
		return 1;
	}
	
}

int Grid::getCell(const struct point *pt, int* idx) {
	int i, j = 0;
	if (!within(pt))  {
		cout << "Point outside grid bounds\n";
		return 0;
	}
	j = getCol(pt->x);
	i = getRow(pt->y);
	idx[0] = j;
	idx[1] = i;
	return 1;
}
	
