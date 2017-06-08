#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <unistd.h>
#include "points2grid/lasfile.hpp"
#include "pct/Grid.hpp"
#include "pct/BBox.hpp"
#include "pct/DTypes.hpp"
#include "pct/Point.hpp"
#include "pct/BlockScheme.hpp"
#include "pct/FileCollection.hpp"
#include "pct/Pixel.hpp"
#include "pct/util.hpp"
#include <float.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char* argv[]) {
	/****************************************************/
	/*   Parameter declaration                          */
	/****************************************************/
	// Index params
	int i,j,k = 0;
	
	// Directory path params
	char input_path[1024] = "/projects/isgs/lidar/champaign/las";
	char scratch_path[1024] = "/gpfs_scratch/ncasler/data/tmp";
	char out_path[1024] = "";
	char tmp_path[1024] = "";

	// MPI Params
	int world_size, world_rank, mpi_err;
	MPI_Comm world_comm = MPI_COMM_WORLD;
	MPI_Info info = MPI_INFO_NULL;
	MPI_Status status;
	MPI_Request request;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(world_comm, &world_size);
	MPI_Comm_rank(world_comm, &world_rank);
	MPI_Errhandler_set(world_comm, MPI_ERRORS_RETURN);
	double starttime, endtime;
	// Set memory limit for grid allocations
	long buffer_lim = 8000000000;
	size_t max_size = buffer_lim / sizeof(Pixel);
	// File specific parameters
	int n_files, file_off, file_blk, file_end = 0;
	char ext[5] = ".las";

	// Metadata
	double g_mins[3] = {DBL_MAX,DBL_MAX,DBL_MAX}; // Global min coord
	double g_maxs[3] = {-DBL_MAX,-DBL_MAX,-DBL_MAX}; // Global max coord
	double l_mins[3] = {DBL_MAX,DBL_MAX,DBL_MAX}; // Local min coord
	double l_maxs[3] = {-DBL_MAX,-DBL_MAX,-DBL_MAX}; // Local max coord

	// File specific params
	FileCollection *g_files = NULL; // Global file list
	int g_n_files = 0; // Global file count
	int l_n_files = 0; // Local file count
	BBox* g_file_bbox = NULL; // Global array of file bboxes
	BBox* l_file_bbox = NULL; // Local array of file bboxes
	// Grid specific params
	struct Point* origin = new Point();
	DType datatype = DT_Float32;
	struct Grid* g_grid = NULL; // Global grid
	int g_cols = 0; // Global grid column count
	int g_rows = 0; // Global grid row count
	
	//Specify resolution > should be parameterized in prod
	double res = 5.0f;
	// Block scheme params
	struct BlockScheme* blk_scheme = NULL;
	int *l_blks = NULL; // Local block array
	int *g_blks = NULL; // Global block array
	int l_n_blks = 0; // Local block count
	int g_n_blks = 0; // Global block count
	int *blk_n_files = NULL; // Array holding block specific file counts
	// Las specific params
	las_file las;
	long long g_n_pts = 0; // Global point count
	long long l_n_pts = 0; // Local point count

	/**************************************************/
	/*       Begin first file scan                    */
	/**************************************************/
	starttime = MPI_Wtime();
	g_files = new FileCollection(&input_path[0], &ext[0]);
	// Check the number of available LAS files
	g_n_files = g_files->countFiles();
	g_file_bbox = new BBox[g_n_files];
	l_file_bbox = new BBox[l_n_files];
	// Create tif output dir
	sprintf(&tmp_path[0], "%s/blocks", scratch_path);
	printf("[%i] Using tmp dir: %s\n", tmp_path);
	struct stat st = {0};
	if (stat(tmp_path, &st) == -1) 
		mkdir(tmp_path, 0700);
	// Set file Block size
	file_blk = ceil((float)g_n_files /(float)world_size);
	file_off = file_blk * world_rank;
	
	
	
	if (file_off + file_blk > g_n_files)
		file_end = g_n_files - file_off;
	else
		file_end = file_off + file_blk;

	// Read subset of file paths from dir
	l_n_files = g_files->getMetadata(file_off, file_end);

	// Scan metadata from files
	for (i = file_off; i < file_end; i++) {
		las.open(g_files->fileList[i]);
		l_n_pts = l_n_pts + (long long) las.points_count();
		compareMin(l_mins, las.minimums());
		compareMax(l_maxs, las.maximums());
		double *tmp_min = las.minimums();
		double *tmp_max = las.maximums();
		j = i - file_off;
		l_file_bbox[j].updateMin(tmp_min[0], tmp_min[1], tmp_min[2]);
		l_file_bbox[j].updateMax(tmp_max[0], tmp_max[1], tmp_max[2]);
		las.close();
	}
	endtime = MPI_Wtime();
	printf("[%i] Metadata gathered in %f seconds\n", world_rank, 
			endtime - starttime);
	MPI_Barrier(world_comm);
	/*****************************************************/
	/*         Gather global min/max point count         */
	/*                 COMMUNICATIONS                    */
	/*****************************************************/
	MPI_Allreduce(&l_mins[0], &g_mins[0], 3, MPI_DOUBLE, MPI_MIN, world_comm);
	MPI_Allreduce(&l_maxs[0], &g_maxs[0], 3, MPI_DOUBLE, MPI_MAX, world_comm);
	MPI_Allreduce(&l_n_pts, &g_n_pts, 1, MPI_LONG_LONG, MPI_SUM, world_comm);
	// Gather the bounding box values for each file
	printf("[%i] Gather bbox values\n", world_rank);
	MPI_Allgather(&l_file_bbox, l_n_files*sizeof(BBox), MPI_BYTE, &g_file_bbox, l_n_files, l_n_files*sizeof(BBox), world_comm);
	double io_time = MPI_Wtime();
	printf("[%i] Communication finished in %f seconds\n", world_rank, io_time - endtime);

	// Create origin for global grid
	origin->update(g_mins[0], g_maxs[0], g_mins[2]);
	g_cols = (int) ceil((g_maxs[0] - g_mins[0]) / res);
	g_rows = (int) ceil((g_maxs[1] - g_mins[1]) / res);
	// Create global grid
	g_grid = new Grid(origin, g_cols, g_rows, datatype, res, res);
	
	// Create global block scheme
	blk_scheme = new BlockScheme(g_grid, max_size, datatype, world_size);
	l_n_blks = blk_scheme->getBlockCount(world_rank); // Local block count
	g_n_blks = blk_scheme->cols * blk_scheme->rows; // Global block count
	l_blks = blk_scheme->getBlocks(world_rank); // Local block id array
	printf("[%i] Block Total: %i,Local: %i, first: %i, last: %i \n", g_n_blks, l_n_blks, l_blks[0], l_blks[l_n_blks-1]);
	
	/***********************************************************/
	/*                 Get file list for block                 */
	/***********************************************************/


	/***********************************************************/
	/*                 Read blocks                             */
	/***********************************************************/

	/***********************************************************/
	/*                Write                                    */
	/***********************************************************/

	/***********************************************************/
	/*        Clean up                                         */
	/***********************************************************/
	g_files->clear();
	printf("[%i]Cleaning up\n", world_rank);
	delete(origin);
	delete(g_files);
	delete[] l_file_bbox;
	delete[] g_file_bbox;
	delete(g_grid);
	delete(blk_scheme);

	MPI_Finalize();
	return 0;
}
