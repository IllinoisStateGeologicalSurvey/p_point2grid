/*
*
COPYRIGHT AND LICENSE

Copyright (c) 2011 The Regents of the University of California.
All rights reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following
conditions are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided
with the distribution.

3. All advertising materials mentioning features or use of this
software must display the following acknowledgement: This product
includes software developed by the San Diego Supercomputer Center.

4. Neither the names of the Centers nor the names of the contributors
may be used to endorse or promote products derived from this
software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*
*
* Based on the notes by Prof. Ramon Arrowsmith(ramon.arrowsmith@asu.edu)
* Authors: Han S Kim (hskim@cs.ucsd.edu), Sriram Krishnan (sriram@sdsc.edu)
*
*/


#include <points2grid/config.h>
#include <points2grid/Interpolation.hpp>
#include <points2grid/Global.hpp>
#include <points2grid/GridPoint.hpp>
#include <points2grid/MpiInterp.hpp>

#include <time.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <mpi.h>


#ifdef HAVE_GDAL
#include "gdal_priv.h"
#include "ogr_spatialref.h"
#endif

MpiInterp::MpiInterp(double dist_x, double dist_y,
                           int size_x, int size_y,
                           double r_sqr,
                           double _min_x, double _max_x,
                           double _min_y, double _max_y,
                           int _window_size,
                           int _rank, int _process_count)
{
    rank = _rank;
    process_count = _process_count;
    comm_done = 0;

    GRID_DIST_X = dist_x;
    GRID_DIST_Y = dist_y;

    GRID_SIZE_X = size_x;
    GRID_SIZE_Y = size_y;

    radius_sqr = r_sqr;

    min_x = _min_x;
    max_x = _max_x;
    min_y = _min_y;
    max_y = _max_y;

    window_size = _window_size;

    cerr << "MpiInterp created successfully" << endl;
}

MpiInterp::~MpiInterp()
{
    free(interp);
}

int MpiInterp::init()
{
    int i, j;
    row_stride = GRID_SIZE_Y / (process_count - 1);
    if (rank != 0)
    {


        w_row_start_index = (rank - 1) * row_stride;
        w_row_end_index =   (rank * row_stride) - 1;
        int left_over = GRID_SIZE_Y % (process_count - 1);
        if (rank == process_count - 1)
        {
            w_row_end_index += left_over;
        }


        int row_count = w_row_end_index - w_row_start_index + 1;
        interp = (GridPoint**) malloc (sizeof(GridPoint *) * row_count);
        if (interp == NULL)
        {
            cerr << "MpiInterp::init() malloc error" << endl;
            return -1;
        }

        for (i = 0; i < row_count; i++)
        {
            interp[i] = (GridPoint *) malloc (sizeof(GridPoint) * GRID_SIZE_Y);
            if (interp[i] == NULL)
            {
                cerr << "MpiInterp::init() malloc error" << endl;
                return -1;
            }
        }

        for (i = 0; i < row_count; i++)
        {
            for (j = 0; j < GRID_SIZE_X; j++)
            {
                interp[i][j].Zmin = DBL_MAX;
                interp[i][j].Zmax = -DBL_MAX;
                interp[i][j].Zmean = 0;
                interp[i][j].count = 0;
                interp[i][j].Zidw = 0;
                interp[i][j].sum = 0;
                interp[i][j].Zstd = 0;
                interp[i][j].Zstd_tmp = 0;
                interp[i][j].empty = 0;
                interp[i][j].filled = 0;
            }
        }
    }
    cerr << "MpiInterp::init() done" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;

}

int MpiInterp::update(double data_x, double data_y, double data_z)
{
    double x;
    double y;

    int lower_grid_x;
    int lower_grid_y;

    //cerr << GRID_DIST_X << " " << GRID_DIST_Y;
    lower_grid_x = (int)floor((double)data_x/GRID_DIST_X);
    lower_grid_y = (int)floor((double)data_y/GRID_DIST_Y);

    if(lower_grid_x > GRID_SIZE_X || lower_grid_y > GRID_SIZE_Y)
    {
        cerr << "larger at (" << lower_grid_x << "," << lower_grid_y << "): ("<< data_x << ", " << data_y << ")" << endl;
        return 0;
    }

    //printf("lower_grid_x: %d, grid_y: %d, arrX: %.2f, arrY: %.2f\n", lower_grid_x, lower_grid_y, arrX[i], arrY[i]);
    x = (data_x - (lower_grid_x) * GRID_DIST_X);
    y = (data_y - (lower_grid_y) * GRID_DIST_Y);

    //cerr << "(" << data_x << " " << data_y << ")=(" << lower_grid_x << ", " << lower_grid_y << "): ";
    //printf("(%f %f) = (%d, %d): ", data_x, data_y, lower_grid_x, lower_grid_y);

    //if(lower_grid_y == 30 && data_y > GRID_DIST_Y * lower_grid_y)
    //printf("(%f %f) = (%d, %d)\n", data_x, data_y, lower_grid_x, lower_grid_y);

    update_first_quadrant(data_z, lower_grid_x+1, lower_grid_y+1, GRID_DIST_X - x, GRID_DIST_Y - y);
    update_second_quadrant(data_z, lower_grid_x, lower_grid_y+1, x, GRID_DIST_Y - y);
    update_third_quadrant(data_z, lower_grid_x, lower_grid_y, x, y);
    update_fourth_quadrant(data_z, lower_grid_x+1, lower_grid_y, GRID_DIST_X - x, y);

    //cerr << "test" << endl;
    return 0;
}

int MpiInterp::finish(char *outputName, int outputFormat, unsigned int outputType)
{
  return finish(outputName, outputFormat, outputType, 0, 0);
}

int MpiInterp::finish(char *outputName, int outputFormat, unsigned int outputType, double *adfGeoTransform, const char* wkt)
{
    int rc;
    int i,j;

    //struct tms tbuf;
    clock_t t0, t1;

    for(i = 0; i < GRID_SIZE_X; i++)
        for(j = 0; j < GRID_SIZE_Y; j++)
        {
            if(interp[i][j].Zmin == DBL_MAX) {
                //		interp[i][j].Zmin = NAN;
                interp[i][j].Zmin = 0;
            }

            if(interp[i][j].Zmax == -DBL_MAX) {
                //interp[i][j].Zmax = NAN;
                interp[i][j].Zmax = 0;
            }

            if(interp[i][j].count != 0) {
                interp[i][j].Zmean /= interp[i][j].count;
                interp[i][j].empty = 1;
            } else {
                //interp[i][j].Zmean = NAN;
                interp[i][j].Zmean = 0;
            }

	    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
            if(interp[i][j].count != 0) {
  	        interp[i][j].Zstd = interp[i][j].Zstd / (interp[i][j].count);
		interp[i][j].Zstd = sqrt(interp[i][j].Zstd);
	    } else {
                interp[i][j].Zstd = 0;
            }


            if(interp[i][j].sum != 0 && interp[i][j].sum != -1)
                interp[i][j].Zidw /= interp[i][j].sum;
            else if (interp[i][j].sum == -1) {
                // do nothing
            } else {
                //interp[i][j].Zidw = NAN;
                interp[i][j].Zidw = 0;
            }
        }

    // Sriram's edit: Fill zeros using the window size parameter
    if (window_size != 0) {
        int window_dist = window_size / 2;
        for (int i = 0; i < GRID_SIZE_X; i++)
            for (int j = 0; j < GRID_SIZE_Y; j++)
            {
                if (interp[i][j].empty == 0) {
                    double new_sum=0.0;
                    for (int p = i - window_dist; p <= i + window_dist; p++) {
                        for (int q = j - window_dist; q <= j + window_dist; q++) {
                            if ((p >= 0) && (p < GRID_SIZE_X) && (q >=0) && (q < GRID_SIZE_Y)) {
                                if ((p == i) && (q == j))
                                    continue;

                                if (interp[p][q].empty != 0) {
                                    double distance = max(abs(p-i), abs(q-j));
                                    interp[i][j].Zmean += interp[p][q].Zmean/(pow(distance,Interpolation::WEIGHTER));
                                    interp[i][j].Zidw += interp[p][q].Zidw/(pow(distance,Interpolation::WEIGHTER));
				    interp[i][j].Zstd += interp[p][q].Zstd/(pow(distance,Interpolation::WEIGHTER));
				    interp[i][j].Zstd_tmp += interp[p][q].Zstd_tmp/(pow(distance,Interpolation::WEIGHTER));
                                    interp[i][j].Zmin += interp[p][q].Zmin/(pow(distance,Interpolation::WEIGHTER));
                                    interp[i][j].Zmax += interp[p][q].Zmax/(pow(distance,Interpolation::WEIGHTER));

                                    new_sum += 1/(pow(distance,Interpolation::WEIGHTER));
                                }
                            }
                        }
                    }
                    if (new_sum > 0) {
                        interp[i][j].Zmean /= new_sum;
                        interp[i][j].Zidw /= new_sum;
                        interp[i][j].Zstd /= new_sum;
                        interp[i][j].Zstd_tmp /= new_sum;
                        interp[i][j].Zmin /= new_sum;
                        interp[i][j].Zmax /= new_sum;
                        interp[i][j].filled = 1;
                    }
                }
            }
    }

    t0 = clock();

    if((rc = outputFile(outputName, outputFormat, outputType, adfGeoTransform, wkt)) < 0)
    {
        cerr << "MpiInterp::finish outputFile error" << endl;
        return -1;
    }

    t1 = clock();

    cerr << "Output Execution time: " << (double)(t1 - t0)/ CLOCKS_PER_SEC << std::endl;


    return 0;
}

//////////////////////////////////////////////////////
// Private Methods
//////////////////////////////////////////////////////

void MpiInterp::update_first_quadrant(double data_z, int base_x, int base_y, double x, double y)
{
    int i;
    int j;
    //double temp;

    //printf("radius: %f ", radius_sqrt);

    for(i = base_x; i < GRID_SIZE_X; i++)
    {
        for(j = base_y; j < GRID_SIZE_Y; j++)
        {
            /*
              temp = (   ((i - base_x)*GRID_DIST + x) * ((i - base_x)*GRID_DIST + x) +
              ((j - base_y)*GRID_DIST + y) * ((j - base_y)*GRID_DIST + y)) ;
              printf("%f ", temp);
            */

            double distance = 	((i - base_x)*GRID_DIST_X + x) * ((i - base_x)*GRID_DIST_X + x) +
                                ((j - base_y)*GRID_DIST_Y + y) * ((j - base_y)*GRID_DIST_Y + y) ;

            if(distance <= radius_sqr)
            {
                //printf("(%d %d) ", i, j);
                //interp[i][j]++;

                // update GridPoint
                //printf("----------------------------i %i worker size x %i rank %i\n", i, GRID_STRIDE_X, rank);
                int target_rank = get_target_rank(i);
                //printf("----------------------------i %i worker size x %i rank %i\n", i, GRID_STRIDE_X, rank);
                updateGridPointSend(target_rank, i, j, data_z, sqrt(distance));

            } else if(j == base_y) {
                //printf("return ");
                return;
            } else {
                //printf("break ");
                break;
            }
        }
    }

    //cerr << "test2" << endl;
}


void MpiInterp::update_second_quadrant(double data_z, int base_x, int base_y, double x, double y)
{
    int i;
    int j;

    for(i = base_x; i >= 0; i--)
    {
        for(j = base_y; j < GRID_SIZE_Y; j++)
        {
            double distance = 	((base_x - i)*GRID_DIST_X + x) * ((base_x - i)*GRID_DIST_X + x) +
                                ((j - base_y)*GRID_DIST_Y + y) * ((j - base_y)*GRID_DIST_Y + y);

            if(distance <= radius_sqr)
            {
                //printf("(%d %d) ", i, j);
                //interp[i][j]++;
                int target_rank = get_target_rank(i);
                updateGridPointSend(target_rank, i, j, data_z, sqrt(distance));


            } else if(j == base_y) {
                return;
            } else {
                break;
            }
        }
    }

}


void MpiInterp::update_third_quadrant(double data_z, int base_x, int base_y, double x, double y)
{
    int i;
    int j;

    for(i = base_x; i >= 0; i--)
    {
        for(j = base_y; j >= 0; j--)
        {
            double distance = 	((base_x - i)*GRID_DIST_X + x) * ((base_x - i)*GRID_DIST_X + x) +
                                ((base_y - j)*GRID_DIST_Y + y) * ((base_y - j)*GRID_DIST_Y + y);

            if(distance <= radius_sqr)
            {
                //if(j == 30)
                //printf("(%d %d)\n", i, j);
                //interp[i][j]++;
                int target_rank = get_target_rank(i);
                updateGridPointSend(target_rank, i, j, data_z, sqrt(distance));
            } else if(j == base_y) {
                return;
            } else {
                break;
            }
        }
    }
}

void MpiInterp::update_fourth_quadrant(double data_z, int base_x, int base_y, double x, double y)
{
    int i, j;

    for(i = base_x; i < GRID_SIZE_X; i++)
    {
        for(j = base_y; j >= 0; j--)
        {
            double distance = 	((i - base_x)*GRID_DIST_X + x) * ((i - base_x)*GRID_DIST_X + x) +
                                ((base_y - j)*GRID_DIST_Y + y) * ((base_y - j)*GRID_DIST_Y + y);

            if(distance <= radius_sqr)
            {
                //printf("(%d %d) ", i, j);
                //interp[i][j]++;
                int target_rank = get_target_rank(i);
                updateGridPointSend(target_rank, i, j, data_z, sqrt(distance));
            } else if (j == base_y) {
                return ;
            } else {
                break;
            }
        }
    }
}

// this can only be called by reader processes, currently only rank 0 is a reader
// todo, modify reader_count to work with more than one reader

int MpiInterp::get_target_rank(int row_index){
    int reader_count = 1;
    int rtn = row_index/row_stride + reader_count;
    if (rtn < process_count)
    {
        return rtn;
    }
    else if (rtn == process_count  && row_index <= w_row_end_index)
    {
        return rtn - 1; // this handles the "left over" rows assigned to the last worker process
    }

    else
    {
        fprintf(stderr, "MpiInterp::get_target_rank: worker process overflow, logic error");
        return -1;
    }
}


void MpiInterp::updateGridPointSend(int target_rank, int x, int y, double data_z, double distance)
{

    comm_x = x;
    comm_y = y;
    comm_data_z = data_z;
    comm_distance = distance;

    //printf("-----------------------------rank %i, target_rank %i\n", rank, target_rank);
    //printf ("before comm_done %i, rank %i\n", comm_done, rank);
    MPI_Send(&comm_done, 1, MPI_INT, target_rank, 1, MPI_COMM_WORLD);
    //printf ("comm_done %i, rank %i\n", comm_done, rank);
    if (!comm_done)
    {
        MPI_Send (&comm_x, 1, MPI_INT, target_rank, 1, MPI_COMM_WORLD);
        MPI_Send (&comm_y, 1, MPI_INT, target_rank, 1, MPI_COMM_WORLD);
        MPI_Send (&comm_data_z, 1, MPI_DOUBLE, target_rank, 1, MPI_COMM_WORLD);
        MPI_Send (&comm_distance, 1, MPI_DOUBLE, target_rank, 1,
                  MPI_COMM_WORLD);
    }

}


void MpiInterp::updateGridPointRecv()
{
    MPI_Status status;

    while(!comm_done){
        MPI_Recv(&comm_done, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        //printf ("comm_done %i, rank %i\n", comm_done, rank);
        if(comm_done)
        {
            break;
        }
        else
        {
            MPI_Recv(&comm_x, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&comm_y, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&comm_data_z, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&comm_distance, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
            updateGridPoint(comm_x, comm_y, comm_data_z, comm_distance);
        }
    }

}


void MpiInterp::updateGridPoint(int x, int y, double data_z, double distance)
{
    //todo update to work with more than one reader
    int reader_count = 1;
    y -= (rank-reader_count)*row_stride;
    //printf ("rank %i, x %i, y %i, stride %i\n", rank, x, y, GRID_STRIDE_X);


    if(interp[x][y].Zmin > data_z)
        interp[x][y].Zmin = data_z;
    if(interp[x][y].Zmax < data_z)
        interp[x][y].Zmax = data_z;

    interp[x][y].Zmean += data_z;
    interp[x][y].count++;

    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
    double delta = data_z - interp[x][y].Zstd_tmp;
    interp[x][y].Zstd_tmp += delta/interp[x][y].count;
    interp[x][y].Zstd += delta * (data_z - interp[x][y].Zstd_tmp);

    double dist = pow(distance, Interpolation::WEIGHTER);

    if(interp[x][y].sum != -1) {
        if(dist != 0) {
            interp[x][y].Zidw += data_z/dist;
            interp[x][y].sum += 1/dist;
        } else {
            interp[x][y].Zidw = data_z;
            interp[x][y].sum = -1;
        }
    } else {
        // do nothing
    }
}

void MpiInterp::printArray()
{
    int i, j;

    for(i = 0; i < GRID_SIZE_X; i++)
    {
        for(j = 1; j < GRID_SIZE_Y; j++)
        {
            cerr << interp[i][j].Zmin << ", " << interp[i][j].Zmax << ", ";
            cerr << interp[i][j].Zmean << ", " << interp[i][j].Zidw << endl;
            //printf("%.20f ", interp[i][j].Zmax);
        }
        //printf("\n");
    }
    cerr << endl;
}

int MpiInterp::outputFile(char *outputName, int outputFormat, unsigned int outputType, double *adfGeoTransform, const char* wkt)
{
    int i,j,k;

    //FILE **arcFiles;
    MPI_File *arcFiles;
    char arcFileName[1024];
    MPI_Offset *arcFileOffsets;
    char** arcFileBuffers;

    //FILE **gridFiles;
    MPI_File *gridFiles;
    char gridFileName[1024];
    MPI_Offset *gridFileOffsets;
    char** gridFileBuffers;

    int buffer_size = 1000000;


    const char *ext[6] = {".min", ".max", ".mean", ".idw", ".den", ".std"};
    unsigned int type[6] = {OUTPUT_TYPE_MIN, OUTPUT_TYPE_MAX, OUTPUT_TYPE_MEAN, OUTPUT_TYPE_IDW, OUTPUT_TYPE_DEN, OUTPUT_TYPE_STD};
    int numTypes = 6;



    // open ArcGIS files
    if(outputFormat == OUTPUT_FORMAT_ARC_ASCII || outputFormat == OUTPUT_FORMAT_ALL)
    {
        if((arcFiles = (MPI_File *)malloc(sizeof(MPI_File) *  numTypes)) == NULL)
        {
            cerr << "Arc MPI_File malloc error: " << endl;
            return -1;
        }
        arcFileOffsets = (MPI_Offset *) = malloc(sizeof(MPI_Offset) * numTypes);
        arcFileBuffers = (char **) malloc(sizeof(char *) * numTypes);

        for(i = 0; i < numTypes; i++)
        {
            if(outputType & type[i])
            {
                strncpy(arcFileName, outputName, sizeof(arcFileName));
                strncat(arcFileName, ext[i], strlen(ext[i]));
                strncat(arcFileName, ".asc", strlen(".asc"));
                if(rank == 0)
                {
                    unlink(arcFileName);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_File_open(MPI_COMM_WORLD, arcFileName, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &(arcFiles[i]));
                arcFileBuffers[i] = (char *) malloc(sizeof(char) * buffer_size);

                //if((arcFiles[i] = fopen(arcFileName, "w+")) == NULL)
                //{
                //    cerr << "File open error: " << arcFileName << endl;
                //    return -1;
                //}
            }
            else
            {
                arcFiles[i] = NULL;
            }
        }
    }
    else
    {
        arcFiles = NULL;
    }

    // open Grid ASCII files
    if(outputFormat == OUTPUT_FORMAT_GRID_ASCII || outputFormat == OUTPUT_FORMAT_ALL)
    {
        if((gridFiles = (MPI_File *)malloc(sizeof(MPI_File) * numTypes)) == NULL)
        {
            cerr << "File array allocation error" << endl;
            return -1;
        }
        gridFileOffsets = (MPI_Offset *) = malloc(sizeof(MPI_Offset) * numTypes);
        gridFileBuffers = (char **) malloc(sizeof(char *) * numTypes);

        for(i = 0; i < numTypes; i++)
        {
            if(outputType & type[i])
            {
                strncpy(gridFileName, outputName, sizeof(arcFileName));
                strncat(gridFileName, ext[i], strlen(ext[i]));
                strncat(gridFileName, ".grid", strlen(".grid"));
                if (rank == 0)
                {
                    unlink (gridFileName);
                }
                MPI_Barrier (MPI_COMM_WORLD);
                MPI_File_open (MPI_COMM_WORLD, gridFileName,
                               MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL,
                               &(gridFiles[i]));

                gridFileBuffers[i] = (char *) malloc(sizeof(char) * buffer_size);
                //if((gridFiles[i] = fopen(gridFileName, "w+")) == NULL)
                //{
                //    cerr << "File open error: " << gridFileName << endl;
                //    return -1;
                //}

            }
            else
            {
                gridFiles[i] = NULL;
            }
        }
    }
    else
    {
        gridFiles = NULL;
    }



    if (rank == 0)
    {
        // print ArcGIS headers
        if (arcFiles != NULL)
        {
            for (i = 0; i < numTypes; i++)
            {
                if (arcFiles[i] != NULL)
                {
                    fprintf (arcFiles[i], "ncols %d\n", GRID_SIZE_X);
                    fprintf (arcFiles[i], "nrows %d\n", GRID_SIZE_Y);
                    fprintf (arcFiles[i], "xllcorner %f\n", min_x);
                    fprintf (arcFiles[i], "yllcorner %f\n", min_y);
                    fprintf (arcFiles[i], "cellsize %f\n", GRID_DIST_X);
                    fprintf (arcFiles[i], "NODATA_value -9999\n");
                }
                MPI_File_get_position(arcFiles[i], &(arcFileOffsets[i]);
            }
        }

        // print Grid headers
        if (gridFiles != NULL)
        {
            for (i = 0; i < numTypes; i++)
            {
                if (gridFiles[i] != NULL)
                {
                    fprintf (gridFiles[i], "north: %f\n", max_y);
                    fprintf (gridFiles[i], "south: %f\n", min_y);
                    fprintf (gridFiles[i], "east: %f\n", max_x);
                    fprintf (gridFiles[i], "west: %f\n", min_x);
                    fprintf (gridFiles[i], "rows: %d\n", GRID_SIZE_Y);
                    fprintf (gridFiles[i], "cols: %d\n", GRID_SIZE_X);
                }
                MPI_File_get_position(gridFiles[i], &(gridFileOffsets[i]);
            }
        }
    } // if(rank == 0)

    MPI_Barrier (MPI_COMM_WORLD);

    if (arcFiles != NULL)
                {
                    for (i = 0; i < numTypes; i++)
                    {
                        if (arcFiles[i] != NULL)
                        {
                            MPI_Bcast (&(arcFileOffsets[i]), 1, MPI_OFFSET, 0,
                            MPI_COMM_WORLD)
                        }
                    }
                }


            if (gridFiles != NULL)
            {
                for (i = 0; i < numTypes; i++)
                {
                    if (gridFiles[i] != NULL)
                    {
                        MPI_Bcast (&(gridFileOffsets[i]), 1, MPI_OFFSET, 0,
                        MPI_COMM_WORLD)
                    }
                }
            }








    for(i = GRID_SIZE_Y - 1; i >= 0; i--)
    {
        for(j = 0; j < GRID_SIZE_X; j++)
        {
            if(arcFiles != NULL)
            {
                // Zmin
                if(arcFiles[0] != NULL)
                {
                    if(interp[j][i].empty == 0 &&
                            interp[j][i].filled == 0)
                        fprintf(arcFiles[0], "-9999 ");
                    else
                        fprintf(arcFiles[0], "%f ", interp[j][i].Zmin);
                }

                // Zmax
                if(arcFiles[1] != NULL)
                {
                    if(interp[j][i].empty == 0 &&
                            interp[j][i].filled == 0)
                        fprintf(arcFiles[1], "-9999 ");
                    else
                        fprintf(arcFiles[1], "%f ", interp[j][i].Zmax);
                }

                // Zmean
                if(arcFiles[2] != NULL)
                {
                    if(interp[j][i].empty == 0 &&
                            interp[j][i].filled == 0)
                        fprintf(arcFiles[2], "-9999 ");
                    else
                        fprintf(arcFiles[2], "%f ", interp[j][i].Zmean);
                }

                // Zidw
                if(arcFiles[3] != NULL)
                {
                    if(interp[j][i].empty == 0 &&
                            interp[j][i].filled == 0)
                        fprintf(arcFiles[3], "-9999 ");
                    else
                        fprintf(arcFiles[3], "%f ", interp[j][i].Zidw);
                }

                // count
                if(arcFiles[4] != NULL)
                {
                    if(interp[j][i].empty == 0 &&
                            interp[j][i].filled == 0)
                        fprintf(arcFiles[4], "-9999 ");
                    else
                        fprintf(arcFiles[4], "%d ", interp[j][i].count);
                }

		// count
                if(arcFiles[5] != NULL)
                {
                    if(interp[j][i].empty == 0 &&
                            interp[j][i].filled == 0)
                        fprintf(arcFiles[5], "-9999 ");
                    else
                        fprintf(arcFiles[5], "%f ", interp[j][i].Zstd);
                }
	    }

            if(gridFiles != NULL)
            {
                // Zmin
                if(gridFiles[0] != NULL)
                {
                    if(interp[j][i].empty == 0 &&
                            interp[j][i].filled == 0)
                        fprintf(gridFiles[0], "-9999 ");
                    else
                        fprintf(gridFiles[0], "%f ", interp[j][i].Zmin);
                }

                // Zmax
                if(gridFiles[1] != NULL)
                {
                    if(interp[j][i].empty == 0 &&
                            interp[j][i].filled == 0)
                        fprintf(gridFiles[1], "-9999 ");
                    else
                        fprintf(gridFiles[1], "%f ", interp[j][i].Zmax);
                }

                // Zmean
                if(gridFiles[2] != NULL)
                {
                    if(interp[j][i].empty == 0 &&
                            interp[j][i].filled == 0)
                        fprintf(gridFiles[2], "-9999 ");
                    else
                        fprintf(gridFiles[2], "%f ", interp[j][i].Zmean);
                }

                // Zidw
                if(gridFiles[3] != NULL)
                {
                    if(interp[j][i].empty == 0 &&
                            interp[j][i].filled == 0)
                        fprintf(gridFiles[3], "-9999 ");
                    else
                        fprintf(gridFiles[3], "%f ", interp[j][i].Zidw);
                }

                // count
                if(gridFiles[4] != NULL)
                {
                    if(interp[j][i].empty == 0 &&
                            interp[j][i].filled == 0)
                        fprintf(gridFiles[4], "-9999 ");
                    else
                        fprintf(gridFiles[4], "%d ", interp[j][i].count);
		}

                // count
                if(gridFiles[5] != NULL)
                {
                    if(interp[j][i].empty == 0 &&
                            interp[j][i].filled == 0)
                        fprintf(gridFiles[5], "-9999 ");
                    else
                        fprintf(gridFiles[5], "%f ", interp[j][i].Zstd);
                }
            }
        }
        if(arcFiles != NULL)
            for(k = 0; k < numTypes; k++)
            {
                if(arcFiles[k] != NULL)
                    fprintf(arcFiles[k], "\n");
            }
        if(gridFiles != NULL)
            for(k = 0; k < numTypes; k++)
            {
                if(gridFiles[k] != NULL)
                    fprintf(gridFiles[k], "\n");
            }
    }

#ifdef HAVE_GDAL
    GDALDataset **gdalFiles;
    char gdalFileName[1024];

    // open GDAL GeoTIFF files
    if(outputFormat == OUTPUT_FORMAT_GDAL_GTIFF || outputFormat == OUTPUT_FORMAT_ALL)
    {
        GDALAllRegister();

        if((gdalFiles = (GDALDataset **)malloc(sizeof(GDALDataset *) *  numTypes)) == NULL)
        {
            cerr << "File array allocation error" << endl;
            return -1;
        }

        for(i = 0; i < numTypes; i++)
        {
            if(outputType & type[i])
            {
                strncpy(gdalFileName, outputName, sizeof(gdalFileName));
                strncat(gdalFileName, ext[i], strlen(ext[i]));
                strncat(gdalFileName, ".tif", strlen(".tif"));

                char **papszMetadata;
                const char *pszFormat = "GTIFF";
                GDALDriver* tpDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

                if (tpDriver)
                {
                    papszMetadata = tpDriver->GetMetadata();
                    if (CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATE, FALSE))
                    {
                        char **papszOptions = NULL;
                        gdalFiles[i] = tpDriver->Create(gdalFileName, GRID_SIZE_X, GRID_SIZE_Y, 1, GDT_Float32, papszOptions);
                        if (gdalFiles[i] == NULL)
                        {
                            cerr << "File open error: " << gdalFileName << endl;
                            return -1;
                        } else {
                            if (adfGeoTransform)
                                gdalFiles[i]->SetGeoTransform(adfGeoTransform);
                            if (wkt)
                                gdalFiles[i]->SetProjection(wkt);
                        }
                    }
                }
            } else {
                gdalFiles[i] = NULL;
            }
        }
    } else {
      gdalFiles = NULL;
    }

    if (gdalFiles != NULL)
    {
        for (i = 0; i < numTypes; i++)
        {
            if (gdalFiles[i] != NULL)
            {
                float *poRasterData = new float[GRID_SIZE_X*GRID_SIZE_Y];
                for (int j = 0; j < GRID_SIZE_X*GRID_SIZE_Y; j++)
                {
                    poRasterData[j] = 0;
                }

                for(j = GRID_SIZE_Y - 1; j >= 0; j--)
                {
                    for(k = 0; k < GRID_SIZE_X; k++)
                    {
                        int index = j * GRID_SIZE_X + k;

                        if(interp[k][j].empty == 0 &&
                                interp[k][j].filled == 0)
                        {
                            poRasterData[index] = -9999.f;
                        } else {
                            switch (i)
                            {
                                case 0:
                                    poRasterData[index] = interp[k][j].Zmin;
                                    break;

                                case 1:
                                    poRasterData[index] = interp[k][j].Zmax;
                                    break;

                                case 2:
                                    poRasterData[index] = interp[k][j].Zmean;
                                    break;

                                case 3:
                                    poRasterData[index] = interp[k][j].Zidw;
                                    break;

                                case 4:
                                    poRasterData[index] = interp[k][j].count;
                                    break;

                                case 5:
                                    poRasterData[index] = interp[k][j].Zstd;
                                    break;
                            }
                        }
                    }
                }
                GDALRasterBand *tBand = gdalFiles[i]->GetRasterBand(1);
                tBand->SetNoDataValue(-9999.f);

                if (GRID_SIZE_X > 0 && GRID_SIZE_Y > 0)
                    tBand->RasterIO(GF_Write, 0, 0, GRID_SIZE_X, GRID_SIZE_Y, poRasterData, GRID_SIZE_X, GRID_SIZE_Y, GDT_Float32, 0, 0);
                GDALClose((GDALDatasetH) gdalFiles[i]);
                delete [] poRasterData;
            }
        }
    }
#endif // HAVE_GDAL

    // close files
    if(gridFiles != NULL)
    {
        for(i = 0; i < numTypes; i++)
        {
            if(gridFiles[i] != NULL)
                fclose(gridFiles[i]);
        }
    }

    if(arcFiles != NULL)
    {
        for(i = 0; i < numTypes; i++)
        {
            if(arcFiles[i] != NULL)
                fclose(arcFiles[i]);
        }
    }

    return 0;
}

