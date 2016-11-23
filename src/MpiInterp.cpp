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
#include <points2grid/debug.hpp>

#include <time.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <mpi.h>
#include <sstream>
#include <assert.h>

#ifdef HAVE_GDAL
//#include "gdal_priv.h"
//#include "ogr_spatialref.h"
#include <ogr_api.h>
#include <ogr_spatialref.h>
#include <gdal.h>
#include <gdal_priv.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#endif

MpiInterp::MpiInterp(double dist_x, double dist_y,
                           int size_x, int size_y,
                           double r_sqr,
                           double _min_x, double _max_x,
                           double _min_y, double _max_y,
                           int _window_size,
                           int _rank, int _process_count, int _reader_count, long _buffer_size, mpi_times *_timer = NULL)
{
    rank = _rank;
    process_count = _process_count;
    reader_count = _reader_count;
    buffer_size = _buffer_size;
    timer = _timer;
    is_reader = 0;
    is_writer = 0;
    readers = (int *) malloc(sizeof(int)*process_count);
    writers = (int *) malloc(sizeof(int)*process_count);
    read_done = (int *) malloc(sizeof(int)*process_count);

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

    bigtiff = 1;
    epsg_code = 0;

    const int nitems=5;
    int          blocklengths[5] = {1,1,1,1,1};
    MPI_Datatype types[5] = {MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint     offsets[5];

    offsets[0] = offsetof(grid_point_info, comm);
    offsets[1] = offsetof(grid_point_info, x);
    offsets[2] = offsetof(grid_point_info, y);
    offsets[3] = offsetof(grid_point_info, data_z);
    offsets[4] = offsetof(grid_point_info, distance);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_grid_point_info);
    MPI_Type_commit(&mpi_grid_point_info);

    const int nitems2 = 10;
    int blocklengths2[10] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    MPI_Datatype types2[10] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
            MPI_UNSIGNED, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
            MPI_INT, MPI_INT };
    MPI_Aint offsets2[10];

    offsets2[0] = offsetof(GridPoint, Zmin);
    offsets2[1] = offsetof(GridPoint, Zmax);
    offsets2[2] = offsetof(GridPoint, Zmean);
    offsets2[3] = offsetof(GridPoint, count);
    offsets2[4] = offsetof(GridPoint, Zidw);
    offsets2[5] = offsetof(GridPoint, Zstd);
    offsets2[6] = offsetof(GridPoint, Zstd_tmp);
    offsets2[7] = offsetof(GridPoint, sum);
    offsets2[8] = offsetof(GridPoint, empty);
    offsets2[9] = offsetof(GridPoint, filled);

    MPI_Type_create_struct (nitems2, blocklengths2, offsets2, types2,
                            &mpi_grid_point);
    MPI_Type_commit (&mpi_grid_point);



    //cerr << "MpiInterp created successfully" << endl;
}

MpiInterp::~MpiInterp()
{
    free(interp);
}

int MpiInterp::init()
{
    int i, j;
    writer_count = 0;
    w_row_start_index = w_row_end_index = 0;

    if (rank < reader_count)
    {
        is_reader = 1;

    }

    row_stride = GRID_SIZE_Y/(process_count-reader_count) + 1;
    // possibly reset row_stride to ensure window filling works correctly
    // this should be a very rare case, i.e., a small DEM with large number of processes...
    if(window_size && row_stride<(window_size/2))
    {
        row_stride = window_size/2 + 1;
    }

    dbg(5, "window_size %i, row_stride %i, rank %i", window_size, row_stride, rank);


    if(rank >= reader_count)
    {
        // idea here is that last writer will always have row_count < row_stride
        // all earlier rank writers will have row_count = row_stride -1
        w_row_start_index = (rank - reader_count) * row_stride;
        w_row_end_index =   ((rank + 1) - reader_count) * row_stride -1;
        dbg(5, "before if s %i, e %i, window_size %i, row_stride %i, is_writer %i, rank %i", w_row_start_index, w_row_end_index,window_size, row_stride, is_writer, rank);

        if(w_row_end_index <= GRID_SIZE_Y-1)
        {
            is_writer = 1;
        }
        else if(w_row_start_index < GRID_SIZE_Y-1 && w_row_end_index > GRID_SIZE_Y-1)
        {
            is_writer = 1;
            w_row_end_index = GRID_SIZE_Y-1;
        }
        else // left over processes that are neither readers or writers
        {
            is_writer = w_row_start_index = w_row_start_index = 0;
        }
        dbg(5, "s %i, e %i, window_size %i, row_stride %i, is_writer %i, rank %i", w_row_start_index, w_row_end_index,window_size, row_stride, is_writer, rank);

        if (is_writer)
        {
            int row_count = w_row_end_index - w_row_start_index + 1;
            interp = (GridPoint**) malloc (sizeof(GridPoint *) * row_count);

            if (interp == NULL)
            {
                cerr << "MpiInterp::init() malloc error" << endl;
                return -1;
            }
            for (i = 0; i < row_count; i++)
            {
                interp[i] = (GridPoint *) malloc (
                        sizeof(GridPoint) * GRID_SIZE_X);
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
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // all processes need to know which are readers and which are writers
    MPI_Allgather(&is_reader, 1, MPI_INT, readers, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&is_writer, 1, MPI_INT, writers, 1, MPI_INT, MPI_COMM_WORLD);
    for(i=0; i<process_count; i++)
    {
        // writer count can be less than process_count - reader_count due to the way rows are partitioned to writers
        // that is, there can be processes with rank >=  reader_count + writer_count that have nothing to do
        if(writers[i])
        {

            writer_count++;

        }
        // readers will set read_done[rank] = 1 when they are finished reading
        read_done[i] = 0;
    }
    // allocate the buffers to send and receive points, use allocation max of 1000 Megs on readers
    // limit writer allocation to 100 Megs
    mpi_point_buffer_count = 1000000000 / sizeof(grid_point_info)/ writer_count;
    if((mpi_point_buffer_count*sizeof(grid_point_info))> 100000000)
    {
        mpi_point_buffer_count =  100000000 / sizeof(grid_point_info);
    }


    dbg(5, "mpi_point_buffer_count - %li\n", mpi_point_buffer_count);

    if(is_reader)
    {
        point_buffers = (grid_point_info **) malloc (writer_count * sizeof(grid_point_info *));
        point_buffer_counts = (long *) malloc(writer_count * sizeof(long));
        for(i=0;i<writer_count;i++)
        {
            point_buffers[i] = (grid_point_info *)malloc(mpi_point_buffer_count * sizeof(grid_point_info));
            point_buffer_counts[i] = 0;
        }
    }
    if(is_writer)
    {
        point_buffer = (grid_point_info *)malloc(mpi_point_buffer_count * sizeof(grid_point_info));
        point_buffer_count = 0;
    }


    //for(i=0; i<process_count; i++){
    //    if(readers[i]) printf("reader, rank %i %i\n", i, rank);
    //    if(writers[i]) printf("writer, rank %i %i\n", i, rank);
    //}
    if(DEBUG)
    {
        printf("MPIInterp.init() done, w_row_start_index %i,  w_row_end_index %i, row_stride %i, GRID_SIZE_Y %i, GRID_SIZE_X %i, rank %i\n", w_row_start_index,w_row_end_index,row_stride,GRID_SIZE_Y,GRID_SIZE_X,rank );
    }
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

int
MpiInterp::finish (char *outputName, int outputFormat, unsigned int outputType,
                   double *adfGeoTransform, const char* wkt)
{
    int rc;
    int i, j;
    MPI_Barrier (MPI_COMM_WORLD);
    if(timer)
    {
        if(rank == reader_count)printf("Writers processing cells...\n");
        timer->process_start = time(NULL);
    }

    //printf("finish starts, rank %i\n", rank);
    // reader_count is the first writer rank, reader ranks are 0 through reader_count-1
    int first_writer_rank = reader_count;
    int last_writer_rank = first_writer_rank + writer_count - 1;
    int row_count = w_row_end_index - w_row_start_index + 1;
    GridPoint **rows_before = NULL;
    GridPoint **rows_after = NULL;
    MPI_Status status;
    //struct tms tbuf;
    clock_t t0, t1;
    if (is_writer)
    {
        for (i = 0; i < row_count; i++)
        {
            for (j = 0; j < GRID_SIZE_X; j++)
            {
                if (interp[i][j].Zmin == DBL_MAX)
                {
                    //		interp[i][j].Zmin = NAN;
                    interp[i][j].Zmin = 0;
                }

                if (interp[i][j].Zmax == -DBL_MAX)
                {
                    //interp[i][j].Zmax = NAN;
                    interp[i][j].Zmax = 0;
                }

                if (interp[i][j].count != 0)
                {
                    interp[i][j].Zmean /= interp[i][j].count;
                    interp[i][j].empty = 1;
                }
                else
                {
                    //interp[i][j].Zmean = NAN;
                    interp[i][j].Zmean = 0;
                }

                // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
                if (interp[i][j].count != 0)
                {
                    interp[i][j].Zstd = interp[i][j].Zstd
                            / (interp[i][j].count);
                    interp[i][j].Zstd = sqrt (interp[i][j].Zstd);
                }
                else
                {
                    interp[i][j].Zstd = 0;
                }

                if (interp[i][j].sum != 0 && interp[i][j].sum != -1)
                    interp[i][j].Zidw /= interp[i][j].sum;
                else if (interp[i][j].sum == -1)
                {
                    // do nothing
                }
                else
                {
                    //interp[i][j].Zidw = NAN;
                    interp[i][j].Zidw = 0;
                }
            }
        }

        if (window_size != 0)
        {
            int window_dist = window_size / 2;

            // reader_count is the first writer rank, reader ranks are 0 through reader_count-1
            int first_writer_rank = reader_count;
            int last_writer_rank = first_writer_rank + writer_count - 1;

            if (rank == first_writer_rank && rank == last_writer_rank)
            {
                // do nothing
            }
            else if (rank == first_writer_rank && rank != last_writer_rank)
            {
                // alloc and get window_dist rows from next rank for rows_after
                rows_after = allocRows (window_dist);
                // MPI send last window_dist rows to rank +1
                //MPI recv first window_dist rows from rank + 1 into rows_after
            }
            else if (rank != first_writer_rank && rank != last_writer_rank)
            {
                // alloc and get window_dist rows from previous rank for rows_before
                rows_before = allocRows (window_dist);
                // alloc and get window_dist rows from next rank for rows_after
                rows_after = allocRows (window_dist);
            }
            else if (rank != first_writer_rank && rank == last_writer_rank)
            {
                // alloc and get window_dist rows from previous rank for rows_before
                rows_before = allocRows (window_dist);
            }
        }
    }
    MPI_Barrier (MPI_COMM_WORLD);

    if (is_writer && window_size != 0)
    {
        int i;
        int window_dist = window_size / 2;
        if (rank == first_writer_rank && rank == last_writer_rank)
        {
            // one writer process, do nothing
        }
        else if (rank == first_writer_rank && rank != last_writer_rank)
        {
           // first writer, send last rows and recv first rows from next writer
           for(i=window_dist-1; i>=0; i--)
           {
               MPI_Send (interp[w_row_end_index-w_row_start_index-i], GRID_SIZE_X, mpi_grid_point,
                         rank + 1, 1, MPI_COMM_WORLD);
           }
           for(i=0; i<window_dist; i++){
               MPI_Recv (rows_after[i], GRID_SIZE_X, mpi_grid_point, rank + 1,
                         1, MPI_COMM_WORLD, &status);
           }

        }
        else if (rank != first_writer_rank && rank != last_writer_rank)
        {
            // mid writer, recv last rows of previous writer and send first rows to previous writer
            for (i = 0; i < window_dist; i++)
            {
                MPI_Recv (rows_before[i], GRID_SIZE_X, mpi_grid_point, rank - 1,
                          1, MPI_COMM_WORLD, &status);
            }
            for (i = 0; i < window_dist; i++)
            {
                MPI_Send (interp[i], GRID_SIZE_X, mpi_grid_point, rank - 1, 1,
                          MPI_COMM_WORLD);
            }
            // mid writer, send last to next writer rows and recv first rows from next writer
            for (i = window_dist - 1; i >= 0; i--)
            {
                MPI_Send (interp[w_row_end_index - w_row_start_index - i],
                          GRID_SIZE_X, mpi_grid_point, rank + 1, 1,
                          MPI_COMM_WORLD);
            }
            for (i = 0; i < window_dist; i++)
            {
                MPI_Recv (rows_after[i], GRID_SIZE_X, mpi_grid_point, rank + 1,
                          1, MPI_COMM_WORLD, &status);
            }

        }
        else if (rank != first_writer_rank && rank == last_writer_rank)
        {
            // last writer, recv last rows of previous writer and send first rows to previous writer
            for(i=0; i<window_dist; i++)
            {
                MPI_Recv (rows_before[i], GRID_SIZE_X, mpi_grid_point,
                          rank - 1, 1, MPI_COMM_WORLD, &status);
            }
            for(i=0; i<window_dist; i++)
            {
                MPI_Send(interp[i], GRID_SIZE_X,
                          mpi_grid_point, rank - 1, 1, MPI_COMM_WORLD);
            }

        }
    }
    // Sriram's edit: Fill zeros using the window size parameter
    MPI_Barrier (MPI_COMM_WORLD);
    if (is_writer)
    {
        if (window_size != 0)
        {
            int window_dist = window_size / 2;
            for (int i = 0; i < row_count; i++)
            {
                for (int j = 0; j < GRID_SIZE_X; j++)
                {
                    if (interp[i][j].empty == 0)
                    {
                        double new_sum = 0.0;
                        for (int p = i - window_dist; p <= i + window_dist; p++)
                        {
                            for (int q = j - window_dist; q <= j + window_dist; q++)
                            {

                                if ((p == i) && (q == j))
                                {
                                    continue;
                                }

                                if ((p >= 0) && (p < row_count) && (q >= 0)
                                        && (q < GRID_SIZE_X))
                                {
                                    if (interp[p][q].empty != 0)
                                    {
                                        double distance = max (abs (p - i),
                                                               abs (q - j));
                                        interp[i][j].Zmean +=
                                                interp[p][q].Zmean
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zidw +=
                                                interp[p][q].Zidw
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zstd +=
                                                interp[p][q].Zstd
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zstd_tmp +=
                                                interp[p][q].Zstd_tmp
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zmin +=
                                                interp[p][q].Zmin
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zmax +=
                                                interp[p][q].Zmax
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));

                                        new_sum +=
                                                1
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                    }
                                }

                                else if ((p < 0) && (q >= 0) && (q < GRID_SIZE_X) && rows_before != NULL)
                                {
                                    int p2 = p + window_dist;
                                    if (rows_before[p2][q].empty != 0)
                                    {
                                        //printf("rows before %i %i %i %i %i %i %f %i %i \n", i, j, p2, p, q, rank, rows_before[p2][q].Zmean, w_row_start_index, w_row_end_index);
                                        double distance = max (abs (p - i),
                                                               abs (q - j));
                                        interp[i][j].Zmean +=
                                                rows_before[p2][q].Zmean
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zidw +=
                                                rows_before[p2][q].Zidw
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zstd +=
                                                rows_before[p2][q].Zstd
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zstd_tmp +=
                                                rows_before[p2][q].Zstd_tmp
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zmin +=
                                                rows_before[p2][q].Zmin
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zmax +=
                                                rows_before[p2][q].Zmax
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));

                                        new_sum +=
                                                1
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                    }
                                }

                                else if ((p >= row_count) && (q >= 0)
                                        && (q < GRID_SIZE_X)
                                        && rows_after != NULL)
                                {
                                    int p2 = p-row_count;
                                    if (rows_after[p2][q].empty != 0)
                                    {
                                        //printf("rows after %i % i %i %i %i %i %f %i %i\n", i, j, p2, p, q, rank, rows_after[p2][q].Zmean, w_row_start_index, w_row_end_index);
                                        double distance = max (abs (p - i),
                                                               abs (q - j));
                                        interp[i][j].Zmean +=
                                                rows_after[p2][q].Zmean
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zidw +=
                                                rows_after[p2][q].Zidw
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zstd +=
                                                rows_after[p2][q].Zstd
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zstd_tmp +=
                                                rows_after[p2][q].Zstd_tmp
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zmin +=
                                                rows_after[p2][q].Zmin
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                        interp[i][j].Zmax +=
                                                rows_after[p2][q].Zmax
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));

                                        new_sum +=
                                                1
                                                        / (pow (distance,
                                                                Interpolation::WEIGHTER));
                                    }
                                }




                            }
                        }
                        if (new_sum > 0)
                        {
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
        }
    }

    if(timer)timer->process_end = time(NULL);
    t0 = clock ();
    MPI_Barrier (MPI_COMM_WORLD);


    if(timer)
    {
        if(rank == reader_count)printf("Writers writing cells...\n");
        timer->output_start = time(NULL);
    }

    if ((rc = outputFile (outputName, outputFormat, outputType, adfGeoTransform,
                          wkt)) < 0)
    {
        cerr << "MpiInterp::finish outputFile error" << endl;
        return -1;
    }
    if(timer)timer->output_end = time(NULL);
    t1 = clock ();

    //cerr << "Output Execution time: " << (double)(t1 - t0)/ CLOCKS_PER_SEC << std::endl;

    return 0;
}

//////////////////////////////////////////////////////
// Private Methods
//////////////////////////////////////////////////////

GridPoint **
MpiInterp::allocRows (int cnt)
{

    //printf("alloc rows called %i, %i\n", cnt, rank);
    GridPoint **rows;
    rows = (GridPoint**) malloc (sizeof(GridPoint *) * cnt);

    if (rows == NULL)
    {
        cerr << "MpiInterp::getRows() malloc error" << endl;
        return NULL;
    }
    int i;
    for (i = 0; i < cnt; i++)
    {
        rows[i] = (GridPoint *) malloc (sizeof(GridPoint) * GRID_SIZE_X);
        if (rows[i] == NULL)
        {
            cerr << "MpiInterp::getRows() malloc error" << endl;
            return NULL;
        }
    }
    return rows;
}




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
                int target_rank = get_target_rank(j);
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
                int target_rank = get_target_rank(j);
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
                int target_rank = get_target_rank(j);
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
                int target_rank = get_target_rank(j);
                updateGridPointSend(target_rank, i, j, data_z, sqrt(distance));
            } else if (j == base_y) {
                return ;
            } else {
                break;
            }
        }
    }
}

int MpiInterp::get_target_rank(int row_index){

    int rtn = row_index/row_stride + reader_count;
    if (rtn < process_count)
    {
        return rtn;
    }
    else
    {
        fprintf(stderr, "MpiInterp::get_target_rank: worker process overflow, logic error, row_index %i row stride %i, w_row_end_index %i, rank %i\n", row_index, row_stride, w_row_end_index, rank);
        return -1;
    }
}


void MpiInterp::updateGridPointSend(int target_rank, int x, int y, double data_z, double distance)
{

    int i = target_rank - reader_count;
    grid_point_info info;
    info.comm = read_done[rank];
    info.x = x;
    info.y = y;
    info.data_z = data_z;
    info.distance = distance;

    point_buffers[i][point_buffer_counts[i]] = info;
    point_buffer_counts[i]++;

    if(info.comm || (!info.comm && point_buffer_counts[i]==mpi_point_buffer_count))
    {
        MPI_Send(&point_buffer_counts[i], 1, MPI_INT, target_rank, 1, MPI_COMM_WORLD);
        MPI_Send(point_buffers[i], point_buffer_counts[i], mpi_grid_point_info, target_rank, 1, MPI_COMM_WORLD);
        point_buffer_counts[i] = 0;
    }

}


void MpiInterp::updateGridPointRecv()
{
    MPI_Status status;
    int comm_done;
    int i;
    grid_point_info info;

    while(true){

        MPI_Recv(&point_buffer_count, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(point_buffer, point_buffer_count, mpi_grid_point_info, status.MPI_SOURCE, 1, MPI_COMM_WORLD, &status);

        for(i=0; i<point_buffer_count; i++)
        {
            if(point_buffer[i].comm)
            {
                read_done[status.MPI_SOURCE] = 1;
            }
            else
            {
                updateGridPoint(point_buffer[i].x, point_buffer[i].y, point_buffer[i].data_z, point_buffer[i].distance);
            }
        }

        // determine if all the readers are done sending points
        comm_done = 1;
        for(i=0; i<reader_count;i++)
        {
            if(read_done[i] == 0){
                comm_done = 0;
            }
        }
        if(comm_done)
        {
            break;
        }
    }

}


void MpiInterp::updateGridPoint(int x, int y, double data_z, double distance)
{

    static int count = 0;

    //count++;

    //if(!(count%100000)){
    //    printf("updateGridPoint, rank %i called %i times\n", rank, count);
    //}

    //printf("update starts, rank %i, count %i\n", rank, count);
    //printf ("updateGridPoint rank %i, x %i, y %i, count %i\n", rank, x, y, count);

    //todo update to work with more than one reader
    //int reader_count = 1;

    y -= (rank-reader_count)*row_stride;

    if(x<0 || x>=GRID_SIZE_X|| y<0 || y>= (w_row_end_index - w_row_start_index + 1 ))
    {
        printf ("ERROR in updateGridPoint rank %i, y %i, x %i\n", rank, y, x);
    }



    if(interp[y][x].Zmin > data_z)
        interp[y][x].Zmin = data_z;
    if(interp[y][x].Zmax < data_z)
        interp[y][x].Zmax = data_z;

    interp[y][x].Zmean += data_z;
    interp[y][x].count++;

    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
    double delta = data_z - interp[y][x].Zstd_tmp;
    interp[y][x].Zstd_tmp += delta/interp[y][x].count;
    interp[y][x].Zstd += delta * (data_z - interp[y][x].Zstd_tmp);

    double dist = pow(distance, Interpolation::WEIGHTER);

    if(interp[y][x].sum != -1) {
        if(dist != 0) {
            interp[y][x].Zidw += data_z/dist;
            interp[y][x].sum += 1/dist;
        } else {
            interp[y][x].Zidw = data_z;
            interp[y][x].sum = -1;
        }
    } else {
        // do nothing
    }

   // printf("update ends, rank %i\n", rank);
}

void MpiInterp::printArray()
{
    int i, j;
    if (is_writer)
    {
        for (i = 0; i < w_row_end_index - w_row_start_index + 1; i++)
        {
            for (j = 1; j < GRID_SIZE_X; j++)
            {
                cerr << interp[i][j].Zmin << ", " << interp[i][j].Zmax << ", ";
                cerr << interp[i][j].Zmean << ", " << interp[i][j].Zidw << endl;
                //printf("%.20f ", interp[i][j].Zmax);
            }
            //printf("\n");
        }
    }
    cerr << endl;
}


void MpiInterp::flushMpiBuffers (MPI_File *arcFiles, MPI_File *gridFiles,
                            int numTypes, long max_buffer_size = 0)
{
    int k = 0;
    MPI_Status status;

    if (arcFiles != NULL)
    {
        for (k = 0; k < numTypes; k++)
        {
            if (arcFiles[k] != NULL)
            {
                if (arc_file_mpi_count[k] > max_buffer_size)
                {
                    MPI_File_write (arcFiles[k], arc_file_mpi_buffer[k],
                                    arc_file_mpi_count[k], MPI_CHAR, &status);

                    arc_file_mpi_buffer[k][0] = 0;
                    arc_file_mpi_count[k] = 0;

                }
            }
        }
    }

    if (gridFiles != NULL)
    {
        for (k = 0; k < numTypes; k++)
        {
            if (gridFiles[k] != NULL)
            {
                if (grid_file_mpi_count[k] > max_buffer_size)
                {
                    MPI_File_write (gridFiles[k], grid_file_mpi_buffer[k],
                                    grid_file_mpi_count[k], MPI_CHAR, &status);
                    //printf("grid %i %lu %i\n", grid_file_mpi_count[k], strlen(grid_file_mpi_buffer[k]), rank);
                    grid_file_mpi_buffer[k][0] = 0;
                    grid_file_mpi_count[k] = 0;
                }
            }
        }
    }

}


int
MpiInterp::outputFile (char *outputName, int outputFormat,
                       unsigned int outputType, double *adfGeoTransform,
                       const char* wkt)
{
    int i, j, k;
    // reader_count is the first writer rank, reader ranks are 0 through reader_count-1
    int first_writer_rank = reader_count;
    int last_writer_rank = first_writer_rank + writer_count - 1;
    //FILE **arcFiles;
    MPI_File *arcFiles;
    char arcFileName[1024];

    MPI_File *gridFiles;
    char gridFileName[1024];

    char buf[1024];
    MPI_Status status;

    const char *ext[6] = { ".min", ".max", ".mean", ".idw", ".den", ".std" };
    unsigned int type[6] = { OUTPUT_TYPE_MIN, OUTPUT_TYPE_MAX, OUTPUT_TYPE_MEAN,
            OUTPUT_TYPE_IDW, OUTPUT_TYPE_DEN, OUTPUT_TYPE_STD };
    int numTypes = 6;


    if (outputFormat == OUTPUT_FORMAT_ARC_ASCII
            || outputFormat == OUTPUT_FORMAT_ALL)
    {
        if ((arcFiles = (MPI_File *) malloc (sizeof(MPI_File) * numTypes))
                == NULL)
        {
            cerr << "Arc MPI_File malloc error: " << endl;
            return -1;
        }
        arc_file_mpi_offset = (MPI_Offset *) malloc (
                sizeof(MPI_Offset) * numTypes); // write position after header
        arc_file_mpi_buffer = (char **) malloc (sizeof(char *) * numTypes);
        arc_file_mpi_count = (long *) malloc (sizeof(long) * numTypes);
        arc_file_mpi_size = (long *) malloc (
                sizeof(long) * numTypes);
        arc_file_mpi_sizes = (long **) malloc (
                sizeof(long *) * numTypes);

        for (i = 0; i < numTypes; i++)
        {
            if (outputType & type[i])
            {
                strncpy (arcFileName, outputName, sizeof(arcFileName));
                strncat (arcFileName, ext[i], strlen (ext[i]));
                strncat (arcFileName, ".asc", strlen (".asc"));
                MPI_File_open (MPI_COMM_WORLD, arcFileName,
                MPI_MODE_CREATE | MPI_MODE_WRONLY,
                               MPI_INFO_NULL, &(arcFiles[i]));
                MPI_File_set_size (arcFiles[i], 0);
                arc_file_mpi_offset[i] = 0;
                arc_file_mpi_buffer[i] = (char *) malloc (
                        sizeof(char) * buffer_size);
                arc_file_mpi_sizes[i] = (long *) malloc (
                        sizeof(long) * process_count);
                arc_file_mpi_count[i] = 0;
                arc_file_mpi_size[i] = 0;
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

    if (outputFormat == OUTPUT_FORMAT_GRID_ASCII
            || outputFormat == OUTPUT_FORMAT_ALL)
    {
        if ((gridFiles = (MPI_File *) malloc (sizeof(MPI_File) * numTypes))
                == NULL)
        {
            cerr << "File array allocation error" << endl;
            return -1;
        }
        grid_file_mpi_offset = (MPI_Offset *) malloc (
                sizeof(MPI_Offset) * numTypes);
        grid_file_mpi_buffer = (char **) malloc (sizeof(char *) * numTypes);
        grid_file_mpi_count = (long *) malloc (sizeof(long) * numTypes);
        grid_file_mpi_size = (long *) malloc (
                sizeof(long) * numTypes);
        grid_file_mpi_sizes = (long **) malloc (
                sizeof(long *) * numTypes);

        for (i = 0; i < numTypes; i++)
        {
            if (outputType & type[i])
            {
                strncpy (gridFileName, outputName, sizeof(arcFileName));
                strncat (gridFileName, ext[i], strlen (ext[i]));
                strncat (gridFileName, ".grid", strlen (".grid"));

                MPI_File_open (MPI_COMM_WORLD, gridFileName,
                MPI_MODE_CREATE | MPI_MODE_WRONLY,
                               MPI_INFO_NULL, &(gridFiles[i]));
                MPI_File_set_size (gridFiles[i], 0);
                grid_file_mpi_offset[i] = 0;
                grid_file_mpi_buffer[i] = (char *) malloc (
                        sizeof(char) * buffer_size);
                grid_file_mpi_sizes[i] = (long *) malloc (
                        sizeof(long) * process_count);
                grid_file_mpi_count[i] = 0;
                grid_file_mpi_size[i] = 0;
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
    MPI_Barrier (MPI_COMM_WORLD);

    if (rank == first_writer_rank)
    {
        // print ArcGIS headers
        if (arcFiles != NULL)
        {
            for (i = 0; i < numTypes; i++)
            {
                if (arcFiles[i] != NULL)
                {
                    buf[0] = 0;
                    sprintf (buf, "ncols %d\n", GRID_SIZE_X);
                    MPI_File_write (arcFiles[i], buf, strlen (buf), MPI_CHAR,
                                    &status);
                    buf[0] = 0;
                    sprintf (buf, "nrows %d\n", GRID_SIZE_Y);
                    MPI_File_write (arcFiles[i], buf, strlen (buf), MPI_CHAR,
                                    &status);
                    buf[0] = 0;
                    sprintf (buf, "xllcorner %f\n", min_x);
                    MPI_File_write (arcFiles[i], buf, strlen (buf), MPI_CHAR,
                                    &status);
                    buf[0] = 0;
                    sprintf (buf, "yllcorner %f\n", min_y);
                    MPI_File_write (arcFiles[i], buf, strlen (buf), MPI_CHAR,
                                    &status);
                    buf[0] = 0;
                    sprintf (buf, "cellsize %f\n", GRID_DIST_X);
                    MPI_File_write (arcFiles[i], buf, strlen (buf), MPI_CHAR,
                                    &status);
                    buf[0] = 0;
                    sprintf (buf, "NODATA_value -9999\n");
                    MPI_File_write (arcFiles[i], buf, strlen (buf), MPI_CHAR,
                                    &status);
                }
                MPI_File_get_position (arcFiles[i], &(arc_file_mpi_offset[i]));
            }
        }

        // print Grid headers
        if (gridFiles != NULL)
        {
            for (i = 0; i < numTypes; i++)
            {
                if (gridFiles[i] != NULL)
                {
                    buf[0] = 0;
                    sprintf (buf, "north: %f\n", max_y);
                    MPI_File_write (gridFiles[i], buf, strlen (buf), MPI_CHAR,
                                    &status);
                    buf[0] = 0;
                    sprintf (buf, "south: %f\n", min_y);
                    MPI_File_write (gridFiles[i], buf, strlen (buf), MPI_CHAR,
                                    &status);
                    buf[0] = 0;
                    sprintf (buf, "east: %f\n", max_x);
                    MPI_File_write (gridFiles[i], buf, strlen (buf), MPI_CHAR,
                                    &status);
                    buf[0] = 0;
                    sprintf (buf, "west: %f\n", min_x);
                    MPI_File_write (gridFiles[i], buf, strlen (buf), MPI_CHAR,
                                    &status);
                    buf[0] = 0;
                    sprintf (buf, "rows: %d\n", GRID_SIZE_Y);
                    MPI_File_write (gridFiles[i], buf, strlen (buf), MPI_CHAR,
                                    &status);
                    buf[0] = 0;
                    sprintf (buf, "cols: %d\n", GRID_SIZE_X);
                    MPI_File_write (gridFiles[i], buf, strlen (buf), MPI_CHAR,
                                    &status);

                }
                MPI_File_get_position (gridFiles[i],
                                       &(grid_file_mpi_offset[i]));
            }
        }
    } // if(rank == first_writer_rank)

    MPI_Barrier (MPI_COMM_WORLD);

    if (arcFiles != NULL)
    {
        for (i = 0; i < numTypes; i++)
        {
            if (arcFiles[i] != NULL)
            {
                MPI_Bcast (&(arc_file_mpi_offset[i]), 1, MPI_OFFSET, first_writer_rank,
                MPI_COMM_WORLD);
            }
        }
    }

    if (gridFiles != NULL)
    {
        for (i = 0; i < numTypes; i++)
        {
            if (gridFiles[i] != NULL)
            {
                MPI_Bcast (&(grid_file_mpi_offset[i]), 1, MPI_OFFSET, first_writer_rank,
                MPI_COMM_WORLD);
            }
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);

    if (arcFiles != NULL)
    {
        for (i = 0; i < numTypes; i++)
        {
            if (arcFiles[i] != NULL)
            {
                MPI_File_seek (arcFiles[i], arc_file_mpi_offset[i],
                MPI_SEEK_SET);
            }
        }
    }

    if (gridFiles != NULL)
    {
        for (i = 0; i < numTypes; i++)
        {
            if (gridFiles[i] != NULL)
            {
                MPI_File_seek (gridFiles[i], grid_file_mpi_offset[i],
                MPI_SEEK_SET);
            }
        }
    }
    for (i = 0; i < 6; i++)
    {
        //printf("arc and grid_file_mpi_offset %lli %lli %i\n", arc_file_mpi_offset[i], grid_file_mpi_offset[i], rank);
    }

    if (is_writer)
    {
        //**************************** Loop 1, calculate total write size for each worker and set file offset *********
        // initialize mpi_size for each file
        if (arcFiles != NULL)
            for (k = 0; k < numTypes; k++)
            {
                if (arcFiles[k] != NULL)
                    arc_file_mpi_size[k] = 0;
            }
        if (gridFiles != NULL)
            for (k = 0; k < numTypes; k++)
            {
                if (gridFiles[k] != NULL)
                    grid_file_mpi_size[k] = 0;
            }

        // add up the write size for all files for this worker
        int row_count = w_row_end_index - w_row_start_index + 1;
        for (i = row_count - 1; i >= 0; i--)
        {
            for (j = 0; j < GRID_SIZE_X; j++)
            {
                if (arcFiles != NULL)
                {
                    // Zmin
                    if (arcFiles[0] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            arc_file_mpi_size[0] += sprintf (buf, "-9999 ");
                        }
                        else
                        {
                            arc_file_mpi_size[0] += sprintf (buf, "%f ",
                                                             interp[i][j].Zmin);
                        }
                    }

                    // Zmax
                    if (arcFiles[1] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            arc_file_mpi_size[1] += sprintf (buf, "-9999 ");
                        }
                        else
                        {
                            arc_file_mpi_size[1] += sprintf (buf, "%f ",
                                                             interp[i][j].Zmax);
                        }
                    }

                    // Zmean
                    if (arcFiles[2] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            arc_file_mpi_size[2] += sprintf (buf, "-9999 ");
                        }
                        else
                        {
                            arc_file_mpi_size[2] += sprintf (
                                    buf, "%f ", interp[i][j].Zmean);
                        }
                    }

                    // Zidw
                    if (arcFiles[3] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            arc_file_mpi_size[3] += sprintf (buf, "-9999 ");
                        }
                        else
                        {
                            arc_file_mpi_size[3] += sprintf (buf, "%f ",
                                                             interp[i][j].Zidw);
                        }
                    }

                    // count
                    if (arcFiles[4] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            arc_file_mpi_size[4] += sprintf (buf, "-9999 ");
                        }
                        else
                        {
                            arc_file_mpi_size[4] += sprintf (
                                    buf, "%d ", interp[i][j].count);
                        }
                    }

                    // count
                    if (arcFiles[5] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            arc_file_mpi_size[5] += sprintf (buf, "-9999 ");
                        }
                        else
                        {
                            arc_file_mpi_size[5] += sprintf (buf, "%f ",
                                                             interp[i][j].Zstd);
                        }
                    }
                }

                if (gridFiles != NULL)
                {
                    // Zmin
                    if (gridFiles[0] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            grid_file_mpi_size[0] += sprintf (buf, "-9999 ");
                        }
                        else
                        {
                            grid_file_mpi_size[0] += sprintf (
                                    buf, "%f ", interp[i][j].Zmin);
                        }
                    }

                    // Zmax
                    if (gridFiles[1] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            grid_file_mpi_size[1] += sprintf (buf, "-9999 ");
                        }
                        else
                        {
                            grid_file_mpi_size[1] += sprintf (
                                    buf, "%f ", interp[i][j].Zmax);
                        }
                    }

                    // Zmean
                    if (gridFiles[2] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            grid_file_mpi_size[2] += sprintf (buf, "-9999 ");
                        }
                        else
                        {
                            grid_file_mpi_size[2] += sprintf (
                                    buf, "%f ", interp[i][j].Zmean);
                        }
                    }

                    // Zidw
                    if (gridFiles[3] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            grid_file_mpi_size[3] += sprintf (buf, "-9999 ");
                        }
                        else
                        {
                            grid_file_mpi_size[3] += sprintf (
                                    buf, "%f ", interp[i][j].Zidw);
                        }
                    }

                    // count
                    if (gridFiles[4] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            grid_file_mpi_size[4] += sprintf (buf, "-9999 ");
                        }
                        else
                        {
                            grid_file_mpi_size[4] += sprintf (
                                    buf, "%d ", interp[i][j].count);
                        }
                    }

                    // count
                    if (gridFiles[5] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            grid_file_mpi_size[5] += sprintf (buf, "-9999 ");
                        }
                        else
                        {
                            grid_file_mpi_size[5] += sprintf (
                                    buf, "%f ", interp[i][j].Zstd);
                        }
                    }
                }

                if (j == GRID_SIZE_X - 1) // put a newline after each row
                {
                    for (k = 0; k < numTypes; k++)
                    {
                        if (arcFiles != NULL && arcFiles[k] != NULL)
                        {
                            arc_file_mpi_size[k] += sprintf(buf, "\n");
                        }
                        if (gridFiles != NULL && gridFiles[k] != NULL)
                        {
                            grid_file_mpi_size[k] += sprintf(buf, "\n");
                        }
                    }
                }
            }
        }

    } // end if(is_writer)

    // ********** Gather, calculate, and set each worker's offset **********
    MPI_Barrier (MPI_COMM_WORLD);

    if (arcFiles != NULL)
    {
        for (i = 0; i < numTypes; i++)
        {
            if (arcFiles[i] != NULL)
            {
                MPI_Allgather (&(arc_file_mpi_size[i]), 1, MPI_LONG,
                               arc_file_mpi_sizes[i], 1, MPI_LONG,
                               MPI_COMM_WORLD);
            }
        }
    }

    if (gridFiles != NULL)
    {
        for (i = 0; i < numTypes; i++)
        {
            if (gridFiles[i] != NULL)
            {
                MPI_Allgather (&(grid_file_mpi_size[i]), 1, MPI_LONG,
                               grid_file_mpi_sizes[i], 1, MPI_LONG,
                               MPI_COMM_WORLD);
            }
        }
    }

    if (is_writer)
    {
        if (arcFiles != NULL)
        {
            for (i = 0; i < numTypes; i++)
            {
                if (arcFiles[i] != NULL)
                {
                    int j;
                    for (j = last_writer_rank; j > rank; j--)
                    {
                        arc_file_mpi_offset[i] += arc_file_mpi_sizes[i][j];
                    }
                    MPI_File_seek (arcFiles[i], arc_file_mpi_offset[i],
                    MPI_SEEK_SET);
                }
            }
        }

        if (gridFiles != NULL)
        {
            for (i = 0; i < numTypes; i++)
            {
                if (gridFiles[i] != NULL)
                {
                    int j;
                    for (j = last_writer_rank; j > rank; j--)
                    {
                        grid_file_mpi_offset[i] += grid_file_mpi_sizes[i][j];
                    }
                    MPI_File_seek (gridFiles[i], grid_file_mpi_offset[i],
                    MPI_SEEK_SET);
                    //printf ("grid_file_mpi_sizes[i] %i %u, %lli,  rank %i\n", i,
                    //        grid_file_mpi_sizes[i][rank],
                    //        grid_file_mpi_offset[i], rank);
                }
            }
        }

        // ********** Write the file **********
        // initialize mpi_size for each file
        if (arcFiles != NULL)
        {
            for (k = 0; k < numTypes; k++)
            {
                if (arcFiles[k] != NULL)
                {
                    arc_file_mpi_buffer[k][0] = 0;
                    arc_file_mpi_count[k] = 0;
                }
            }
        }
        if (gridFiles != NULL)
        {
            for (k = 0; k < numTypes; k++)
            {
                if (gridFiles[k] != NULL)
                {
                    grid_file_mpi_buffer[k][0] = 0;
                    grid_file_mpi_count[k] = 0;
                }
            }
        }
        // write the file
        int row_count = w_row_end_index - w_row_start_index + 1;
        for (i = row_count - 1; i >= 0; i--)
        {
            for (j = 0; j < GRID_SIZE_X; j++)
            {

                if (arcFiles != NULL)
                {
                    // Zmin
                    if (arcFiles[0] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            arc_file_mpi_count[0] += sprintf (
                                    arc_file_mpi_buffer[0]+arc_file_mpi_count[0], "-9999 ");
                        }
                        else
                        {
                            arc_file_mpi_count[0] += sprintf (
                                    arc_file_mpi_buffer[0]+ arc_file_mpi_count[0], "%f ",
                                    interp[i][j].Zmin);
                        }

                        //printf("%s\n", arc_file_mpi_buffer[0]);
                    }

                    // Zmax
                    if (arcFiles[1] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            arc_file_mpi_count[1] += sprintf (
                                    arc_file_mpi_buffer[1]+arc_file_mpi_count[1], "-9999 ");
                        }
                        else
                        {
                            arc_file_mpi_count[1] += sprintf (
                                    arc_file_mpi_buffer[1]+arc_file_mpi_count[1], "%f ",
                                    interp[i][j].Zmax);
                        }
                    }

                    // Zmean
                    if (arcFiles[2] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            arc_file_mpi_count[2] += sprintf (
                                    arc_file_mpi_buffer[2]+arc_file_mpi_count[2], "-9999 ");
                        }
                        else
                        {
                            arc_file_mpi_count[2] += sprintf (
                                    arc_file_mpi_buffer[2]+arc_file_mpi_count[2], "%f ",
                                    interp[i][j].Zmean);
                        }
                    }

                    // Zidw
                    if (arcFiles[3] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            arc_file_mpi_count[3] += sprintf (
                                    arc_file_mpi_buffer[3]+arc_file_mpi_count[3], "-9999 ");
                        }
                        else
                        {
                            arc_file_mpi_count[3] += sprintf (
                                    arc_file_mpi_buffer[3]+arc_file_mpi_count[3], "%f ",
                                    interp[i][j].Zidw);
                        }
                    }

                    // count
                    if (arcFiles[4] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            arc_file_mpi_count[4] += sprintf (
                                    arc_file_mpi_buffer[4]+arc_file_mpi_count[4], "-9999 ");
                        }
                        else
                        {
                            arc_file_mpi_count[4] += sprintf (
                                    arc_file_mpi_buffer[4]+arc_file_mpi_count[4], "%d ",
                                    interp[i][j].count);
                        }
                    }

                    // count
                    if (arcFiles[5] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            arc_file_mpi_count[5] += sprintf (
                                    arc_file_mpi_buffer[5]+arc_file_mpi_count[5], "-9999 ");
                        }
                        else
                        {
                            arc_file_mpi_count[5] += sprintf (
                                    arc_file_mpi_buffer[5]+arc_file_mpi_count[5], "%f ",
                                    interp[i][j].Zstd);
                        }
                    }
                }

                if (gridFiles != NULL)
                {
                    // Zmin
                    if (gridFiles[0] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            grid_file_mpi_count[0] += sprintf (
                                    grid_file_mpi_buffer[0]+grid_file_mpi_count[0],"-9999 ");
                        }
                        else
                        {
                            grid_file_mpi_count[0] += sprintf (
                                    grid_file_mpi_buffer[0]+grid_file_mpi_count[0],"%f ",
                                    interp[i][j].Zmin);
                        }
                    }

                    // Zmax
                    if (gridFiles[1] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            grid_file_mpi_count[1] += sprintf (
                                    grid_file_mpi_buffer[1]+grid_file_mpi_count[1], "-9999 ");
                        }
                        else
                        {
                            grid_file_mpi_count[1] += sprintf (
                                    grid_file_mpi_buffer[1]+grid_file_mpi_count[1], "%f ",
                                    interp[i][j].Zmax);
                        }
                    }

                    // Zmean
                    if (gridFiles[2] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            grid_file_mpi_count[2] += sprintf (
                                    grid_file_mpi_buffer[2]+grid_file_mpi_count[2], "-9999 ");
                        }
                        else
                        {
                            grid_file_mpi_count[2] += sprintf (
                                    grid_file_mpi_buffer[2]+grid_file_mpi_count[2], "%f ",
                                    interp[i][j].Zmean);
                        }
                    }

                    // Zidw
                    if (gridFiles[3] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            grid_file_mpi_count[3] += sprintf (
                                    grid_file_mpi_buffer[3]+grid_file_mpi_count[3], "-9999 ");
                        }
                        else
                        {
                            grid_file_mpi_count[3] += sprintf (
                                    grid_file_mpi_buffer[3]+grid_file_mpi_count[3], "%f ",
                                    interp[i][j].Zidw);
                        }
                    }

                    // count
                    if (gridFiles[4] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            grid_file_mpi_count[4] += sprintf (
                                    grid_file_mpi_buffer[4]+grid_file_mpi_count[4], "-9999 ");
                        }
                        else
                        {
                            grid_file_mpi_count[4] += sprintf (
                                    grid_file_mpi_buffer[4]+grid_file_mpi_count[4], "%d ",
                                    interp[i][j].count);
                        }
                    }

                    // count
                    if (gridFiles[5] != NULL)
                    {
                        if (interp[i][j].empty == 0 && interp[i][j].filled == 0)
                        {
                            grid_file_mpi_count[5] += sprintf (
                                    grid_file_mpi_buffer[5]+grid_file_mpi_count[5], "-9999 ");
                        }
                        else
                        {
                            grid_file_mpi_count[5] += sprintf (
                                    grid_file_mpi_buffer[5]+grid_file_mpi_count[5], "%f ",
                                    interp[i][j].Zstd);
                        }
                    }
                }


                if (j == GRID_SIZE_X - 1)
                {
                    for (k = 0; k < numTypes; k++)
                    {

                        if (arcFiles != NULL && arcFiles[k] != NULL)
                        {
                            arc_file_mpi_count[k] += sprintf (
                                    arc_file_mpi_buffer[k] + arc_file_mpi_count[k],
                                    "\n");
                        }
                        if (gridFiles != NULL && gridFiles[k] != NULL)
                        {
                            grid_file_mpi_count[k] += sprintf (
                                    grid_file_mpi_buffer[k]
                                            + grid_file_mpi_count[k],
                                    "\n");
                        }

                    }
                }
                flushMpiBuffers (arcFiles, gridFiles, numTypes, buffer_size - 1024);
            }
        }

        flushMpiBuffers (arcFiles, gridFiles, numTypes, 0);

    } // if (is_writer)

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef HAVE_GDAL
    MPI_File *tifFiles = NULL;

    if (outputFormat == OUTPUT_FORMAT_GDAL_GTIFF
            || outputFormat == OUTPUT_FORMAT_ALL)
    {

        char tifFileName[1024];
        char tifTemplateName[1024];

        strcpy (tifTemplateName, outputName);
        strcat (tifTemplateName, ".template.tif");


        unsigned char *tifTemplateContents;
        unsigned long  tifTemplateCount;

        tifFiles = (MPI_File *) malloc (sizeof(MPI_File) * numTypes);
        tif_file_mpi_offset = (MPI_Offset *) malloc (
                sizeof(MPI_Offset) * numTypes); // write position after header
        tif_file_mpi_size = (long *) malloc (sizeof(long) * numTypes);
        tif_file_mpi_sizes = (long **) malloc (sizeof(long *) * numTypes);

        for (i = 0; i < numTypes; i++)
        {
            if (outputType & type[i])
            {
                strncpy (tifFileName, outputName, sizeof(tifFileName));
                strncat (tifFileName, ext[i], strlen (ext[i]));
                strncat (tifFileName, ".tif", strlen (".tif"));
                MPI_File_open (MPI_COMM_WORLD, tifFileName,
                MPI_MODE_CREATE | MPI_MODE_WRONLY,
                               MPI_INFO_NULL, &(tifFiles[i]));
                MPI_File_set_size (tifFiles[i], 0);
                tif_file_mpi_offset[i] = 0;
                tif_file_mpi_sizes[i] = (long *) malloc (
                        sizeof(long) * process_count);
                tif_file_mpi_size[i] = 0;
            }
            else
            {
                tifFiles[i] = NULL;
            }
        }

        // create tifTemplate containing tiff header and directory, correct the strip offsets and sizes, all with first writer process

        if (rank == first_writer_rank)
        {
            GDALAllRegister ();
            char **papszMetadata;
            const char *pszFormat = "GTIFF";
            GDALDriver* tpDriver = GetGDALDriverManager ()->GetDriverByName (
                    pszFormat);
            papszMetadata = tpDriver->GetMetadata ();
            if (CSLFetchBoolean (papszMetadata, GDAL_DCAP_CREATE, FALSE))
            {
                char **options = NULL;
                options = CSLSetNameValue (options, "SPARSE_OK", "YES");
                // set whether bigtiff
                if (bigtiff)
                {
                    options = CSLSetNameValue (options, "BIGTIFF", "YES");
                }
                GDALDataset *gdal = tpDriver->Create (tifTemplateName,
                                                      GRID_SIZE_X, GRID_SIZE_Y,
                                                      1, GDT_Float32, options);
                assert(gdal != NULL);


                // the following geo transform is valid for a "bottom first" row format for tif output
                // double transform[6];
                // transform[0] = min_x;
                // transform[1] = GRID_DIST_X;
                // transform[2] = 0;
                // transform[3] = min_y;
                // transform[4] = 0;
                // transform[5] = GRID_DIST_Y;

                // set the geo transform for "top first" format for tif output
                double transform[6];
                transform[0] = min_x;
                transform[1] = GRID_DIST_X;
                transform[2] = 0;
                transform[3] = min_y + GRID_SIZE_Y;
                transform[4] = 0;
                transform[5] = -GRID_DIST_Y;
                gdal->SetGeoTransform (transform);
                // set the epsg projection wkt
                if (epsg_code)
                {
                    OGRSpatialReference *ogr_sr = new OGRSpatialReference();
                    ogr_sr->importFromEPSG (epsg_code);
                    char *ogr_wkt = NULL;
                    ogr_sr->exportToWkt (&ogr_wkt);
                    gdal->SetProjection (ogr_wkt);
                    dbg(4, "%s\n", ogr_wkt);
                }

                // set raster band and no data value
                GDALRasterBand *tBand = gdal->GetRasterBand (1);
                tBand->SetNoDataValue (-9999.f);
                GDALClose ((GDALDatasetH) gdal);

                // get the size of the template and allocate space for its contents
                struct stat stat_buf;
                stat (tifTemplateName, &stat_buf);
                off_t st_size = stat_buf.st_size;
                unsigned char *b = (unsigned char *) malloc (
                        st_size * (sizeof(unsigned char)));

                // determine tiff_version, read the template and update offsets based tiff_version
                int fd = open (tifTemplateName, O_RDWR);
                lseek (fd, 2, SEEK_SET);
                read (fd, b, 2);
                unsigned long tiff_version = parse_ushort (b);
                dbg(4, "tiff_version = %li", tiff_version);

                unsigned long bytes_per_row;
                unsigned long rows_per_strip;
                unsigned long strip_count;
                unsigned long strip_offsets;
                unsigned long strip_byte_counts;
                unsigned long bytes_per_strip;
                unsigned long bytes_last_strip;

                if (tiff_version == 42) // Regular TIFF
                {
                    // read directory offset
                    lseek (fd, 4, SEEK_SET);
                    read (fd, b, 4);
                    unsigned long directory_offset = parse_uint (b);

                    // read entry count
                    lseek (fd, directory_offset, SEEK_SET);
                    read (fd, b, 2);
                    unsigned long entry_cnt = parse_ushort (b);

                    dbg(4, "directory_offset %li, entry_cnt %li", directory_offset, entry_cnt);

                    for (i = 0; i < entry_cnt; i++)
                    {
                        read (fd, b, 12);
                        unsigned long tag_id = parse_ushort (b);
                        unsigned long data_type = parse_ushort (b + 2);
                        unsigned long data_count = parse_uint (b + 4);
                        unsigned long data_offset = parse_uint (b + 8);
                        dbg(3,
                            "tag_id %li, data_type %li, data_count %li, data_offset %li",
                            tag_id, data_type, data_count, data_offset);

                        if (tag_id == 273) // TIFFTAG_STRIPOFFSETS
                        {
                            strip_count = data_count;
                            strip_offsets = data_offset;
                        }
                        if (tag_id == 279) // TIFFTAG_STRIPBYTECOUNTS
                        {
                            strip_byte_counts = data_offset;
                        }

                        if (tag_id == 278) // TIFFTAG_ROWSPERSTRIP
                        {
                            rows_per_strip = data_offset;
                        }

                    }
                }
                else
                { // tiff_version == 43, BIGTIFF

                    // read directory offset
                    lseek (fd, 8, SEEK_SET);
                    read (fd, b, 8);
                    unsigned long directory_offset = parse_ulong (b);

                    // read entry count
                    lseek (fd, directory_offset, SEEK_SET);
                    read (fd, b, 8);
                    unsigned long entry_cnt = parse_ulong (b);
                    dbg(4, "directory_offset %li, entry_cnt %li", directory_offset, entry_cnt);

                    for (i = 0; i < entry_cnt; i++)
                    {
                        read (fd, b, 20);
                        unsigned long tag_id = parse_ushort (b);
                        unsigned long data_type = parse_ushort (b + 2);
                        unsigned long data_count = parse_ulong (b + 4);
                        unsigned long data_offset = parse_ulong (b + 12);
                        dbg(4,
                            "tag_id %li, data_type %li, data_count %li, data_offset %li",
                            tag_id, data_type, data_count, data_offset);

                        if (data_type == 2)
                        {
                            dbg(4, "data_offset %s", b + 12);
                        }

                        if (tag_id == 273) // TIFFTAG_STRIPOFFSETS
                        {
                            strip_count = data_count;
                            strip_offsets = data_offset;
                        }
                        if (tag_id == 279) // TIFFTAG_STRIPBYTECOUNTS
                        {
                            strip_byte_counts = data_offset;
                        }

                        if (tag_id == 278) // TIFFTAG_ROWSPERSTRIP
                        {
                            rows_per_strip = data_offset;
                        }

                    }
                }


                unsigned long integer_bytes = (tiff_version == 42)?4:8;

                bytes_per_row = GRID_SIZE_X * 4;
                bytes_per_strip = bytes_per_row * rows_per_strip;
                bytes_last_strip = bytes_per_strip;
                if (GRID_SIZE_Y % rows_per_strip)
                {
                    bytes_last_strip = (GRID_SIZE_Y % rows_per_strip)
                            * bytes_per_row;
                }

                //printf("strip offset count and offset %lu %lu\n", strip_count, strip_offsets);
                //printf("strip byte counts offset %lu\n", strip_byte_counts);
                //printf("rows per strip  %lu\n", rows_per_strip);

                // now finally update the strip offsets and strip byte counts
                unsigned long cur_write_offset = st_size;
                for (i = 0; i < strip_count; i++)
                {
                    //write cur_strip
                    lseek (fd, strip_offsets + (integer_bytes * i), SEEK_SET);
                    write (fd, &cur_write_offset, integer_bytes);
                    cur_write_offset += bytes_per_strip;
                    lseek (fd, strip_byte_counts + (integer_bytes * i), SEEK_SET);
                    if (i != strip_count - 1)
                    {
                        write (fd, &bytes_per_strip, integer_bytes);
                    }
                    else
                    {
                        write (fd, &bytes_last_strip, integer_bytes);
                    }
                }
                // seek back to beginning and store the complete template contents for later write by first writer process only
                lseek (fd, 0, SEEK_SET);
                read (fd, b, st_size);
                tifTemplateContents = b;
                tifTemplateCount = st_size;
                close (fd);
            }
        }
        // all writers need the size of the template to add to the starting offset of their contents
        MPI_Bcast (&tifTemplateCount, 1, MPI_UNSIGNED_LONG, first_writer_rank, MPI_COMM_WORLD);

        // allocate and fill the writer's buffer with each writer's set of raster rows
        float *poRasterData = NULL;
        int row_count = w_row_end_index - w_row_start_index + 1;

        if (tifFiles != NULL)
        {
            for (i = 0; i < numTypes; i++)
            {
                if (tifFiles[i] != NULL)
                {
                    //gdalFiles[i] = open_raster (gdalFileNames[i]);
                    if (is_writer)
                    {
                        if (poRasterData == NULL)
                        {
                            poRasterData = new float[row_count * GRID_SIZE_X];
                        }
                        for (j = 0; j < row_count * GRID_SIZE_X; j++)
                        {
                            poRasterData[j] = 0;
                        }
                        // write rows in "top first" format for tif output
                        for (k = row_count-1; k >= 0; k--)
                        {
                            for (j = 0; j < GRID_SIZE_X; j++)
                            {
                                // set index for "top first" format for tif output
                                int index = (row_count-1 - k) * (GRID_SIZE_X) + j;

                                if (interp[k][j].empty == 0
                                        && interp[k][j].filled == 0)
                                {
                                    poRasterData[index] = -9999.f;
                                }
                                else
                                {
                                    switch (i)
                                    {
                                        case 0:
                                            poRasterData[index] =
                                                    interp[k][j].Zmin;
                                            break;

                                        case 1:
                                            poRasterData[index] =
                                                    interp[k][j].Zmax;
                                            break;

                                        case 2:
                                            poRasterData[index] =
                                                    interp[k][j].Zmean;
                                            break;

                                        case 3:
                                            poRasterData[index] =
                                                    interp[k][j].Zidw;
                                            break;

                                        case 4:
                                            poRasterData[index] =
                                                    interp[k][j].count;
                                            break;

                                        case 5:
                                            poRasterData[index] =
                                                    interp[k][j].Zstd;
                                            break;
                                    }
                                }
                            }
                        }
                        // write the template header and all writer's buffers
                        //dbg(5, "tifTemplateName %s, %u, %i, %i, %i\n", tifTemplateName, tifTemplateCount, row_stride, row_count, rank);
                        if (rank == first_writer_rank)
                        {
                            MPI_File_write_at (tifFiles[i], 0,
                                               tifTemplateContents,
                                               tifTemplateCount, MPI_BYTE,
                                               &status);
                        }
                        // The following "bottom first" rows was replaced with the "top first" section below
                        // MPI_Offset offset = tifTemplateCount
                        //        + ((rank - first_writer_rank) * row_stride
                        //                * GRID_SIZE_X * 4);

                        // write the tiff with "top first rows"
                        int last_writer_rank = first_writer_rank + writer_count - 1;
                        MPI_Offset offset = 0;
                        if(rank==last_writer_rank)
                        {
                            offset = tifTemplateCount;
                        }
                        else
                        {
                            // this accounts for the tif being written upside down, and uses the fact that all writer ranks except the last have row_count == row_stride
                            offset = ((unsigned long)(GRID_SIZE_Y - ((rank-first_writer_rank + 1)*row_stride) ))  *   ((unsigned long)GRID_SIZE_X) * 4 + tifTemplateCount;
                        }

                        MPI_File_write_at (tifFiles[i], offset, poRasterData,
                                           (row_count * GRID_SIZE_X * 4),
                                           MPI_BYTE,
                                           &status);

                    }
                }

            }
        }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    if (tifFiles != NULL)
    {
        for (i = 0; i < numTypes; i++)
        {
            if (tifFiles[i] != NULL)
                MPI_File_close (&(tifFiles[i]));
        }
        //unlink (tifTemplateName);
    }



#endif // HAVE_GDAL

    // close files
    MPI_Barrier(MPI_COMM_WORLD);
    if (arcFiles != NULL)
    {
        for (i = 0; i < numTypes; i++)
        {
            if (arcFiles[i] != NULL)
                MPI_File_close (&(arcFiles[i]));
        }
    }

    if (gridFiles != NULL)
    {
        for (i = 0; i < numTypes; i++)
        {
            if (gridFiles[i] != NULL)
                MPI_File_close (&(gridFiles[i]));
        }
    }



    return 0;
}

unsigned long MpiInterp::parse_ulong(unsigned char *buffer)
{
   unsigned long result = 0;
   unsigned long temp = 0;
   temp = buffer[7];
   result |= temp<<56;
   temp = buffer[6];
   result |= temp<<48;
   temp = buffer[5];
   result |= temp<<40;
   temp = buffer[4];
   result |= temp<<32;
   temp = buffer[3];
   result |= temp<<24;
   temp = buffer[2];
   result |= temp<<16;
   temp = buffer[1];
   result |= temp<<8;
   temp = buffer[0];
   result |= temp<<0;

   return result;
}



unsigned int MpiInterp::parse_uint(unsigned char *buffer)
{
   unsigned int result = 0;
   unsigned int temp = 0;
   temp = buffer[3];
   result |= temp<<24;
   temp = buffer[2];
   result |= temp<<16;
   temp = buffer[1];
   result |= temp<<8;
   temp = buffer[0];
   result |= temp<<0;

   return result;
}

unsigned short MpiInterp::parse_ushort(unsigned char *buffer)
{
   unsigned short result = 0;
   unsigned short temp = 0;

   temp = buffer[1];
   result |= temp<<8;
   temp = buffer[0];
   result |= temp<<0;

   return result;
}



