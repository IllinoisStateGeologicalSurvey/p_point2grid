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

#include <string.h>
#include <math.h>
#include <float.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>

#include <points2grid/lasfile.hpp>
#include <points2grid/debug.hpp>

#include <boost/scoped_ptr.hpp>


#include <fstream>  // std::ifstream
#include <iostream> // std::cerr
#include <string.h>

#include "mpi.h"

/////////////////////////////////////////////////////////////
// Public Methods
/////////////////////////////////////////////////////////////

Interpolation::Interpolation(double x_dist, double y_dist, double radius,
                             int _window_size, int _interpolation_mode = INTERP_AUTO,
                             int _rank = 0, int _process_count = 1, int _reader_count = 1, long _buffer_size = 10000, mpi_times *_timer = NULL) : GRID_DIST_X (x_dist), GRID_DIST_Y(y_dist)
{
    rank = _rank;
    process_count = _process_count;
    reader_count = _reader_count;
    buffer_size = _buffer_size;
    timer = _timer;

    radius_sqr = radius * radius;
    window_size = _window_size;
    interpolation_mode = _interpolation_mode;
    input_files = NULL;
    input_file_count = 0;


    min_x = DBL_MAX;
    min_y = DBL_MAX;

    max_x = -DBL_MAX;
    max_y = -DBL_MAX;

}

Interpolation::~Interpolation()
{
    delete interp;
}

int Interpolation::init(char **inputNames, int inputNamesSize, int inputFormat, int bigtiff, int epsg_code, double *bbox)
{

    clock_t t0, t1;

    t0 = clock();

    char *inputName = inputNames[0];

    if(inputName == NULL)
    {
        cerr << "Wrong Input File Name" << endl;
        return -1;
    }

    if (interpolation_mode != INTERP_MPI)
    {
        // Only the first file is used for non mpi runs, set input_files[0].name is needed during interpolation
        // set it now...
        input_files = (input_file_info *) malloc (sizeof(input_file_info));
        input_file_count = 1;
        input_files[0].name = (char *) malloc ((strlen(inputNames[0]) + 1) * sizeof(char));
        strcpy (input_files[0].name, inputNames[0]);

        if (inputFormat == INPUT_ASCII)
        {
            FILE *fp;
            char line[1024];
            double data_x, data_y;

            if ((fp = fopen (inputName, "r")) == NULL)
            {
                cerr << "file open error" << endl;
                return -1;
            }

            // throw the first line away - it contains the header
            fgets (line, sizeof(line), fp);

            // read the data points to find min and max values
            while (fgets (line, sizeof(line), fp) != NULL)
            {
                data_x = atof (strtok (line, ",\n"));
                if (min_x > data_x)
                    min_x = data_x;
                if (max_x < data_x)
                    max_x = data_x;

                data_y = atof (strtok (NULL, ",\n"));
                if (min_y > data_y)
                    min_y = data_y;
                if (max_y < data_y)
                    max_y = data_y;

            }

            fclose (fp);
        }
        else
        { // las input

            las_file las;
            las.open (inputName);

            min_x = las.minimums ()[0];
            min_y = las.minimums ()[1];
            max_x = las.maximums ()[0];
            max_y = las.maximums ()[1];

            las.close ();

        }
    }
    else if (interpolation_mode == INTERP_MPI)
    {
        if (rank < reader_count)
        {
            input_files = (input_file_info *) malloc (
                    inputNamesSize * sizeof(input_file_info));
            input_file_count = inputNamesSize;
            for (int i = 0; i < input_file_count; i++)
            {
                input_files[i].name = (char *) malloc (
                        (strlen (inputNames[i]) + 1) * sizeof(char));
                strcpy (input_files[i].name, inputNames[i]);

                input_files[i].min_x = 0;
                input_files[i].min_y = 0;
                input_files[i].max_x = 0;
                input_files[i].max_y = 0;
                input_files[i].point_count = 0;
                input_files[i].peek_rank = -1;
            }

            for (int i = 0; i < input_file_count; i++)
            {
                if (rank == i % reader_count)
                {
                    las_file las;
                    las.open (input_files[i].name);
                    input_files[i].min_x = las.minimums ()[0];
                    input_files[i].min_y = las.minimums ()[1];
                    input_files[i].max_x = las.maximums ()[0];
                    input_files[i].max_y = las.maximums ()[1];
                    input_files[i].point_count = las.points_count ();
                    input_files[i].peek_rank = rank;
                    las.close ();
                    //printf("rank = %i i = %i point_count = %li\n", rank, i, input_files[i].point_count);

                }
                //printf("rank = %i i = %i point_count = %li input_files[i].rank = %i\n", rank, i, input_files[i].point_count, input_files[i].peek_rank);
            }
            for (int i = 0; i < input_file_count; i++)
            {
                if (input_files[i].peek_rank != -1)
                {
                    for (int j = 0; j < reader_count; j++)
                    {
                        if (rank != j)
                        {

                            MPI_Send (input_files[i].name, strlen(input_files[i].name)+1, MPI_CHAR, j, 1, MPI_COMM_WORLD);
                            MPI_Send (&(input_files[i].min_x), 1, MPI_DOUBLE, j, 2, MPI_COMM_WORLD);
                            MPI_Send (&(input_files[i].min_y), 1, MPI_DOUBLE, j, 3, MPI_COMM_WORLD);
                            MPI_Send (&(input_files[i].max_x), 1, MPI_DOUBLE, j, 4, MPI_COMM_WORLD);
                            MPI_Send (&(input_files[i].max_y), 1, MPI_DOUBLE, j, 5, MPI_COMM_WORLD);
                            MPI_Send (&(input_files[i].point_count), 1, MPI_LONG, j, 6, MPI_COMM_WORLD);
                            MPI_Send (&(input_files[i].peek_rank), 1, MPI_INT, j, 7, MPI_COMM_WORLD);
                        }
                    }

                }
            }

            for (int i = 0; i < input_file_count; i++)
            {
                if (input_files[i].peek_rank == -1)
                {
                    MPI_Recv (input_files[i].name, strlen(input_files[i].name)+1, MPI_CHAR, i % reader_count, 1, MPI_COMM_WORLD, NULL);
                    MPI_Recv (&(input_files[i].min_x), 1, MPI_DOUBLE, i % reader_count, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv (&(input_files[i].min_y), 1, MPI_DOUBLE, i % reader_count, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv (&(input_files[i].max_x), 1, MPI_DOUBLE, i % reader_count, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv (&(input_files[i].max_y), 1, MPI_DOUBLE, i % reader_count, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv (&(input_files[i].point_count), 1, MPI_LONG, i % reader_count, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv (&(input_files[i].peek_rank), 1, MPI_INT, i % reader_count, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            min_x = min_y = DBL_MAX;
            max_x = max_y = DBL_MIN;

            for (int i = 0; i < input_file_count; i++)
            {
                if (input_files[i].min_x < min_x)
                {
                    min_x = input_files[i].min_x;
                }
                if (input_files[i].min_y < min_y)
                {
                    min_y = input_files[i].min_y;
                }
                if (input_files[i].max_x > max_x)
                {
                    max_x = input_files[i].max_x;
                }
                if (input_files[i].max_y > max_y)
                {
                    max_y = input_files[i].max_y;
                }
            }
        }
    }

    t1 = clock();

    //////////////////////////////////////////////////////////////////////
    // Intialization Step excluding min/max searching
    //////////////////////////////////////////////////////////////////////

    // send input file info collected above by the readers to the writers so all have it
    if(rank == 0){
        for(int i = reader_count; i < process_count; i++){
            MPI_Send (&min_x, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            MPI_Send (&min_y, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            MPI_Send (&max_x, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            MPI_Send (&max_y, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        }
    }
    if(rank>=reader_count){
        MPI_Recv(&min_x, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&min_y, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&max_x, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&max_y, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // assign min max values to bbox, return if there is no overlap
    if (bbox != NULL)
    {
        if (rectangles_overlap (min_x, min_y, max_x, max_y, bbox[0], bbox[1],
                                bbox[2], bbox[3]))
        {
            min_x = bbox[0];
            min_y = bbox[1];
            max_x = bbox[2];
            max_y = bbox[3];
        }
        else
        {
            cerr << "BBOX and LAS input data do not overlap" << endl;
            return -1;
        }

    }

    // all processes set raster size
    GRID_SIZE_X = (int)(ceil((max_x - min_x)/GRID_DIST_X)) + 1;
    GRID_SIZE_Y = (int)(ceil((max_y - min_y)/GRID_DIST_Y)) + 1;

    //dbg(3, ("GRID_SIZE_x %i, grid_size y %i, rank %i\n", GRID_SIZE_X, GRID_SIZE_Y, rank);

    // Construct an Interp based on the mode
    if (interpolation_mode == INTERP_AUTO) {
        // if the size is too big to fit in memory,
        // then construct out-of-core structure
        if(GRID_SIZE_X * GRID_SIZE_Y > MEM_LIMIT) {
            interpolation_mode= INTERP_OUTCORE;
        } else {
            interpolation_mode = INTERP_INCORE;
        }
    }

    if (interpolation_mode == INTERP_OUTCORE) {

        // this is untested and remains for historical reasons
        cerr << "Using out of core interp code" << endl;
        interp = new OutCoreInterp(GRID_DIST_X, GRID_DIST_Y, GRID_SIZE_X, GRID_SIZE_Y, radius_sqr, min_x, max_x, min_y, max_y, window_size);
        if(interp == NULL)
        {
            cerr << "OutCoreInterp construction error" << endl;
            return -1;
        }
        cerr << "Interpolation uses out-of-core algorithm" << endl;

    } else if (interpolation_mode == INTERP_MPI){

        interp = new MpiInterp(GRID_DIST_X, GRID_DIST_Y, GRID_SIZE_X, GRID_SIZE_Y, radius_sqr, min_x, max_x, min_y, max_y, window_size, rank, process_count, reader_count, buffer_size, timer);

    } else {

        interp = new InCoreInterp(GRID_DIST_X, GRID_DIST_Y, GRID_SIZE_X, GRID_SIZE_Y, radius_sqr, min_x, max_x, min_y, max_y, window_size);

    }

    interp->set_bigtiff(bigtiff);
    interp->set_epsg_code(epsg_code);
    if(interp->init() < 0)
    {
        cerr << "inter->init() error" << endl;
        return -1;
    }

    return 0;
}

int Interpolation::interpolation(char *outputName,
                                 int inputFormat,
                                 int outputFormat,
                                 unsigned int outputType)
{
    if(timer)
    {
        if(rank == reader_count)printf("Readers sending points to writers...\n");
        timer->interp_start = time(NULL);
    }

    int rc;

    double data_x, data_y;
    double data_z;
    // This case is for historical preservation only, not used or tested
    if (inputFormat == INPUT_ASCII && interpolation_mode != INTERP_MPI) {
        FILE *fp;
        char line[1024];

        if((fp = fopen(input_files[0].name, "r")) == NULL)
        {
            printf("file open error\n");
            return -1;
        }

        // throw the first line away - it contains the header
        fgets(line, sizeof(line), fp);

        // read every point and generate DEM
        while(fgets(line, sizeof(line), fp) != NULL)
        {
            data_x = atof(strtok(line, ",\n"));
            data_y = atof(strtok(NULL, ",\n"));
            data_z = atof(strtok(NULL, ",\n"));

            data_x -= min_x;
            data_y -= min_y;

            //if((rc = interp->update(arrX[i], arrY[i], arrZ[i])) < 0)
            if((rc = interp->update(data_x, data_y, data_z)) < 0)
            {
                cerr << "interp->update() error while processing " << endl;
                return -1;
            }
        }

        fclose(fp);
    } 

    else
    { // input format is LAS
        if (interpolation_mode == INTERP_MPI)
        {
            if (interp->getIsReader ())
            {
                if (input_file_count == 1) // Single input file
                {
                    if (rectangles_overlap (min_x, min_y, max_x, max_y,
                                            input_files[0].min_x,
                                            input_files[0].min_y,
                                            input_files[0].max_x,
                                            input_files[0].max_y))
                    {
                        las_file las;
                        las.open (input_files[0].name);

                        size_t count = las.points_count ();
                        size_t index (0);
                        size_t stride = count / interp->getReaderCount ();
                        size_t left_over_count = count
                                % interp->getReaderCount ();
                        index = rank * stride;
                        count = (rank + 1) * stride;
                        if (rank == interp->getReaderCount () - 1)
                        {
                            count += left_over_count;
                        }

                        while (index < count)
                        {
                            data_x = las.getX (index);
                            data_y = las.getY (index);
                            data_z = las.getZ (index);
                            if (point_contained (data_x, data_y, min_x, min_y,
                                                 max_x, max_y))
                            {
                                data_x -= min_x;
                                data_y -= min_y;
                                //
                                //cerr << "calling update rank "<< rank << endl;
                                if ((rc = interp->update (data_x, data_y,
                                                          data_z)) < 0)
                                {
                                    cerr
                                            << "interp->update() error while processing "
                                            << endl;
                                    return -1;
                                }
                            }
                            index++;
                        }
                    }
                    interp->getReadDone ()[rank] = 1;

                    for (int i = 0; i < process_count; i++)
                    {
                        if (interp->getWriters ()[i])
                        {
                            interp->updateGridPointSend (i, 1, 1, 1, 1);
                        }
                    }
                }

                else // (input_file_count > 1)
                {

                    for (int i = 0; i < input_file_count; i++)
                    {
                        if (input_files[i].peek_rank == rank)
                        {
                            if (rectangles_overlap (min_x, min_y, max_x, max_y,
                                                    input_files[i].min_x,
                                                    input_files[i].min_y,
                                                    input_files[i].max_x,
                                                    input_files[i].max_y))
                            {
                                dbg(5,
                                    "start read and send rank %i, peek_rank %i, name %s, point_count %li, global and file min_x, miny %lf, %lf, %lf, %lf\n",
                                    rank, input_files[i].peek_rank,
                                    input_files[i].name,
                                    input_files[i].point_count, min_x,
                                    input_files[i].min_x, min_y,
                                    input_files[i].min_y);

                                las_file las;
                                las.open (input_files[i].name);

                                size_t count = input_files[i].point_count;
                                size_t index (0);

                                while (index < count)
                                {
                                    data_x = las.getX (index);
                                    data_y = las.getY (index);
                                    data_z = las.getZ (index);

                                    if (point_contained (data_x, data_y, min_x,
                                                         min_y, max_x, max_y))
                                    {
                                        data_x -= min_x;
                                        data_y -= min_y;
                                        if ((rc = interp->update (data_x,
                                                                  data_y,
                                                                  data_z)) < 0)
                                        {
                                            cerr
                                                    << "interp->update() error while processing "
                                                    << endl;
                                            return -1;
                                        }
                                    }
                                    index++;
                                }
                                las.close ();
                            }
                            dbg(2, "*** rank %i, input file %s done with read and send ***", rank, input_files[i].name);
                        }
                    }

                    interp->getReadDone ()[rank] = 1;

                    for (int i = 0; i < process_count; i++)
                    {
                        if (interp->getWriters ()[i])
                        {
                            interp->updateGridPointSend (i, 1, 1, 1, 1);
                        }
                    }

                }
            }
            else if (interp->getIsWriter ())            // rank is a writer
            {
                interp->updateGridPointRecv ();
            }
            //MPI_Barrier (MPI_COMM_WORLD);
        }
        else // interpolation_mode is other than INTERP_MPI,
             // single input file only supported, bbox supported

        {
            if (rectangles_overlap (min_x, min_y, max_x, max_y,
                                    input_files[0].min_x, input_files[0].min_y,
                                    input_files[0].max_x, input_files[0].max_y))
            {

                las_file las;
                las.open (input_files[0].name);

                size_t count = las.points_count ();
                size_t index (0);
                while (index < count)
                {
                    data_x = las.getX (index);
                    data_y = las.getY (index);
                    data_z = las.getZ (index);
                    if (point_contained (data_x, data_y, min_x, min_y, max_x,
                                         max_y))
                    {

                        data_x -= min_x;
                        data_y -= min_y;

                        if ((rc = interp->update (data_x, data_y, data_z)) < 0)
                        {
                            cerr << "interp->update() error while processing "
                                    << endl;
                            return -1;
                        }
                    }
                    index++;
                }
            }
        }
    }
    if(timer)timer->interp_end = time(NULL);

    if((rc = interp->finish(outputName, outputFormat, outputType)) < 0)
    {
        cerr << "interp->finish() error" << endl;
        return -1;
    }


    return 0;
}

void Interpolation::setRadius(double r)
{
    radius_sqr = r * r;
}

unsigned int Interpolation::getGridSizeX()
{
    return GRID_SIZE_X;
}

unsigned int Interpolation::getGridSizeY()
{
    return GRID_SIZE_Y;
}

// rx1, ry1 = lower left, rx2, ry2 = upper right
int Interpolation::point_contained(double x, double y, double rx1, double ry1, double rx2, double ry2)
{
    return (x >= rx1 && x <= rx2 && y >= ry1 && y <= ry2);
}
// rx1, ry1 = lower left, rx2, ry2 = upper right of first rectangle, sx1, sy1 = lower left, sx2, sy2 = upper right of second rectangle
int Interpolation::rectangles_overlap (double rx1, double ry1, double rx2, double ry2, double sx1, double sy1, double sx2, double sy2)
{
    return sx1 < rx2 && sx2 > rx1 && sy1 < ry2 && sy2 > ry1;
}



