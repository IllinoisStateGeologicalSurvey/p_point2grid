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
    data_count = 0;
    radius_sqr = radius * radius;
    window_size = _window_size;
    interpolation_mode = _interpolation_mode;
    input_files = NULL;
    input_file_count = 0;


    min_x = DBL_MAX;
    min_y = DBL_MAX;

    max_x = -DBL_MAX;
    max_y = -DBL_MAX;
    //printf ("mode %i\n", interpolation_mode);
}

Interpolation::~Interpolation()
{
    delete interp;
}

int Interpolation::init(char **inputNames, int inputNamesSize, int inputFormat)
{


    //unsigned int i;

    //struct tms tbuf;
    clock_t t0, t1;


    //////////////////////////////////////////////////////////////////////
    // MIN/MAX SEARCHING
    // TODO: the size of data, min, max of each coordinate
    //       are required to implement streaming processing....
    //
    // This code can be eliminated if database can provide these values
    //////////////////////////////////////////////////////////////////////

    //t0 = times(&tbuf);
    t0 = clock();

    char *inputName = inputNames[0];

    if(inputName == NULL)
    {
        cerr << "Wrong Input File Name" << endl;
        return -1;
    }

    //printf("inputName: '%s'\n", inputName);
    if (interpolation_mode != INTERP_MPI)
    {
        if (inputFormat == INPUT_ASCII)
        {
            FILE *fp;
            char line[1024];
            double data_x, data_y;
            //double data_z;

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

                data_count++;

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

            data_count = las.points_count ();

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
                    //printf("rank = %i i = %i data_count = %li\n", rank, i, input_files[i].point_count);

                }
                //printf("rank = %i i = %i data_count = %li input_files[i].rank = %i\n", rank, i, input_files[i].point_count, input_files[i].peek_rank);
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
            data_count = 0;
            for (int i = 0; i < input_file_count; i++){
                if(input_files[i].min_x < min_x){ min_x = input_files[i].min_x; }
                if(input_files[i].min_y < min_y){ min_y = input_files[i].min_y; }
                if(input_files[i].max_x > max_x){ max_x = input_files[i].max_x; }
                if(input_files[i].max_y > max_y){ max_y = input_files[i].max_y; }
                data_count += input_files[i].point_count;
                //printf("rank = %i data_count = %li\n", rank, input_files[i].point_count);
            }

        }

    }

    t1 = clock();

    //////////////////////////////////////////////////////////////////////
    // Intialization Step excluding min/max searching
    //////////////////////////////////////////////////////////////////////

    GRID_SIZE_X = (int)(ceil((max_x - min_x)/GRID_DIST_X)) + 1;
    GRID_SIZE_Y = (int)(ceil((max_y - min_y)/GRID_DIST_Y)) + 1;

    if(rank == 0){
        for(int i = reader_count; i < process_count; i++){
            MPI_Send(&GRID_SIZE_X, 1, MPI_UNSIGNED, i, 1, MPI_COMM_WORLD);
            MPI_Send(&GRID_SIZE_Y, 1, MPI_UNSIGNED, i, 1, MPI_COMM_WORLD);
            MPI_Send (&min_x, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            MPI_Send (&min_y, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            MPI_Send (&max_x, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            MPI_Send (&max_y, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        }
    }
    if(rank>=reader_count){
        MPI_Recv(&GRID_SIZE_X, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&GRID_SIZE_Y, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&min_x, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&min_y, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&max_x, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&max_y, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //printf("GRID_SIZE_x %i, grid_size y %i, rank %i\n", GRID_SIZE_X, GRID_SIZE_Y, rank);

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
        cerr << "Using out of core interp code" << endl;;

        interp = new OutCoreInterp(GRID_DIST_X, GRID_DIST_Y, GRID_SIZE_X, GRID_SIZE_Y, radius_sqr, min_x, max_x, min_y, max_y, window_size);
        if(interp == NULL)
        {
            cerr << "OutCoreInterp construction error" << endl;
            return -1;
        }

        cerr << "Interpolation uses out-of-core algorithm" << endl;

    } else if (interpolation_mode == INTERP_MPI){
        // printf("before constructor rank %i  pc %i rc %i \n", rank, process_count, reader_count);

        interp = new MpiInterp(GRID_DIST_X, GRID_DIST_Y, GRID_SIZE_X, GRID_SIZE_Y, radius_sqr, min_x, max_x, min_y, max_y, window_size, rank, process_count, reader_count, buffer_size, timer);


    } else {
        cerr << "Using incore interp code" << endl;

        interp = new InCoreInterp(GRID_DIST_X, GRID_DIST_Y, GRID_SIZE_X, GRID_SIZE_Y, radius_sqr, min_x, max_x, min_y, max_y, window_size);

        cerr << "Interpolation uses in-core algorithm" << endl;
    }
    // printf("before init rank %i\n", rank);
    if(interp->init() < 0)
    {
        cerr << "inter->init() error" << endl;
        return -1;
    }
    // printf("after init rank %i\n", rank);
    //cerr << "Interpolation::init() done successfully" << endl;

    return 0;
}

int Interpolation::interpolation(char *inputName,
                                 char *outputName,
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

    if(rank == 0)
    {
        //printf("Interpolation Starts, rank %i\n", rank);
    }

    if (inputFormat == INPUT_ASCII && interpolation_mode != INTERP_MPI) {
        FILE *fp;
        char line[1024];

        if((fp = fopen(inputName, "r")) == NULL)
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
                if (input_file_count == 1)
                {
                    printf("input file count == 1\n");
                    strcpy(inputName, input_files[0].name);
                    las_file las;
                    las.open (inputName);

                    size_t count = las.points_count ();
                    size_t index (0);
                    size_t stride = count / interp->getReaderCount ();
                    size_t left_over_count = count % interp->getReaderCount ();
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

                        data_x -= min_x;
                        data_y -= min_y;
                        //
                        //cerr << "calling update rank "<< rank << endl;
                        if ((rc = interp->update (data_x, data_y, data_z)) < 0)
                        {
                            cerr << "interp->update() error while processing "
                                    << endl;
                            return -1;
                        }
                        index++;
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

                            dbg(5, "start read and send rank %i, peek_rank %i, name %s, point_count %li, global and file min_x, miny %lf, %lf, %lf, %lf\n", rank, input_files[i].peek_rank, input_files[i].name, input_files[i].point_count, min_x, input_files[i].min_x, min_y, input_files[i].min_y);

                            las_file las;
                            las.open (input_files[i].name);

                            size_t count = input_files[i].point_count;
                            size_t index (0);

                            while (index < count)
                            {
                                data_x = las.getX (index);
                                data_y = las.getY (index);
                                data_z = las.getZ (index);
                                assert(data_x >= input_files[i].min_x && data_x <= input_files[i].max_x &&
                                       data_y >= input_files[i].min_y && data_y <= input_files[i].max_y);
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
                                index++;
                            }
                            las.close();
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
        else // interpolation_mode is other than INTERP_MPI
        {

            las_file las;
            las.open(inputName);


            size_t count = las.points_count();
            size_t index(0);
            while (index < count) {
                data_x = las.getX(index);
                data_y = las.getY(index);
                data_z = las.getZ(index);

                data_x -= min_x;
                data_y -= min_y;

                if ((rc = interp->update(data_x, data_y, data_z)) < 0) {
                    cerr << "interp->update() error while processing " << endl;
                    return -1;
                }
                index++;
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

unsigned long Interpolation::getDataCount()
{
    return data_count;
}

unsigned int Interpolation::getGridSizeX()
{
    return GRID_SIZE_X;
}

unsigned int Interpolation::getGridSizeY()
{
    return GRID_SIZE_Y;
}

