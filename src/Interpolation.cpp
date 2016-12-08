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
#include "signal.h"

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

    shapes = NULL;
    shape_count = 0;

}

Interpolation::~Interpolation()
{
    delete interp;
}

int Interpolation::init(char **inputNames, int inputNamesSize, int inputFormat, int bigtiff, int epsg_code, double *bbox,
                        int *_classifications, int _classification_count, int _first_returns, int _last_returns, SHPHandle _shape_filter)
{

    clock_t t0, t1;

    t0 = clock();

    char *inputName = inputNames[0];
    classifications = _classifications;
    classification_count = _classification_count;
    first_returns = _first_returns;
    last_returns = _last_returns;
    shape_filter = _shape_filter;



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
                        size_t points_processed(0);
                        while (index < count)
                        {
                            data_x = las.getX (index);
                            data_y = las.getY (index);
                            data_z = las.getZ (index);
                            if (point_contained (data_x, data_y, min_x, min_y,
                                                 max_x, max_y))
                            {
                                if (pass_filter(las, index))
                                {
                                    points_processed ++;
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
                            }
                            index++;
                        }
                        dbg(5, "rank = %i, index = %li, points_processed = %li", rank, index, points_processed);
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
                                size_t points_processed (0);
                                while (index < count)
                                {
                                    data_x = las.getX (index);
                                    data_y = las.getY (index);
                                    data_z = las.getZ (index);

                                    if (point_contained (data_x, data_y, min_x,
                                                         min_y, max_x, max_y))
                                    {
                                        if (pass_filter (las, index))
                                        {
                                            points_processed++;
                                            data_x -= min_x;
                                            data_y -= min_y;
                                            if ((rc = interp->update (data_x,
                                                                      data_y,
                                                                      data_z))
                                                    < 0)
                                            {
                                                cerr
                                                        << "interp->update() error while processing "
                                                        << endl;
                                                return -1;
                                            }
                                        }
                                    }
                                    index++;
                                }
                                dbg(3, "rank = %i, index = %li, points_processed = %li", rank, index, points_processed);
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
             // single input file only supported, bbox supported, filters supported

        {
            if (rectangles_overlap (min_x, min_y, max_x, max_y,
                                    input_files[0].min_x, input_files[0].min_y,
                                    input_files[0].max_x, input_files[0].max_y))
            {

                las_file las;
                las.open (input_files[0].name);

                size_t count = las.points_count ();
                size_t index (0);
                size_t points_processed(0);
                while (index < count)
                {
                    data_x = las.getX (index);
                    data_y = las.getY (index);
                    data_z = las.getZ (index);
                    if (point_contained (data_x, data_y, min_x, min_y, max_x,
                                         max_y))
                    {
                        if (pass_filter (las, index))
                        {
                            points_processed++;
                            data_x -= min_x;
                            data_y -= min_y;

                            if ((rc = interp->update (data_x, data_y, data_z))
                                    < 0)
                            {
                                cerr
                                        << "interp->update() error while processing "
                                        << endl;
                                return -1;
                            }
                        }
                    }
                    index++;
                }
                dbg(3, "index = %li, points_processed = %li", index, points_processed);
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

int Interpolation::pass_filter (las_file &las, size_t index)
{


    if(!in_shape(las.getX(index), las.getY(index)))
    {
        return 0;
    }

    if(classification_count == 0 && first_returns == 0 && last_returns == 0)
    {
        return 1;
    }
    // At this point there either is no shape_filter shape_filter has been passed, check the point attr filters
    int point_classification = las.getClassification(index);
    for(int i=0; i < classification_count; i++)
    {
        if(point_classification == classifications[i])
        {
            return 1;
        }
    }
    if (first_returns && las.getReturnNumber(index)==1)
    {
        return 1;
    }
    if (last_returns && las.getReturnNumber(index) == las.getReturnCount(index))
    {
        return 1;
    }

    return 0;

}

int Interpolation::in_shape(double x, double y)
{
    //printf("in_shape\n");
    if(shape_filter == NULL)
    {
        return 1;
    }
    //printf("in_shape, before init\n");
    if(shapes == NULL)
    {
        init_shape_filter_index();
    }
    int cross_count = 0;

    geos::geom::Coordinate las_point;
    las_point.x = x;
    las_point.y = y;

    geos::geom::Envelope env(x-1, max_x + 1, y - 1, y + 1);

    std::vector<void*> found_segments;

    shape_filter_index.query(&env, found_segments);


    geos::algorithm::RayCrossingCounter rcc(las_point);


    for (std::size_t i = 1; i < found_segments.size (); i++)
    {

        ShapeSegment *segment = (ShapeSegment *)(found_segments[i]);
        double x1 =
                shapes[segment->shape]->padfX[segment->vertex_1];
        double y1 =
                shapes[segment->shape]->padfY[segment->vertex_1];
        double x2 =
                shapes[segment->shape]->padfX[segment->vertex_1 + 1];
        double y2 =
                shapes[segment->shape]->padfY[segment->vertex_1 + 1];

        geos::geom::Coordinate p1 (x1, y1);
        geos::geom::Coordinate p2 (x2, y2);
        printf("%lf %lf   %lf   %lf %lf   %lf\n", x1,x2, x,  y1, y2, y );
        cross_count = countSegment (p1, p2, las_point);

        //    if ( rcc.isOnSegment() )
        //            return rcc.getLocation();
    }


    return cross_count % 2;


    /*
    double min_y = y - max_short_segment_length/2;
    double max_y = y + max_short_segment_length/2;


    std::vector<FilterIndex>::iterator begin;
    std::vector<FilterIndex>::iterator end;

    FilterIndex lower;
    lower.Y = min_y;
    lower.index = 0;
    FilterIndex upper;
    upper.Y = max_y;
    upper.index = 0;

    begin = std::lower_bound(shape_filter_short_segments.begin(), shape_filter_short_segments.end(), lower, FilterIndexLess);
    end = std::upper_bound(shape_filter_short_segments.begin(), shape_filter_short_segments.end(), upper, FilterIndexLess);
*/

}




int Interpolation::init_shape_filter_index()
{
    shape_count = shape_filter->nRecords;
    shapes = (SHPObject **) malloc(shape_count * sizeof(SHPObject *));

    for(int i=0; i<shape_count; i++)
    {
        shapes[i] =  SHPReadObject(shape_filter, i);
        int part_cnt = shapes[i]->nParts;
        for(int j=0; j<part_cnt; j++)
        {
            int start_vertex = shapes[i]->panPartStart[j];
            int end_vertex;
            if(j < part_cnt - 1)
            {
                end_vertex = shapes[i]->panPartStart[j+1] - 1;
            }
            else
            {
                end_vertex =  shapes[i]->nVertices - 1;
            }
            for(int k = start_vertex; k < end_vertex; k++)
                // end_vertex is same as first vertex in a ring, so segment count is end_vertex - start_vertex - 1
            {
                double x1,x2,y1,y2;
                x1=shapes[i]->padfX[k];
                x2=shapes[i]->padfX[k+1];
                y1=shapes[i]->padfY[k];
                y2=shapes[i]->padfY[k+1];
                //if(fabs(x1-x2)<1.0) x2+=1.0;
                //else if(fabs(y1-y2)<1.0) y2+=1.0;

                geos::geom::Envelope  env(x1,x2,y1,y2);

                ShapeSegment *segment = (ShapeSegment *)malloc(sizeof(ShapeSegment));
                segment->shape = i;
                segment->vertex_1 = k;
                shape_filter_index.insert(&env, segment);
                printf("filter insert %i\n", k);
            }

        }
    }
    printf("%s\n", shape_filter_index.toString().c_str());


    return 1;
}




int
Interpolation::countSegment(const geos::geom::Coordinate& p1,
                                 const geos::geom::Coordinate& p2, const geos::geom::Coordinate& point)
{


        int cross = 0;
        // For each segment, check if it crosses
        // a horizontal ray running from the test point in
        // the positive x direction.

        // check if the segment is strictly to the left of the test point
        if (p1.x < point.x && p2.x < point.x)
                return 0;

        // check if the point is equal to the current ring vertex
        if (point.x == p2.x && point.y == p2.y)
        {
                //isPointOnSegment = true;
                return 0;
        }

        // For horizontal segments, check if the point is on the segment.
        // Otherwise, horizontal segments are not counted.
        if (p1.y == point.y && p2.y == point.y)
        {
                //double minx = p1.x;
                //double maxx = p2.x;

                //if (minx > maxx)
                //{
                //        minx = p2.x;
                //        maxx = p1.x;
                //}

                //if (point.x >= minx && point.x <= maxx)
                        //isPointOnSegment = true;

                return 0;
        }

        // Evaluate all non-horizontal segments which cross a horizontal ray
        // to the right of the test pt.
        // To avoid double-counting shared vertices, we use the convention that
        // - an upward edge includes its starting endpoint, and excludes its
        //   final endpoint
        // - a downward edge excludes its starting endpoint, and includes its
        //   final endpoint
        if (((p1.y > point.y) && (p2.y <= point.y)) ||
                ((p2.y > point.y) && (p1.y <= point.y)) )
        {
                // For an upward edge, orientationIndex will be positive when p1->p2
                // crosses ray. Conversely, downward edges should have negative sign.







                int sign = orientationIndex(p1, p2, point);
                if (sign == 0)
                {
                        //isPointOnSegment = true;
                        return 0;
                }

                if (p2.y < p1.y)
                        sign = -sign;

                // The segment crosses the ray if the sign is strictly positive.
                if (sign > 0)
                        cross++;
        }
        return cross;
}




int
Interpolation::orientationIndex(const geos::geom::Coordinate& p1,
         const geos::geom::Coordinate& p2, const geos::geom::Coordinate& q)
{
        // travelling along p1->p2, turn counter clockwise to get to q return 1,
        // travelling along p1->p2, turn clockwise to get to q return -1,
        // p1, p2 and q are colinear return 0.
        double dx1=p2.x-p1.x;
        double dy1=p2.y-p1.y;
        double dx2=q.x-p2.x;
        double dy2=q.y-p2.y;
        return geos::algorithm::RobustDeterminant::signOfDet2x2(dx1,dy1,dx2,dy2);
}





/*
double Interpolation::get_standard_deviation(double a[], int size)
{
    double sum = 0.0, mean;
    int i;
    for(i=0; i<size; i++)
    {
        sum += a[i];
    }
    mean = sum/size;
    sum = 0;
    for(i=0; i<size; i++)
    {
        sum += pow(a[i] - mean, 2);
    }
    return sqrt(sum/size);
}
*/


/*
    int shape_filter_segment_count = shape_filter_object->nVertices -1; // point is same as first, don't record it
    double *segment_lengths = (double *) malloc(shape_filter_segment_count * sizeof(double));


    for (int i = 0; i < shape_filter_segment_count; i++)
    {
        segment_lengths[i] = fabs(shape_filter_object->padfY[i] - shape_filter_object->padfY[i + 1]);
    }
    double segment_length_std_deviation = get_standard_deviation(segment_lengths, shape_filter_segment_count);
    max_short_segment_length = 4 * segment_length_std_deviation;
    int short_segments = 0;
    int long_segments = 0;
    int *segment_is_short  = (int *) malloc(shape_filter_segment_count * sizeof(int));
    for (int i = 0; i < shape_filter_segment_count; i++)
    {

        if(segment_lengths[i] <= max_short_segment_length)
        {
           segment_is_short[i] = 1;
           short_segments++;
        }
        else
        {
            segment_is_short[i] = 0;
            long_segments++;
        }
    }
    dbg(3, "before malloc %i %i %lf", short_segments, long_segments, max_short_segment_length);

    FilterIndex *shape_filter_short_indices = (FilterIndex *) malloc(short_segments * sizeof(FilterIndex));
    FilterIndex *shape_filter_long_indices = (FilterIndex *) malloc(long_segments * sizeof(FilterIndex));
    int short_index = 0;
    int long_index = 0;
    for(int i=0; i<shape_filter_segment_count; i++)
    {
        if (segment_is_short)
        {
            shape_filter_short_indices[short_index].Y = shape_filter_object->padfY[i];
            shape_filter_short_indices[short_index].index = i;
            short_index++;
        }
        else
        {
            shape_filter_long_indices[long_index].Y = shape_filter_object->padfY[i];
            shape_filter_long_indices[long_index].index = i;
            long_index++;
        }
    }



    if(short_segments > 0)
    {
        shape_filter_short_segments = std::vector<FilterIndex>(shape_filter_short_indices, shape_filter_short_indices + short_segments);
    }
    if(long_segments >0)
    {
        shape_filter_long_segments = std::vector<FilterIndex>(shape_filter_long_indices, shape_filter_long_indices + long_segments);
    }

    for (unsigned i = 0; i < shape_filter_short_segments.size (); i++)
    {
        std::cout << shape_filter_short_segments[i].Y << " " << shape_filter_short_segments[i].index << " ";
    }
    std::cout << '\n';

    std::sort (shape_filter_short_segments.begin (), shape_filter_short_segments.end (), FilterIndexLess);

    for (unsigned i = 0; i < shape_filter_short_segments.size (); i++)
    {
        std::cout << shape_filter_short_segments[i].Y << " " << shape_filter_short_segments[i].index << " ";
    }
    std::cout << '\n';
*/


