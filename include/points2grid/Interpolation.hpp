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

#pragma once

#include <string>
#include <iostream>

using namespace std;

#include <points2grid/export.hpp>
#include <points2grid/GridPoint.hpp>
#include <points2grid/CoreInterp.hpp>
#include <points2grid/OutCoreInterp.hpp>
#include <points2grid/InCoreInterp.hpp>
#include <points2grid/MpiInterp.hpp>
#include <points2grid/export.hpp>
#include <points2grid/Global.hpp>
#include <points2grid/lasfile.hpp>
#include <shapefil.h>
#include <geos/index/quadtree/Quadtree.h>
#include <geos/index/intervalrtree/SortedPackedIntervalRTree.h>
#include <geos/geom/Envelope.h>

#include <geos/algorithm/RayCrossingCounter.h>
#include <geos/algorithm/RobustDeterminant.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Location.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/index/ItemVisitor.h>




class P2G_DLL Interpolation
{
public:
    Interpolation(double x_dist, double y_dist, double radius,
                  int _window_size, int _interpolation_mode, int _rank, int _process_count, int _reader_count, long _buffer_size, mpi_times *_timer);
    ~Interpolation();

    int init(char **inputNames, int inputNamesSize, int inputFormat, int bigtiff, int epsg_code, double *bbox,
             int *classifications, int classification_count, int first_returns, int last_returns, SHPHandle shape_filter);
    int interpolation(char *outputName, int inputFormat, int outputFormat, unsigned int type);
//    unsigned long getDataCount();

    unsigned int getGridSizeX();
    unsigned int getGridSizeY();

    // for debug
    void printArray();

    // deprecated
    void setRadius(double r);

    int point_contained(double x, double y, double rx1, double ry1, double rx2, double ry2);
    int rectangles_overlap (double rx1, double ry1, double rx2, double ry2, double sx1, double sy1, double sx2, double sy2);
    int in_shape(double x, double y);

    int pass_filter(las_file &las, size_t index);

    int countSegment(const geos::geom::Coordinate& p1,
                                     const geos::geom::Coordinate& p2, const geos::geom::Coordinate& point);
    int
    orientationIndex(const geos::geom::Coordinate& p1,
             const geos::geom::Coordinate& p2, const geos::geom::Coordinate& q);

public:
    double GRID_DIST_X;
    double GRID_DIST_Y;

//    static const int MAX_POINT_SIZE = 16000000;
    static const int WEIGHTER = 2;

    // update this to the maximum grid that will fit in memory
    // as a rule of thumb, memory requirement = MEM_LIMIT*55 bytes
    static const unsigned int MEM_LIMIT = 200000000;
    CoreInterp *getInterp(){return interp;}


private:
    double min_x;
    double min_y;
    double max_x;
    double max_y;

    unsigned int GRID_SIZE_X;
    unsigned int GRID_SIZE_Y;

    double radius_sqr;
    int window_size;
    int interpolation_mode;

    int *classifications;
    int classification_count;
    int first_returns;
    int last_returns;

    SHPHandle shape_filter;

    geos::index::intervalrtree::SortedPackedIntervalRTree shape_filter_index;

    struct ShapeSegment
    {
        int shape;
        int vertex_1;
    };

    class ShapeItemVisitor : public geos::index::ItemVisitor
    {
    public:
        std::vector<void *> items;
        void visitItem(void *item)
        {
            items.push_back(item);
        }
    };

    SHPObject **shapes;
    int shape_count;

    int init_shape_filter_index();


    // mpi variables
    int process_count;
    int rank;
    int reader_count;
    long buffer_size;
    input_file_info *input_files;
    int input_file_count;
    mpi_times *timer;
    // end, mpi variables

    CoreInterp *interp;
};

