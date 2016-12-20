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

#include <iostream>
#include <points2grid/GridPoint.hpp>
#include <points2grid/CoreInterp.hpp>
#include <points2grid/GridFile.hpp>
#include <points2grid/Global.hpp>


#include <mpi.h>

using namespace std;

class P2G_DLL MpiInterp : public CoreInterp
{
public:
    MpiInterp() {};
    MpiInterp(double dist_x, double dist_y,
                 int size_x, int size_y,
                 double r_sqr,
                 double _min_x, double _max_x,
                 double _min_y, double _max_y,
                 int _window_size,
                 int _rank, int _process_count, int _reader_count, long _buffer_size,
                 mpi_times *_timer);
    ~MpiInterp();

    virtual int init();
    virtual int update(double data_x, double data_y, double data_z);
    virtual int finish(char *outputName, int outputFormat, unsigned int outputType);
    virtual int finish(char *outputName, int outputFormat, unsigned int outputType, double *adfGeoTransform, const char* wkt);
    virtual void updateGridPointSend(int target_rank, int x, int y, double data_z, double distance);
    virtual void updateGridPointRecv();
    virtual int getIsReader(){return is_reader;};
    virtual int getIsWriter(){return is_writer;};
    virtual int* getReaders(){return readers;};
    virtual int* getWriters(){return writers;};
    virtual int* getReadDone(){ return read_done; };
    virtual int  getReaderCount(){ return reader_count; };
    virtual int  getWriterCount(){ return writer_count; };



private:
    int rank;
    int process_count;
    long buffer_size;
    mpi_times *timer;
    MPI_Datatype mpi_grid_point_info;
    MPI_Datatype mpi_grid_point;
    GridPoint **interp;
    double radius_sqr;
    int is_fill_empty_cells;
    void fill_empty_cells();
    void fill_empty_cell(int row, int col, int begin, int reset, int last_call);

private:
    GridPoint ** allocRows(int cnt);
    void update_first_quadrant(double data_z, int base_x, int base_y, double x, double y);
    void update_second_quadrant(double data_z, int base_x, int base_y, double x, double y);
    void update_third_quadrant(double data_z, int base_x, int base_y, double x, double y);
    void update_fourth_quadrant(double data_z, int base_x, int base_y, double x, double y);

    void updateGridPoint(int x, int y, double data_z, double distance);
    void printArray();
    int outputFile(char *outputName, int outputFormat, unsigned int outputType, double *adfGeoTransform, const char* wkt);
    void flushMpiBuffers (MPI_File *arcFiles, MPI_File *gridFiles, int numTypes, long max_buffer_size);
    unsigned long parse_ulong(unsigned char *buffer);
    unsigned int parse_uint(unsigned char *buffer);
    unsigned short parse_ushort(unsigned char *buffer);


    int row_stride;  // set for all processes, this is the y grid stride between workers, last worker my have smaller row size
    int w_row_start_index; // set only for workers aka writers
    int w_row_end_index; // set only for workers aka writers
    int get_target_rank(int grid_index);
    int reader_count;
    int writer_count;
    int is_reader;
    int is_writer;
    int *readers; //  array process_count long that holds 0 or 1 flags indicating whether a process is a reader, rank order
    int *writers; //  array process_count long that holds 0 or 1 flags indicating whether a process is a writer, rank order
    int *read_done; //array process_count long that holds 0 or 1 flags indicating whether a process is done reading points,
                    // used by readers to indicate they are done sending points to writers

    // each write process uses these to recv points from a reader
    grid_point_info *point_buffer;
    long point_buffer_count;
    // each read porocess uses these to send points to a specific writer
    grid_point_info **point_buffers;
    long *point_buffer_counts;
    long mpi_point_buffer_count; // holds the number of points to be allocated in the point buffer
                                 // based on number of writers and assumed memory availability of 500 Meg/core
                                 // calculated now in init, tbd, make configurable...


    // writing
    long mpi_buffer_size;

    MPI_Offset *arc_file_mpi_offset;
    char** arc_file_mpi_buffer;   // each element is a buffer for mpi writes of mpi_buffer_size
    long *arc_file_mpi_count;      // each element is the current data size of cooresponding element of arc_file_mpi_buffers
    long *arc_file_mpi_size;        // each element is the total size of data that will be written by a worker
                                            // used to set write offsets calculated by first sprintf loop
    long **arc_file_mpi_sizes;      // used in an Allgather to calculate the write offset of THIS process


    MPI_Offset *grid_file_mpi_offset;
    char** grid_file_mpi_buffer;   // each element is a buffer for mpi writes of mpi_buffer_size
    long *grid_file_mpi_count;      // each element is the current data size of cooresponding element of arc_file_mpi_buffers
    long *grid_file_mpi_size;        // each element is the total size of data that will be written by a worker
                                            // used to set write offsets calculated by first sprintf loop
    long **grid_file_mpi_sizes;      // used in an Allgather to calculate the write offset of THIS process

    MPI_Offset *tif_file_mpi_offset;
    long *tif_file_mpi_size;        // each element is the total size of data that will be written by a worker
    long **tif_file_mpi_sizes;      // used in an Allgather to calculate the write offset of THIS process


};

