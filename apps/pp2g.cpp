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
* LIDAR DEM Generation
* This program takes a set of randomly distributed points (ASCII or LAS)
* as input, and generates gridded points (DEM)
* Based on the notes by Prof. Ramon Arrowsmith(ramon.arrowsmith@asu.edu)
* Authors: Han S Kim (hskim@cs.ucsd.edu), Sriram Krishnan (sriram@sdsc.edu)
*
*/

/*

 *
 * FileName: main.cpp
 * Author: Han S Kim (hskim@cs.ucsd.edu), Sriram Krishnan (sriram@sdsc.edu)
 * Description: this program starts from this file
 *
 */
#include <points2grid/config.h>
#include <points2grid/Interpolation.hpp>
#include <points2grid/Global.hpp>

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdexcept>
#include <boost/program_options.hpp>
#include <mpi.h>

#ifdef CURL_FOUND
#include <curl/curl.h>
#endif

namespace po = boost::program_options;

const std::string appName("points2grid");


int main(int argc, char **argv)
{


    int rank = 0;
    int process_count = 1;
    int reader_count = 1;
    int buffer_size = 10000;
    int is_mpi = 0;
    mpi_times *timer = NULL;
    // end, mpi variables
    clock_t t0, t1;

    // parameters

    char **inputNames;
    int inputNamesSize;
    char inputName[1024] = {0};
    char inputURL[2048] = {0};
    char outputName[1024] = {0};

    int input_format = INPUT_LAS;
    int interpolation_mode = INTERP_AUTO;
    int output_format = 0;
    unsigned int type = 0x00000000;
    double GRID_DIST_X = 6.0;
    double GRID_DIST_Y = 6.0;
    double searchRadius = (double) sqrt(2.0) * GRID_DIST_X;
    int window_size = 0;

    // argument processing..
    po::options_description general("General options"),
       df("Data file(s)"),
       ot("Output Type"),
       res("Resolution"),
       nf("Null Filling"),
       desc;

    general.add_options()
    ("help", "produce a help message")
    ("output_file_name,o", po::value<std::string>(), "required. name of output file without extension, i.e. if you want the output file to be test.asc, this parameter shoud be \"test\"")
    ("search_radius,r", po::value<float>(), "specifies the search radius. The default value is square root 2 of horizontal distance in a grid cell")
    ("output_format", po::value<std::string>(), "'all' generates every possible format,\n"
     "'arc' for ArcGIS format,\n"
     "'grid' for Ascii GRID format,\n"
     "'tif' for GDAL GeoTIFF format,\n"
     "the default value is --all")
    ("input_format", po::value<std::string>(), "'ascii' expects input point cloud in ASCII format\n"
     "'las' expects input point cloud in LAS format (default)")
    ("interpolation_mode,m", po::value<std::string>()->default_value("auto"), "'incore' stores working data in memory\n"
     "'outcore' stores working data on the filesystem, 'parallel' uses number of processes specified by mpirun -n\n"
     "'auto' (default) guesses based on the size of the data file")
     ("reader_count,c", po::value<int>(), "when interpolation mode is 'parallel', arg is number of reader processes, default is 1")
     ("buffer_size,b", po::value<int>(), "when interpolation mode is 'parallel', arg is write buffer size in bytes, default is 10000 bytes")
     ("mpi_times,t", "time mpi run");

    df.add_options()
#ifdef CURL_FOUND
    ("data_file_name,i", po::value<std::string>(), "path to unzipped plain text data file")
    ("data_file_url,l", po::value<std::string>(), "URL of unzipped plain text data file"
     "You must specify either a data_file_name or data_file_url.");
#else
//    ("data_file_name,i", po::value<std::string>(), "required. path(s) to unzipped las or plain text data files");
    ("data_file_name,i", po::value< std::vector<string> >()->multitoken(), "required. path(s) to unzipped las or plain text data files");
#endif

    ot.add_options()
    ("min", "the Zmin values are stored")
    ("max", "the Zmax values are stored")
    ("mean", "the Zmean values are stored")
    ("idw", "the Zidw values are stored")
    ("std", "the Zstd values are stored")
    ("den", "the density values are stored")
    ("all", "all the values are stored (default)");

    res.add_options()
    ("resolution", po::value<float>(), "The resolution is set to the specified value. Use square grids.\n"
     "If no resolution options are specified, a 6 unit square grid is used")
    ("resolution-x", po::value<float>(), "The X side of grid cells is set to the specified value")
    ("resolution-y", po::value<float>(), "The Y side of grid cells is set to the specified value");

    nf.add_options()
    ("fill", "fills nulls in the DEM. Default window size is 3.")
    ("fill_window_size", po::value<int>(), "The fill window is set to value. Permissible values are 3, 5 and 7.");

    desc.add(general).add(df).add(ot).add(res).add(nf);

    po::variables_map vm;

    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            cout << "------------------------------------------------------------------------" << endl;
            cout << "   " << appName << " (development version )" << endl;
            cout << "------------------------------------------------------------------------" << endl;
            cout << "Usage: " << appName << " [options]" << endl;
            cout << desc << endl;
            exit(0);
        }

        po::notify(vm);


        if (vm.count("output_format")) {
            std::string of = vm["output_format"].as<std::string>();
            if(of.compare("all") == 0) {
                output_format = OUTPUT_FORMAT_ALL;
            }
            else if(of.compare("arc") == 0)
                output_format = OUTPUT_FORMAT_ARC_ASCII;
            else if(of.compare("grid") == 0)
                output_format = OUTPUT_FORMAT_GRID_ASCII;
            else if(of.compare("tif") == 0)
                            output_format = OUTPUT_FORMAT_GDAL_GTIFF;
            else {
                throw std::logic_error("'" + of + "' is not a recognized output_format");
            }
        }
        // resolution
        if(vm.count("resolution")) {
            float res = vm["resolution"].as<float>();
            GRID_DIST_X = res;
            GRID_DIST_Y = res;
            if(searchRadius == sqrt(2.0) * 6.0)
                searchRadius = (double) sqrt(2.0) * GRID_DIST_X;

            if(GRID_DIST_X == 0) {
                throw std::logic_error("resolution must not be 0");
            }
        }

        if(vm.count("resolution-x")) {
            GRID_DIST_X = vm["resolution-x"].as<float>();
            if(searchRadius == sqrt(2.0) * 6.0)
                searchRadius = (double) sqrt(2.0) * GRID_DIST_X;

            if(GRID_DIST_X == 0) {
                throw std::logic_error("resolution-x must not be 0");
            }
        }

        if(vm.count("resolution-y")) {
            GRID_DIST_Y = vm["resolution-y"].as<float>();
            if(GRID_DIST_Y == 0) {
                throw std::logic_error("resolution-y must not be 0");
            }
        }

        if(vm.count("min")) {
            type |= OUTPUT_TYPE_MIN;
        }

        if(vm.count("max")) {
            type |= OUTPUT_TYPE_MAX;
        }

        if(vm.count("mean")) {
            type |= OUTPUT_TYPE_MEAN;
        }

        if(vm.count("idw")) {
            type |= OUTPUT_TYPE_IDW;
        }

	if(vm.count("std")) {
            type |= OUTPUT_TYPE_STD;
        }

        if(vm.count("den")) {
            type |= OUTPUT_TYPE_DEN;
        }

        if(vm.count("all")) {
            type = OUTPUT_TYPE_ALL;
        }

        if(vm.count("fill")) {
            window_size = 3;
        }

        if(vm.count("fill_window_size")) {
            window_size = vm["fill_window_size"].as<int>();
            if(!((window_size == 3) || (window_size == 5) || (window_size == 7))) {
                throw std::logic_error("window-size must be either 3, 5, or 7");
            }
        }

        if(vm.count("input_format")) {
            std::string inf = vm["input_format"].as<std::string>();

            if(inf.compare("ascii") == 0)
                input_format = INPUT_ASCII;
            else if(inf.compare("las") == 0)
                input_format = INPUT_LAS;
            else {
                throw std::logic_error("'" + inf + "' is not a recognized input_format");
            }
        }


#ifdef CURL_FOUND
        if(vm.count("data_file_name")) {
            strncpy(inputName, vm["data_file_name"].as<std::string>().c_str(), sizeof(inputName));
        }
        if(vm.count("data_file_url")) {
            strncpy(inputURL, vm["data_file_url"].as<std::string>().c_str(), sizeof(inputURL));
        }

        if((inputName == NULL || !strcmp(inputName, "")) &&
                (inputURL == NULL || !strcmp(inputURL, "")))
        {
            throw std::logic_error("you must specify a valid data file");
        }
#else
        if(!vm.count("data_file_name"))
        {
            throw std::logic_error("data_file_name  must be specified");
        }
        else
        {
            inputNames = (char **)malloc(vm["data_file_name"].as<std::vector<string> >().size() * sizeof(char*));
            inputNamesSize = vm["data_file_name"].as<std::vector<string> >().size();

            for(unsigned i = 0; i< vm["data_file_name"].as<std::vector<string> >().size(); i++)
            {
                inputNames[i] = (char *) malloc( ( strlen(vm["data_file_name"].as<std::vector<string> >()[i].c_str()) + 1 ) * sizeof(char));
                strcpy(inputNames[i],vm["data_file_name"].as<std::vector<string> >()[i].c_str());
            }

            strncpy(inputName, inputNames[0], sizeof(inputName));
            if (!strcmp(inputName, ""))
            {
                throw std::logic_error("data_file_name must not be an empty string");
            }
        }
#endif

        if (!vm.count("output_file_name")) {
            throw std::logic_error("output_file_name must be specified");
        }
        else {
            strncpy(outputName, vm["output_file_name"].as<std::string>().c_str(), sizeof(outputName));
        }

        if(vm.count("search_radius")) {
            searchRadius = vm["search_radius"].as<float>();
        }

        if(vm.count("interpolation_mode")) {
            std::string im(vm["interpolation_mode"].as<std::string>());
            if (im.compare("auto") == 0) {
                interpolation_mode = INTERP_AUTO;
            }
            else if (im.compare("incore") == 0) {
                interpolation_mode = INTERP_INCORE;
            }
            else if (im.compare("outcore") == 0) {
                interpolation_mode = INTERP_OUTCORE;
            }
            else if (im.compare("parallel") == 0) {
                 interpolation_mode = INTERP_MPI;
            }
            else {
                throw std::logic_error("'" + im + "' is not a recognized interpolation_mode");
            }
        }

        if(type == 0)
            type = OUTPUT_TYPE_ALL;

        if(vm.count("reader_count")) {
            reader_count = vm["reader_count"].as<int>();
        }
        if(vm.count("buffer_size")) {
            buffer_size = vm["buffer_size"].as<int>();
        }
        if(vm.count("mpi_times")) {
            timer = (mpi_times *)malloc(sizeof(mpi_times));
        }

#ifdef CURL_FOUND
        // download file from URL, and set input name
        if (!((inputURL == NULL || !strcmp(inputURL, "")))) {

            CURL *curl;
            CURLcode res;

            /* get the file name from the URL */
            int i = 0;
            for(i = sizeof(inputURL); i>= 0; i--) {
                if(inputURL[i] == '/')
                    break;
            }
            strncpy(inputName, inputURL+i+1, sizeof(inputName));

            curl = curl_easy_init();
            if (!curl) {
                cout << "Can't initialize curl object to download input from: "
                     << inputURL << endl;
                exit(1);
            }

            /* set URL */
            curl_easy_setopt(curl, CURLOPT_URL, inputURL);
            curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1);

            /* and write to file */
            FILE *fp;
            fp = fopen(inputName, "w");
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);

            /* perform the curl request and clean up */
            res = curl_easy_perform(curl);
            curl_easy_cleanup(curl);
            fclose(fp);

            if (res != 0) {
                cout << "Error while downloading input from: " << inputURL << endl;
                exit(1);
            }
        }
#endif

        if (interpolation_mode == INTERP_MPI)
        {
            if(timer)timer->start = time(NULL);
            MPI_Init (&argc, &argv);
            MPI_Comm_size (MPI_COMM_WORLD, &process_count);
            MPI_Comm_rank (MPI_COMM_WORLD, &rank);
            if (reader_count >= process_count)
            {
                if (rank == 0)
                {
                    printf ("process count = %i is not greater than reader process count = %i.\n",
                            process_count, reader_count);
                    printf ("decrease --parallel(-p) parameter value or add processes.\n");
                }

                MPI_Finalize ();
                return 0;
            }
            MPI_Barrier(MPI_COMM_WORLD);

        }

        if (rank == 0)
        {
            cout << "Parameters ************************" << endl;
            cout << "inputName: '" << inputName << "'" << endl;
            cout << "input_format: " << input_format << endl;
            cout << "outputName: '" << outputName << "'" << endl;
            cout << "GRID_DIST_X: " << GRID_DIST_X << endl;
            cout << "GRID_DIST_Y: " << GRID_DIST_Y << endl;
            cout << "searchRadius: " << searchRadius << endl;
            cout << "output_format: " << output_format << endl;
            cout << "type: " << type << endl;
            cout << "fill window size: " << window_size << endl;
            cout << "parallel interp, reader process count: " << reader_count
                    << endl;
            cout << "parallel interp, writer buffer size: " << buffer_size
                    << endl;
            cout << "************************************" << endl;
        }
        if (interpolation_mode == INTERP_MPI)
        {
            MPI_Barrier (MPI_COMM_WORLD);
        }

    }
    catch (std::exception& e) {
        cerr << "error: " << e.what() << endl;
        cerr << "execute `" << appName << " --help` to see usage information" << endl;
        exit(1);
    }

    t0 = clock();

    if(timer)
    {
        if(rank == reader_count)printf("Allocating memory...\n");
        timer->init_start = time(NULL);
    }
    Interpolation *ip = new Interpolation(GRID_DIST_X, GRID_DIST_Y, searchRadius,
                                          window_size, interpolation_mode, rank, process_count, reader_count, buffer_size, timer);

    if(ip->init(inputNames, inputNamesSize, input_format) < 0)
    {
        fprintf(stderr, "Interpolation::init() error\n");
        return -1;
    }
    if(timer)
    {
        timer->init_end = time(NULL);
    }
    t1 = clock();
    if(rank == 0)
    {
        //printf("Init + Min/Max time: %10.2f\n", (double)(t1 - t0)/CLOCKS_PER_SEC);
    }
    t0 = clock();
    if(ip->interpolation(inputName, outputName, input_format, output_format, type) < 0)
    {
        fprintf(stderr, "Interpolation::interpolation() error\n");
        return -1;
    }

    t1 = clock();
    if (rank == 0 && !timer)
    {
        printf ("DEM generation + Output time: %10.2f\n",
                (double) (t1 - t0) / CLOCKS_PER_SEC);
        printf ("# of data: %d\n", ip->getDataCount ());
        printf ("dimension: %d x %d\n", ip->getGridSizeX (),
                ip->getGridSizeY ());
    }



    if ((interpolation_mode == INTERP_MPI && timer))
    {
        MPI_Barrier (MPI_COMM_WORLD);
        //int first_writer_rank = ip->getInterp ()->getReaderCount ();
        if (timer)
            timer->end = time(NULL);
        /*
         long *interp_start = (long *) malloc(sizeof(long)*process_count);

         MPI_Gather(&(timer->interp_start), 1, MPI_LONG,
         interp_start, 1, MPI_LONG, 0, MPI_COMM_WORLD);
         */
        //int i;
        //if (rank == 0)
        //{
        //for (i = 0; i < process_count; i++)
        //{
        //    printf ("interp_start %li, rank %i\n", interp_start[i], i);
        //}
        //}
        //printf("reader count %i, writer count %i\n",ip->getInterp()->getReaderCount(), ip->getInterp()->getWriterCount());
        //for(i=0; i<process_count; i++){
        //    printf("reader ranks %i writer ranks %i\n",ip->getInterp()->getReaders()[i], ip->getInterp()->getWriters()[i]);
        // }


        if (rank == reader_count)
                {
                    printf("Finished, first writer process times...\n");
                    printf("Total time %li seconds.\n", timer->end - timer->start);

                    printf ("  Allocation: %li seconds.\n",
                            timer->init_end - timer->init_start);
                    printf ("  Read and send: %li seconds.\n",
                            timer->interp_end - timer->interp_start);
                    printf ("  Process cells: %li seconds.\n",
                            timer->process_end - timer->process_start);
                    printf ("  Write cells: %li seconds.\n",
                            timer->output_end -   timer->output_start);

       }


    }

    if (interpolation_mode == INTERP_MPI)
    {
        MPI_Finalize ();
    }
    else
    {
        delete ip;
    }

    return 0;
}
