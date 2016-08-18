p_points2grid 

Background:

The points2grid application was extended with MPI to allow the application 
to be run in parallel on a compute cluster. The goal was an application that 
would scale to arbitrarily large input, limited only by the amount of disk 
space needed to store the input and output files. No intermediate files are 
generated but the size of the grid output is limited by the sum of the memory
of the writer processes. This code was extended from the 
https://github.com/CRREL/points2grid/ repository which is the most actively
developed fork of the original work by OpenTopography.

The algorithm_and_test_results.pdf slide presentation and 
lidar_parallel_processing.pdf paper contain a description of the implementation
and test results on large LAS files run on the Stampede cluster at the 
University of Texas-Austin. 

Dependencies:

An MPI implementation must be installed. OpenMPI 1.6 and 1.8 are known to work
and were used in development and testing. mpic++ must be found in your PATH.
GDAL 1.9.2 or later with internal BIGTIFF support must be installed.
Boost 1.55 or later must be installed. 

Install:

git clone https://github.com/jwend/p_points2grid
cd p_points2grid
make

The p_points2grid executable is in the bin directory.

Test:

mpirun -n 10 bin/p_points2grid -i test.las -o grid -m parallel --reader_count 4
--mpi_time --buffer_size 4000000 --resolution 1 --fill_window 3

Limitations and Supported Features:

p_points2grid works only with LAS version 1.0, 1.1, and 1.2 input. 
Version 1.2 was tested most extensively with up to the 117 GB file size input,
producing grids up to 27000X27000 cells. This cooresponds to 1 Meter resolution
over a nominal 15 minute X 15 minute geographic extent. 

Input flags specific to p_points2grid are:
-m [--interpolation_mode]
-c [--reader_count]
-b [--buffer_size]
-t [--mpi_times]

The output type of the grid is controlled with --min, --max, --mean, --idw,
--std, -den, or --all. The default is all types.
The oupput format is controlled with --output_format [all|arc|grid].
The default is all formats.
 
These flags are documented with:
bin/p_points2grid --help

For example, here is a run on a cluster on a 117 GB LAS 1.2 file over the 
Grand Canyon area that produces all grids of all types at 1 meter resolution:

mpirun -n 1024 ./p_points2grid -i grand_canyon.las -o gc -m parallel 
--reader_count 256 --mpi_time --resolution 1 --buffer_size 400000 
--fill_window 7

Known Bugs:
--fill_window [3|5|7] can cause the program to fail if the grid row count 
for writer process N is very small and causes overlap with process 
N+2 or N-2. This only happens when writer process counts are very high and 
grid row counts are very low, specifically when 
(writer row count / 2) < fill_window.

    







