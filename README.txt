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
The following must also be installed:
GDAL 1.9.2 or later with internal BIGTIFF.
Boost 1.55 or later. 
GEOS 3.4.2 or later.
shapelib 1.3.0 or later.  

Install:

git clone https://github.com/jwend/p_points2grid
cd p_points2grid
make

The p_points2grid executable is in the bin directory.

Test:

mpirun -n 10 bin/p_points2grid -i test.las -o grid -m parallel --reader_count 4
--mpi_time --buffer_size 4000000 --resolution 1 --fill_window 3

Limitations and Supported Features:

p_points2grid works only with LAS version 1.0, 1.1, 1.2, and 1.4 input. 
Version 1.2 was tested most extensively with up to the 117 GB file size input,
producing grids up to 27000X27000 cells. This cooresponds to 1 Meter resolution
over a nominal 15 minute X 15 minute geographic extent. 

Input flags specific to p_points2grid are:
-m [--interpolation_mode]
-c [--reader_count]
-b [--buffer_size]
-t [--mpi_times]
-f [--input_data_file_name ]
--output_bbox
--output_shape
--keep_class
--first_returns
--last_returns
--epsg
--bigtiff
--fill_empty

All flags are documented with:
bin/p_points2grid --help

The output type of the grid is controlled with --min, --max, --mean, --idw,
--std, -den, or --all. The default is all types.
The oupput format is controlled with --output_format [all|arc|grid].
The default is all formats.
 
Examples: 

Here is a run on a cluster on a 117 GB LAS 1.2 file over the 
Grand Canyon area that produces all grids of all types at 1 meter resolution:

mpirun -n 1024 ./p_points2grid -i grand_canyon.las -o gc -m parallel 
--reader_count 256 --mpi_time --resolution 1 --buffer_size 400000 
--fill_window 7

Here is a run on a cluster on 67 LAS 1.4 input files whose paths 
are listed in the file elwha using only points with classification 2, 
producing a 1 meter resolution BIGTIFF with an epsg projection code of 26910. 
Cell values are the mean of point elevation values within the default search 
radius.     

mpirun -n 200 ./p_points2grid -f elwha -o elwha_kc2 -m parallel -c 80 --mean 
--keep_class 2 --output_format tif  -b 1000000 --resolution 1 --epsg 26910  
--mpi_times

Here is a run on a cluster on 67 LAS 1.4 input files whose paths 
are listed in the file elwha using only points with classification 2, 
producing a 1 meter resolution BIGTIFF with an epsg projection code of 26910.
Cell values are the mean of point elevation values within the default search 
radius. The output is clipped to the extent of elwha_urban.shp and empty cells
are filled using row wise linear interpolation.     

mpirun -n 200 ./p_points2grid -f elwha -o elwha_kc2_fill_shape -m parallel 
-c 80 --mean --keep_class 2 --output_format tif  -b 1000000 --resolution 1 
--epsg 26910  --mpi_times --fill_empty --output_shape elwha_urban.shp









