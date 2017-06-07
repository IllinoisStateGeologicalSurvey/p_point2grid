#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include "Point.hpp"
#include "points2grid/lasfile.hpp"
#include "Grid.hpp"
#include "BlockScheme.hpp"
#include "CRS.hpp"
#include "DTypes.h"
#include "mpi.h"

int main(int argc, char* argv[]) {
	char lasTest[60] = "/projects/isgs/dupage/las/POINTS_206.las";
	int mpi_size, mpi_rank, mpi_err;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Info info = MPI_INFO_NULL;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(comm, &mpi_size);
	MPI_Comm_rank(comm, &mpi_rank);
	MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN); /* Return info about errors*/

	//This test should read one or many files and send some point buffers to be written to disk then exit.
	//

	MPI_Finalize();
	return 0;
}
