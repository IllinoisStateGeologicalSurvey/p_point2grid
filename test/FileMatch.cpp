#include <glob.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>

std::string matchToEnv(boost::smatch const& match) {
	const char * s = getenv( match[1].str().c_str() );
	return s == NULL ? std::string("") : std::string(s);
}


std::string expandEnvironmentVariables(const std::string& in) {
	static boost::regex re( "\\$\\{([^}]+)\\}" );
	return boost::regex_replace(in, re, matchToEnv);
}


/**
 * \brief divide_tasks: Subset a series of tasks given a number of processes 
 *
 * \param task_count (int): Number of tasks to subset
 * \param proc_count (int): Number of processes available
 * \param offsets (int *): Array of integers to hold the subset offsets
 * \param blocks (int *): Array of integers holding the size of each subset
 * 
 */
int divide_tasks(int task_count, int proc_count, int* offsets, int* blocks) {
	int i = 0;
	int remainder = 0;
	int blockSize = 0;
	int blockRemainder = 0;
	blockSize = floor(task_count / proc_count);
	remainder = task_count % proc_count;
	if (remainder > 0) {
		blockSize++;
		for (i = 0; i < proc_count; i++) {
			blockRemainder = i - remainder;
			if (blockRemainder >= 0) {
				offsets[i] = (i * blockSize) - blockRemainder;
				blocks[i] = blockSize - 1;
			} else {
				offsets[i] = i * blockSize;
				blocks[i] = blockSize;
			}
			printf("[%i] Task offset: %i, size: %i, remainder: %i\n", i, offsets[i], blocks[i], blockRemainder);
		}
	} else {
		for (i = 0; i < proc_count; i++) {
			offsets[i] = i * blockSize;
			blocks[i] = blockSize;
		}
	}
	return 0;
}


int main(int argc, char **argv) {
	namespace po = boost::program_options;
	
	typedef std::vector<std::string> filelist_t;
	//Container to hold filename list
	filelist_t iPuts;
	int rank;
	int process_count;
	
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Info info = MPI_INFO_NULL;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &process_count);
	MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	
	int t_offsets[process_count];
	int t_blocks[process_count];

	// - Set up
	po::options_description desc("Allowed Options");
	desc.add_options()
		("help,h", "produce help message")
		// Multitoken means -i foo bar baz works
		("input-file,i", po::value<filelist_t>(&iPuts)->multitoken(), "input files");

	po::positional_options_description pDesc;
	pDesc.add("input-file", -1);
	
	// - Parse
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(pDesc).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 0;
	}

	std::cout << "User provided inputs : " << iPuts.size() << std::endl;
	std::cout << "Raw:" << std::endl;
	std::copy(iPuts.begin(),
			iPuts.end(),
			std::ostream_iterator<filelist_t::value_type>(std::cout, "\n"));
	
	// handle expansion of regex expressions, envvars ?
	

	filelist_t expandEnvVars;
	std::transform(iPuts.begin(), iPuts.end(), std::back_inserter(expandEnvVars), expandEnvironmentVariables);
	// Now actual regex?
	// - If we're on UNIX, glob may be an option
	filelist_t::iterator start = expandEnvVars.begin();
	filelist_t::iterator end = expandEnvVars.end();
	
	filelist_t globED;
	int i;
	// Parse globbed path for files 
	while (start != end) {
		glob_t globBuffer;
		glob((*start).c_str(), GLOB_TILDE, NULL, &globBuffer);
		std::cout << "[" << rank << "] " << *start << " total: " << globBuffer.gl_pathc << std::endl;
		// Divide file paths across processess
		divide_tasks(globBuffer.gl_pathc, process_count, &t_offsets[0], &t_blocks[0]);
		// Print the outpust subset definitions
		if (rank == 0)
		{
			for (i = 0; i < process_count; i++) {
				std::cout << "[" << i << "] Block Definition: " << t_offsets[i] << ", " << t_blocks[i] << std::endl;
			}
		}
		// Copy filenames to subsets
		for (i = t_offsets[rank]; i < t_offsets[rank] + t_blocks[rank]; ++i) {
			globED.push_back(globBuffer.gl_pathv[i]);
		}
		
		// Free the glob buffer
		globfree(&globBuffer);
		++start;


	}
	//Generate files containing the per-process subset
	std::ofstream pathfile;
	std::ostringstream outFileName;
	outFileName <<  "las_" << rank << ".txt";
	std::string outName = outFileName.str();
	pathfile.open(outName.c_str());
	// Read paths out to file
	filelist_t::iterator path_start = globED.begin();
	filelist_t::iterator path_end = globED.end();
	while (path_start != path_end) {
		pathfile << "[" << rank << "] "<< *path_start << std::endl;
		++path_start;
	}
	pathfile.close();
	MPI_Finalize();


	//std::copy(globED.begin(),
	//		globED.end(),
	//		std::ostream_iterator<filelist_t::value_type>(std::cout,"\n"));

	return 0;
}
