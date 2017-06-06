#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <unistd.h>
#include <ios>
#include <fstream>
#include <string>
#include <math.h>

double randomVal(int range, int min) {
	int newMin = min;
	bool isNeg = false;
	double randDbl = 0.0;
	if (min < 0) {
		newMin = abs(min);
		isNeg = true;
	} 
	if (isNeg) {
		randDbl = (double) (rand() % range - newMin);
	} else {
		randDbl = (double) (rand() % range + newMin);
	}
	return randDbl;
};


///////////////////////////////////////////////////////////////////////////
// 
// process_mem_usage( double &, double &) - takes two doubles by reference
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB
//
// On failure returns 0.0 0.0

void process_mem_usage(double& vm_usage, double& resident_set)
{
	using std::ios_base;
	using std::ifstream;
	using std::string;

	vm_usage = 0.0;
	resident_set = 0.0;

	// 'file' stat seems to give the most reliable results
	ifstream stat_stream("/proc/self/stat", ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	//
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	
	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
		>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
		>> utime >> stime >> cutime >> cstime >> priority >> nice 
		>> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about these

	stat_stream.close();

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage = vsize / 1024.0;
	resident_set = rss * page_size_kb;
}

double to_degrees(double radians) {
	return radians * (180.0 / M_PI);
}


	
