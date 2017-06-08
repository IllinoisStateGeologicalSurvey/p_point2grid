#ifndef UTIL_HPP
#define UTIL_HPP



double randomVal(int range, int min);

void process_mem_usage(double& vm_usage, double& resident_set);

double to_degrees(double radians);

void compareMin(double* mins, double* tmp);

void compareMax(double* maxs, double* tmp);
#endif
