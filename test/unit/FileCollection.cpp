#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>
#include <dirent.h>
#include <string.h>
#include <errno.h>
#include "FileCollection.hpp"
#include "util.hpp"
#include <sys/time.h>
#include <sys/resource.h>

int main(int argc, char* argv[]) {
	char path[200] = "/projects/isgs/lidar/champaign/las";
	char fileExt[10] = ".las";
	FileCollection* files = new FileCollection(&path[0], &fileExt[0]);
	printf("FileCount: %i\n", files->count);
	int count = files->count;
	int idx =0;
	files->getMetadata(0,1500);
	printf("Filename at index: %i is %s\n", idx, (files->fileList[idx]));
	//process_mem_usage(vm_usage, resident_set);
	delete(files);
	return 0;
}
