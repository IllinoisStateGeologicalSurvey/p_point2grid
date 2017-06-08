#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>
#include <dirent.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <limits.h>
#include "points2grid/Global.hpp"
#include "FileCollection.hpp"
#include "util.hpp"

FileCollection::FileCollection(char* _basePath, char* _fileExt) {
	DIR *dir;
	struct dirent *entry;
	struct stat statbuf;
	strncpy(basePath, _basePath, strlen(_basePath));
	strncpy(fileExt, _fileExt, strlen(_fileExt));
	//printf("BasePath: %s\n", basePath);
	//printf("File ext: %s\n", fileExt);
	count = countFiles();
	fileList = (char**)malloc(sizeof(char*) * count); 
}

int FileCollection::countFiles() {
	DIR *dir;
	struct dirent *entry;
	char * ext;

	int i = 0;
	if ((dir = opendir(basePath)) == NULL)
	{
		fprintf(stderr, "Error: Failed to open base directory\n");
		exit(1);
	}
	chdir(basePath);
	while ((entry = readdir(dir)) != NULL)
	{
		if ((ext = strrchr(entry->d_name, '.')) != NULL)
		{
			if (strcmp(ext, fileExt) == 0)
			{
				i++;
			}
		}
	}
	return i;
}

void FileCollection::clear() {
	int i = 0;
	for (i = 0; i < count; i++) {
		if (fileList[i] != NULL) {
			free(fileList[i]);
			fileList[i] = NULL;
		}
	}
};

int FileCollection::getMetadata(int start, int end) {
	DIR *dir;
	struct dirent *entry;
	struct stat statbuf;
	double vm_usage, resident_set = 0.0;
	int counter = 0;
	char * ext;
	char fname[2056] = "";
	if ((dir = opendir(basePath)) == NULL)
	{
		fprintf(stderr, "IO Error: failed to open base directory: %s, ERRCODE: %s\n", basePath, strerror(errno));
		exit(1);
	}
	chdir(basePath);
	while(((entry = readdir(dir)) != NULL) && (counter < end))
	{
		lstat(entry->d_name, &statbuf);
		if (S_ISDIR(statbuf.st_mode)) {
			/* Found a directory, but ignore links to parent directory */
			if (strcmp(".", entry->d_name) == 0 || strcmp("..", entry->d_name) == 0)
				continue;
			//TODO: Handle directory recursion
		}
		else 
		{
			if (S_ISREG(statbuf.st_mode)) 
			{
				if ((ext = strrchr(entry->d_name, '.')) != NULL) 
				{
					if (strcmp(ext, fileExt) == 0)
					{
						//If the counter is within the read threshold, 
						//allocate the filename and add to the collection
						if (counter >= start) {
							//printf("Appending file index %i for %s\n", counter, entry->d_name);
							memset(fname, 0, sizeof(fname));
							int pathLen = sprintf(&fname[0], "%s/%s",basePath, entry->d_name);
							//strcat(fname, basePath);
							//strcat(fname, "/");
							//strcat(fname, entry->d_name);
							//strcat(fname, '\0');
							//printf("Found file: %s\n", fname);
							// Allocate the file path
							//int pathLen = strlen(fname);
							//int pathLen = strlen(basePath) + strlen(entry->d_name) + 2;
							fileList[counter] = (char*) malloc(sizeof(char) * pathLen);
							strcpy(fileList[counter], fname);
							//printf("Copied path %i: %s\n", counter, fileList[counter]);
						}
						counter++;
						//process_mem_usage(vm_usage, resident_set);
						//printf("Using %d, %d memory\n", vm_usage/1024, resident_set/1024);
					}
					else 
					{ 
						continue;
					}
				}
			}
		}
	}
	chdir("..");
	closedir(dir);
	return counter;
};



FileCollection::~FileCollection() {
	count = 0;
	//Empty string
	memset(basePath,0,strlen(basePath));
	memset(fileExt,0,strlen(fileExt));
	int i = 0;
	for (i = 0; i < count; i++) {
		if (fileList[i] != NULL) {
			free(fileList[i]);
			fileList[i] = NULL;
		}
	}
	free(fileList);
}
