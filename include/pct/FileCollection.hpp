#ifndef FILE_COLLECTION_HPP
#define FILE_COLLECTION_HPP
#include "points2grid/Global.hpp"

typedef struct FileCollection
	{
		int count;
		char basePath[1024];
		char fileExt[20];
		char** fileList;
		FileCollection();
		FileCollection(char* _basePath, char* _fileExt);
		~FileCollection();
		int countFiles();
		void clear();
		int getMetadata(int start, int end);

} FileCollection;






#endif // FILE_COLLECTION_HPP
