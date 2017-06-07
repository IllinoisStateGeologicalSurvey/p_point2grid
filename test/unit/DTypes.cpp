#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "DTypes.h"

int main(int argc, char * argv[]) {
	DType datatype = DT_Int16;
	int buf;
	buf = GetDataTypeSize(datatype);
	printf("%s is %i\n", "16-bit Integer", buf);
	return 0;
}
