#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "DTypes.h"

int GetDataTypeSize(DType datatype) {
	int buf;
	switch(datatype) {
		case(DT_Unknown): buf = 0;
						break;
		case(DT_Byte): buf = 1;
						break;
		case(DT_UInt16): buf = 2;
						break;
		case(DT_Int16): buf = 2;
						break;
		case(DT_CInt16): buf=2;
						break;
		case(DT_UInt32): buf = 4;
						break;
		case(DT_Int32): buf = 4;
						break;
		case(DT_CInt32): buf = 4;
						break;
		case(DT_Float32): buf = 4;
						break;
		case(DT_CFloat32): buf = 4;
						break;
		case(DT_Float64): buf = 8;
						break;
		case(DT_CFloat64): buf = 8;
						break;
		default: buf = 0;
				 break;
	}
	return buf;
}

const char* GetDataTypeName(DType datatype) {
	switch(datatype) {
		case(DT_Byte): return "Byte";
		case(DT_UInt16): return "Unsigned 16-bit Integer";
		case(DT_Int16): return "Signed 16-bit Integer";
		case(DT_CInt16): return "Complex 16-bit Integer";
		case(DT_UInt32): return "Unsigned 32-bit Integer";
		case(DT_Int32): return "Signed 32-bit Integer";
		case(DT_CInt32): return "Complex 32-bit Integer";
		case(DT_Float32): return "32-bit Floating Point";
		case(DT_CFloat32): return "Complex 32-bit Floating Point";
		case(DT_Float64): return "64-bit Floating Point";
		case(DT_CFloat64): return "Complex 64-bit Floating Point";
		default: return "Unknown Datatype";
	}
}

