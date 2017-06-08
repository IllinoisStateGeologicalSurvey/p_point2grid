/** Dtypes: This is a basic struct to allow for parameterization of data types in grids.
 * This will follow the GDALDataType class definition **/

#ifndef DTYPES_H
#define DTYPES_H

typedef enum {
	DT_Unknown = 0,
	DT_Byte = 1,
	DT_UInt16 = 2,
	DT_Int16 = 3,
	DT_UInt32 = 4,
	DT_Int32 = 5,
	DT_Float32 = 6,
	DT_Float64 = 7,
	DT_CInt16 = 8,
	DT_CInt32 = 9,
	DT_CFloat32 = 10,
	DT_CFloat64 = 11,
	DT_TypeCount = 12 /* Maximum type count */
} DType;

/** Returns the size of the data type in bytes (like the sizeof() operator **/
int GetDataTypeSize(DType datatype);

const char* GetDataTypeName(DType datatype);

#endif
