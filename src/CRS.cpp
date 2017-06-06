#include <math.h>
#include <iostream>
#include "Point.hpp"
#include "CRS.hpp"
#include "geocent.h"
/** 
 * Methods taken from https://www.osti.gov/scitech/servlets/purl/110235
 */
/*******************************************
 *			DEFINES
 */
//#define HALF_PI			(M_PI / 2.0e0)
//#define AD_C			1.002600		/* Ralph Toms region 1 constant */


/**********************************************
 *			Global declarations
 */
/* Ellipsoidal parameters, default to WGS 84 */
//double Geocent_a = 6378137.0; /* Semi-major axis of ellipsoid in meters */
//double Geocent_b = 6356752.3142; /* Semi-minor axis of ellipsoid */

//double Geocent_a2 = 40680631590769.0;	/* Square of semi-major axis */
//double Geocent_b2 = 40408299984087.05;	/* Square of semi-minor axis */
//double Geocent_e2 = 0.0066943799901413800; /* Eccentricity squared */
//double Geocent_ep2 = 0.00673949675658690300; /* 2nd eccentricity squared */
/* 
 * These state variables are for optimization purposes. The only function
 * that should modify them is Set_Geocentric_Parameters.
 */


CRS::CRS() {
	radius = 6378137.0;
	inv_flattening = 298.257224;
	flattening = 1/inv_flattening;
	e2 = 2 * flattening - pow(flattening, 2.0);
	eccentricity = sqrt(e2);
}

CRS::~CRS() {
	radius = 0.0;
	inv_flattening = 0.0;
	flattening = 0.0;
	e2 = 0.0;
	eccentricity = 0.0;
}

void CRS::transform(struct point *pt) {
	//double C;
	//double S;
	
	double lon = pt->x * M_PI / 180.0;
	double lat = pt->y * M_PI / 180.0;
	double height = pt->z;
	//std::cout << "e = " << eccentricity << ", C: " << C << ", S: " << S << "\n" << std::flush;
	//C = 1 / (sqrt(pow(cos(lat), 2.0) + (pow((1 - flattening), 2.0) * pow(sin(lat), 2.0))));
	//S = pow((1 - flattening), 2.0) * C;
	Convert_Geodetic_To_Geocentric(lat, lon, height, &pt->x, &pt->y, &pt->z);

	//pt->x = (radius * C + alt) * cos(lat) * cos(lon);
	//pt->y = (radius * C + alt) * cos(lat) * sin(lon);
	//pt->z = (radius * S + alt) * sin(lat);

}


// as defined in https::/www.osti.gov/scitech/servlets/purl/110235/
void CRS::reverse(struct point *pt) {
	double x = pt->x;
	double y = pt->y;
	double z = pt->z;
	
	Convert_Geocentric_To_Geodetic(x, y, z, &pt->y, &pt->x, &pt->z);
	// Convert to degrees
	
	pt->x = pt->x * 180 / M_PI;
	pt->y = pt->y * 180 / M_PI;
	
}



	
