#ifndef ECEF_HPP
#define ECEF_HPP

#include <math.h>
#include "Point.hpp"
#include <iostream>
#include "geocent.h"
class CRS
{
	private:
		double flattening;
		double inv_flattening;
		double radius;
		double eccentricity;
		double e2;
	public:
		CRS();
		~CRS();
		void transform(struct point *pt);
		void reverse(struct point *pt);
};

#endif
