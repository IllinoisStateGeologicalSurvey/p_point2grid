#ifndef ECEF_HPP
#define ECEF_HPP

#include <math.h>
#include "pct/Point.hpp"
#include <iostream>
#include "pct/geocent.h"
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
		void transform(struct Point *pt);
		void reverse(struct Point *pt);
};

#endif
