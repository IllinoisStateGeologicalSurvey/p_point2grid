#include "pct/Pixel.hpp"
Pixel::Pixel() {
	sum = 0.0;
	count = 0;
	filled = 0;
}

Pixel::~Pixel() {
	sum = 0.0;
	count = 0;
	filled = 0;
}

float Pixel::avg() {
	return sum / count;
}
