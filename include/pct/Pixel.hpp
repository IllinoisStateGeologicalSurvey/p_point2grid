#ifndef PIXEL_HPP
#define PIXEL_HPP

typedef struct Pixel {
	int count;
	float sum;
	char filled;
	Pixel();
	~Pixel();
	float avg();

} Pixel;

#endif
