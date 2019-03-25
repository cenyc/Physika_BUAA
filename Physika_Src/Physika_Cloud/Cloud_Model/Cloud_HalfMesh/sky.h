#pragma once
#include <Physika_Render\Color\color.h>
#include "Image.h"
#include "Pixel.h"
using namespace Physika;

class Sky {
public:
	Color4f* img_sky;
	float*  img_grey_sky;
	Sky();
	void Initial();
	void CreateSkyPossion(Image image, Pixel pixel);
	float* GetImg_grey_sky();
};