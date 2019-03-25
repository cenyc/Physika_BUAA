#pragma once
#include "global.h"
#include <Physika_Render\Color\color.h>
#include "Tool.h"
using namespace Physika;

class Image {
public:
	IplImage* pImg;
	int img_width;
	int img_height;
	int img_maxWH;
	float *   img_mat;
	Color4f* img_mat_cor;
	Image();
	void Initial(void);
	void ReadImage(char* filename, Tool tool);
	int GetImg_width();
	int GetImg_height();
	int GetImg_maxWH();
	Color4f* GetImg_mat_cor();
	float* GetImg_mat();

};