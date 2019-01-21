#ifndef _COMMON_H
#define _COMMON_H

#ifdef WIN32
#include <windows.h>
#include "SDL/SDL.h"

#else
#include <SDL/SDL.h>
#endif

#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

#include <stdlib.h>
#include <stdarg.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "IL/il.h"
#include <math.h>

#define PI 3.141592653589f
#define SQR(x) ( (x) * (x) )

#ifdef _MSC_VER
#pragma warning (disable: 4244)
#pragma warning (disable: 4305)
#endif

//The order of includes is very important!
#include "Vector.h"
#include "Color.h"
#include "Matrix.h"
#include "Camera.h"
#include "GLText.h"
#include "VolumeClouds.h"
#include "Frustum.h"
#include "SkyDome.h"

extern CCamera 		Camera;
extern CFrustum 	Frustum;

inline float /*__fastcall*/ carmack_func(float x)
{
	int carmack;
	float isx, halfx;	//Inverse Squareroot of x

	halfx = 0.5f*x;
	carmack = *(int*)&x; 
	carmack = 0x5f3759df - (carmack>>1); 
	isx = *(float*)&carmack; 

	isx = isx*(1.5f-halfx*isx*isx);  //Newton-Rhapson step, add more for accuracy

	return isx;
}

#endif


