#pragma once
#include "global.h"

struct Box//初始化规则体表示的分辨率为350*350*350的规则体
{
	float  x_min;
	float  y_min;
	float  z_min;
	float  x_max;
	float  y_max;
	float  z_max;

	Box()
	{
		x_min = 0;
		y_min = 0;
		z_min = -1;
		x_max = 1;
		y_max = 1;
		z_max = 1;
	}

	void Set(float x_min, float x_max, float y_min, float y_max, float  z_min, float z_max)
	{
		this->x_min = x_min;
		this->y_min = y_min;
		this->z_min = z_min;
		this->x_max = x_max;
		this->y_max = y_max;
		this->z_max = z_max;
	}
	bool IsIn(Vector3f pt)
	{
		if (pt[0] >= x_min&&pt[0] <= x_max&&pt[1] >= y_min&&pt[1] <= y_max&&pt[2] >= z_min&&pt[2] <= z_max)
			return  true;
		else
			return false;

	}

};
class CloudInterior
{
public:
	CloudInterior(void);
	Box  CloudBox;

	float interval_x;
	float interval_y;
	float interval_z;

	int* isInCloud;
	int* isInCloudTemp;
	void Update(const Cylinder& curCylinder);
	virtual float PathLen(Vector3f P, Vector3f direction, const Cylinder& extraLocalVolume);
	bool FindSegment(float intesection[2], Vector3f start, Vector3f end);

	float* heightField;
	void  CreatHeightField();
	float Interpolat(float x, float y);
private:
	bool isInLocalVolume(Vector3f p0, const Cylinder&  localVolume);
	bool IsCloudCube(int x, int y, int z);

public:
	~CloudInterior(void);
};