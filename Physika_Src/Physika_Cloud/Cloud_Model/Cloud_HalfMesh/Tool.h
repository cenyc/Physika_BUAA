#pragma once
#include "Physika_Core\Vectors\vector_3d.h"
#include "global.h"
class Tool {
public:
	float PhaseFunction(Vector3f v1, Vector3f v2);
	float PhaseFunction(float  cosAngle);
	bool isProbabilityGreater(float threahold);
	float Rand(float vaue1, float value2);
	void ComputeTriangleNormal(float normal[3], float PA[3], float PB[3], float PC[3]);
	Vector3f   MatMult(float Mat[9], Vector3f vec);
	void PrintRunnngIfo(char* ifo);
};