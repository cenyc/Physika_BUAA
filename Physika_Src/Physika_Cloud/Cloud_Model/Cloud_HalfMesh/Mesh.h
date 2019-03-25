#pragma once
#include "Pixel.h"
#include "Image.h"


class Mesh {
public:
	float* vertexList;
	int* edgeList;
	int* faceList;
	int ver_number;
	int edge_number;
	int face_number;

	//------2019.1.25-----liyunfei-------
	//因为边界点存储在final_points向量的最前端，所以该变量存储了边界点的数目，即最后一个边界点在向量中的位置（用于最后扣除轮廓以外的面）
	int edge_point_number;
	//存取云内部采样点
	vector<Vector2f> midPoint;
	//存取轮廓点+云内部采样点
	vector<Vector2f> final_points;
	//-----------------------------------

	Mesh();
	void Initial();
	void CreatBaseMesh();

    //*****************
	//for  image mesh
	float *  Cloud_vertexList;
	int*     Cloud_facelist;
	int      Cloud_vertexnumber;
	int      Cloud_facenumber;
};