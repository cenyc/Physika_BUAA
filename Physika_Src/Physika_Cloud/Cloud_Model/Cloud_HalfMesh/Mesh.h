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
	//��Ϊ�߽��洢��final_points��������ǰ�ˣ����Ըñ����洢�˱߽�����Ŀ�������һ���߽���������е�λ�ã��������۳�����������棩
	int edge_point_number;
	//��ȡ���ڲ�������
	vector<Vector2f> midPoint;
	//��ȡ������+���ڲ�������
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