// MeshDeform.cpp : �������̨Ӧ�ó������ڵ㡣
//
//#pragma once
#include "total.h"
#include <Eigen/Sparse>
#include <vector>
#include "MeshDeformation.h"
#include <iostream>
using namespace std;

int main()
{
	MeshDeformation *meshDeformation = new MeshDeformation();

	string basemesh("../basemesh-0001.off");
	meshDeformation->CreateMesh(basemesh);

	cout << "Create Mesh DONE!" << endl;

	cout << "Start Create Optimized Mesh----------" << endl;
	meshDeformation->CreateOptimizedMesh(1);
	cout << "Finish Create Optimized Mesh----------" << endl;
	cout << "Start Create Entire HeightField And Mesh----------" << endl;

	meshDeformation->CreateEntireHeightFieldAndMesh();
	cout << "Finish Create Entire HeightField And Mesh----------" << endl;

	string outmesh("../new-0001.off");
	meshDeformation->OutputDeformedMesh(outmesh);
}

