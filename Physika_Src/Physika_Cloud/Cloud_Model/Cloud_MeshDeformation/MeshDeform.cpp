// MeshDeform.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "MeshDeformation.h"
using namespace std;

//int main()
//{
//	MeshDeformation *meshDeformation = new MeshDeformation();
//
//	meshDeformation->CreateMesh("./basemesh-0001.off");
//	std::cout << "Create Mesh DONE!" << std::endl;
//
//	std::cout << "Start Create Optimized Mesh----------" << std::endl;
//	meshDeformation->CreateOptimizedMesh(1);
//	std::cout << "Finish Create Optimized Mesh----------" << std::endl;
//	std::cout << "Start Create Entire HeightField And Mesh----------" << std::endl;
//
//	meshDeformation->CreateEntireHeightFieldAndMesh();
//	std::cout << "Finish Create Entire HeightField And Mesh----------" << std::endl;
//
//	meshDeformation->OutputDeformedMesh("./new-0001.off");
//}

void MeshDeform(string half_path,string full_path, int optimize_times)
{
	MeshDeformation *meshDeformation = new MeshDeformation();

	meshDeformation->CreateMesh(half_path);
	std::cout << "Create Mesh DONE!" << std::endl;

	std::cout << "Start Create Optimized Mesh----------" << std::endl;
	meshDeformation->CreateOptimizedMesh(optimize_times);
	std::cout << "Finish Create Optimized Mesh----------" << std::endl;
	std::cout << "Start Create Entire HeightField And Mesh----------" << std::endl;

	meshDeformation->CreateEntireHeightFieldAndMesh();
	std::cout << "Finish Create Entire HeightField And Mesh----------" << std::endl;

	meshDeformation->OutputDeformedMesh(full_path);
	delete meshDeformation;
}
