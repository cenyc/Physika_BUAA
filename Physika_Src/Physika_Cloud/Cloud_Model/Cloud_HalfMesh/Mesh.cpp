#pragma once
#include "Mesh.h"
#include"delaunay.h"
#include"edge.h"
#include"numeric.h"
#include"vector2.h"

Mesh::Mesh()
{
	cout << "creat a mesh object" << endl;
}

void Mesh::Initial()
{
	vertexList = NULL;
	edgeList = NULL;
	faceList = NULL;
	ver_number = 0;
	edge_number = 0;
	face_number = 0;

	//------------

}

//怎么写？？？？
void Mesh::CreatBaseMesh()
{
	//从文件读取
	//ifstream fin;
	//fin.open("C:/Users/李云飞/Desktop/boudaryVertexList6.off", std::ofstream::in);
	//if (fin.fail()){
	//	std::cerr << "打开新文件失败！" << std::endl;
	//}
	//else{
	//	string s;
	//	//空过第一行
	//	getline(fin, s);
	//	//读节点数、面片数、边数
	//	getline(fin, s);
	//	size_t pos = s.find(' ');
	//	ver_number = atoi(s.substr(0, pos).c_str());
	//	s = s.substr(pos + 1);
	//	pos = s.find(' ');
	//	face_number = atoi(s.substr(0, pos).c_str());
	//	edge_number = atoi(s.substr(pos + 1).c_str());
	//	vertexList = new float[ver_number * 3];
	//	faceList = new int[face_number * 4];
	//	for (int i = 0; i < ver_number; i++){
	//		float x, y, z;
	//		getline(fin, s);
	//		pos = s.find(' ');
	//		x = atof(s.substr(0, pos).c_str());
	//		s = s.substr(pos + 1);
	//		pos = s.find(' ');
	//		y = atof(s.substr(0, pos).c_str());
	//		z = atof(s.substr(pos + 1).c_str());
	//		vertexList[3 * i + 0] = x;
	//		vertexList[3 * i + 1] = y;
	//		vertexList[3 * i + 2] = z;
	//	}
	//	for (int i = 0; i < face_number; i++){
	//		int num, a, b, c;
	//		getline(fin, s);
	//		pos = s.find(' ');
	//		num = atoi(s.substr(0, pos).c_str());
	//		s = s.substr(pos + 1);
	//		pos = s.find(' ');
	//		a = atoi(s.substr(0, pos).c_str());
	//		s = s.substr(pos + 1);
	//		pos = s.find(' ');
	//		b = atoi(s.substr(0, pos).c_str());
	//		c = atoi(s.substr(pos+1).c_str());
	//		faceList[4 * i + 0] = num;
	//		faceList[4 * i + 1] = a;
	//		faceList[4 * i + 2] = b;
	//		faceList[4 * i + 3] = c;
	//	}
	//}
	//fin.close();

	vector<generate_mesh::Vector2<float>> points;
	ver_number = final_points.size();
	for (int i = 0; i < ver_number; i++){
		float x, y;
		x = final_points[i].x;
		y = final_points[i].y;
		generate_mesh::Vector2<float> position(x, y);
		points.push_back(position);
	}

	//delaunay三角剖分算法生成2维平面mesh网格
	cout << "Generating " << ver_number << " random points" << std::endl;
	Delaunay<float> triangulation;
	const vector<Triangle<float> > triangles = triangulation.triangulate(points);
	cout << triangles.size() << " triangles generated\n";
	const vector<Edge<float> > edges = triangulation.getEdges();

	vertexList = new float[ver_number * 3];

	//传入顶点
	for (int i = 0; i < ver_number; i++){
		vertexList[3 * i + 0] = points[i].x;
		vertexList[3 * i + 1] = points[i].y;
		vertexList[3 * i + 2] = 0.0f;
	}

	//三角形三边找点的索引
	int x = 0;
	vector<Vector3> face_vector;
	for (unsigned i = 0; i < triangles.size(); i++){
		int index1 = 0, index2 = 0, index3 = 0;
		for (int j = 0; j < ver_number; j++){
			generate_mesh::Vector2<float> position(triangles[i].p1.x, triangles[i].p1.y);
			if (points[j] == position){
				index1 = j;
				break;
			}
		}
		for (int j = 0; j < ver_number; j++){
			generate_mesh::Vector2<float> position(triangles[i].p2.x, triangles[i].p2.y);
			if (points[j] == position){
				index2 = j;
				break;
			}
		}
		for (int j = 0; j < ver_number; j++){
			generate_mesh::Vector2<float> position(triangles[i].p3.x, triangles[i].p3.y);
			if (points[j] == position){
				index3 = j;
				break;
			}
		}

		//如果这个面的三个顶点都属于轮廓点，那么舍弃
		if (index1 < edge_point_number && index2 < edge_point_number && index3 < edge_point_number){
			continue;
		}
		else{
			Vector3 face_index;
			face_index.x = index1;
			face_index.y = index3;
			face_index.z = index2;
			//cout << index1 << " " << index3 << " " << index2 << endl;
			face_vector.push_back(face_index);
			x++;
		}

	}

	//存储面
	face_number = face_vector.size();
	faceList = new int[face_vector.size() * 4];
	for (int i = 0; i < face_vector.size(); i++){
		faceList[4 * i + 0] = 3;
		faceList[4 * i + 1] = face_vector[i].x;
		faceList[4 * i + 2] = face_vector[i].y;
		faceList[4 * i + 3] = face_vector[i].z;
		//cout << faceList[4 * i + 0] << " "<<faceList[4 * i + 1] << " " << faceList[4 * i + 2] << " "<<faceList[4 * i + 3] << endl;
	}

}

