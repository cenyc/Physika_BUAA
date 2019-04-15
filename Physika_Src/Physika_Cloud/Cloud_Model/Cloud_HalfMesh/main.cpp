#pragma once
#include "global.h"
#include "Image.h"
#include "Pixel.h"
#include "sun.h"
#include "mesh.h"
#include "Cloud.h"

void main() {
	Tool tool;
	Image image;
	Pixel pixel;
	Sun sun;
	Cloud cloud;
	Sky sky;
	//--------
	Mesh mesh;
	//--------
	image.Initial();
	pixel.Initial();
	sun.Initial();
	//--------
	mesh.Initial();
	//--------
	sky.Initial();
	cloud.Initial();
	image.ReadImage("./test1.png", tool);
	pixel.CreatePixelType(image, tool);
	sun.CreateSunColor(pixel, image);
	pixel.CreatePerfectBoundary(image, mesh);
	float a = 0.3*M_PI;
	cout << a << endl;
	sun.CreateSun(a, a);

	//------2019.1.25-----liyunfei-------
	//输入轮廓点和云内部采集点，使用delaunay三角剖分算法，生成二维mesh网格
	mesh.CreatBaseMesh();
	//-----------------------------------

	sky.CreateSkyPossion(image, pixel);
	cloud.CreatePropagationPath(image, sun, pixel);
	cloud.PropagationCylinders(image, sky, sun);
	cloud.CreateHeightFieldHalf(image, pixel, sky);
	cloud.RemoveOneLoopPixelBoudary(image, pixel, cloud.heightField);
	cloud.CreateCloudMesh(image, pixel, sky, mesh, cloud.heightField);

	//------------------
	cloud.ExportCloudMesh(mesh, "./test.off");

	cout << "finished" << endl;

}