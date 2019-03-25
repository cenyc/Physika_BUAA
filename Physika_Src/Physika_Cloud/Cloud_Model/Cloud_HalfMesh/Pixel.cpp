#include "Pixel.h"
#include "POINT.h"
#include <algorithm>
#include "Mesh.h"
Pixel::Pixel()
{
	cout << "creat an pixel object" << endl;
}

void Pixel::Initial()
{
	skycount = 0;
	cloudPixelCount = 0;
	seg_thresh = 0.479;
	horizontalLine = 280;
	boundaryPixel_number = 0;
	pixelTypeList = NULL;
	boundaryPixelIndexList = NULL;
}
int * Pixel::GetPixelTypeList()
{
	return this -> pixelTypeList;
}

int * Pixel::GetBoundaryList()
{
	return this ->boundaryPixelIndexList;
}

int Pixel::GetBoundaryPixel_num()
{
	return this -> boundaryPixel_number;
}

void Pixel::CreatePixelType(Image image, Tool tool)
{
	tool.PrintRunnngIfo("Create Pixel Type");

	pixelTypeList = new int[image.GetImg_height()*image.GetImg_width()];
	if (pixelTypeList == NULL)
	{
		cout << "allocate memory for pixel type failed!\n";
		exit(1);
	}
	//for the type of pixel, -1---boundary, 1---cloud,0---sky,2--ground,
	for (int i = horizontalLine; i<image.GetImg_height(); i++)
		for (int j = 0; j<image.GetImg_width(); j++)
		{
			pixelTypeList[i*image.GetImg_width() + j] = 2;
		}

	int  height = min(horizontalLine, image.GetImg_height());
	//将图像分割，分离出地面，天空和云
	for (int i = 0; i<height; i++)
		for (int j = 0; j<image.GetImg_width(); j++)
		{

			Color4f color = image.GetImg_mat_cor()[i*image.GetImg_width() + j];

			//---------------------------------------------------------------------------------
			float tempMax = -99999;
			if (tempMax<color.redChannel())
				tempMax = color.redChannel();
			if (tempMax<color.greenChannel())
				tempMax = color.greenChannel();
			if (tempMax<color.blueChannel())
				tempMax = color.blueChannel();

			float tempMin = 99999;
			if (tempMin>color.redChannel())
				tempMin = color.redChannel();
			if (tempMin>color.greenChannel())
				tempMin = color.greenChannel();
			if (tempMin>color.blueChannel())
				tempMin = color.blueChannel();
			//----------------------------------------------------------------------------------

			if ((tempMax - tempMin) / tempMax<seg_thresh)//seg_thresh应该是设置的阈值，分离出云和天空
			{
				pixelTypeList[i*image.GetImg_width() + j] = 1;
				cloudPixelCount++;
			}
			else
			{
				pixelTypeList[i*image.GetImg_width() + j] = 0;
				skycount++;

			}

		}

	bool* CloudBoudaryPixelList = new bool[image.GetImg_height()*image.GetImg_width()];
	if (CloudBoudaryPixelList == NULL)
	{
		cout << "allocate memory for cloud boundary pixel failed!\n";
		exit(1);
	}

	for (int i = 0; i<image.GetImg_height(); i++)
		for (int j = 0; j<image.GetImg_width(); j++)
		{
			if (pixelTypeList[i*image.GetImg_width() + j] == 1)
			{
				if (i - 1>0 && pixelTypeList[(i - 1)*image.GetImg_width() + j] != 1)
				{
					CloudBoudaryPixelList[i*image.GetImg_width() + j] = true;
					continue;
				}
				if (i + 1<image.GetImg_height() && pixelTypeList[(i + 1)*image.GetImg_width() + j] != 1)
				{
					CloudBoudaryPixelList[i*image.GetImg_width() + j] = true;
					continue;
				}
				if (j - 1>0 && pixelTypeList[(i)*image.GetImg_width() + j - 1] != 1)
				{
					CloudBoudaryPixelList[i*image.GetImg_width() + j] = true;
					continue;
				}
				if (j + 1<image.GetImg_width() && pixelTypeList[(i)*image.GetImg_width() + j + 1] != 1)
				{
					CloudBoudaryPixelList[i*image.GetImg_width() + j] = true;
					continue;
				}

				CloudBoudaryPixelList[i*image.GetImg_width() + j] = false;

			}

		}
	//在云的点周围，通过云的四周判断，确定出云的边界

	for (int i = 0; i<image.GetImg_height(); i++)
		for (int j = 0; j<image.GetImg_width(); j++)
		{
			if (pixelTypeList[i*image.GetImg_width() + j] == 1 && CloudBoudaryPixelList[i*image.GetImg_width() + j])
			{
				pixelTypeList[i*image.GetImg_width() + j] = -1;
			}
		}

	delete[] CloudBoudaryPixelList;
	ofstream out("E:/代码/mesh_x/output/pixelifo_test.txt");
	out << 414 << endl;
	for (int i = 0; i<image.GetImg_height(); i++)
		for (int j = 0; j<image.GetImg_width(); j++)
		{
			out << pixelTypeList[i*image.GetImg_width() + j] << " ";

		}
	out << endl;
}

void Pixel::CreatePerfectBoundary(Image image, Mesh &mesh)
{
	if (pixelTypeList == NULL)
		return;


	int* cloudMask = new int[image.GetImg_width()*image.GetImg_height()];

	//---------------
	//内部采点
	int count = 0;
	bool tag = false;
	//---------------

	for (int i = 0; i<image.GetImg_height(); i++)
		for (int j = 0; j<image.GetImg_width(); j++)
		{
			cloudMask[i*image.GetImg_width() + j] = 0;

			if (pixelTypeList[i*image.GetImg_width() + j] == 1) //cloud
			{
				cloudMask[i*image.GetImg_width() + j] = 255;

				tag = true;

				//------------2019.1.25-----liyunfei-------
				//防止图像边界点有白边影响判断，所以除去最外两圈像素
				if (i>2 && i<image.GetImg_height() - 2 && j>2 && j < image.GetImg_width() - 2){
					//采集云内部采样点
					count++;
					//在云内部，每隔五个像素取一次（可更改来调节生成的mesh精细程度）
					if (count % 5 == 0){
						//-------------------------------
						Vector2f point(j*1.0f / image.img_maxWH, (image.GetImg_height() - i)*1.0f / image.img_maxWH);

						//float x, y;
						//y = (image.GetImg_height() - i)*1.0f / image.img_maxWH;
						//x = j*1.0f / image.img_maxWH;
						//Vector2f point;
						//point.x = x;
						//point.y = y;
						mesh.midPoint.push_back(point);
					}
				}
				//------------------------------------------
			}
		}

	IplImage* g_gray = NULL;
	int g_thresh = 120; //设置阈值
	CvMemStorage* g_storage = NULL;//内存存储器是一个可用来存储诸如序列，轮廓，图形,子划分等动态增长数据结构的底层结构。它是由一系列以同等大小的内存块构成，呈列表型 。
	if (g_storage == NULL)
	{
		CvSize  cs;//CvSize，OpenCV的基本数据类型之一。表示矩阵框大小，以像素为精度。与CvPoint结构类似，但数据成员是integer类型的width和height。
		cs.width = image.GetImg_width();
		cs.height = image.GetImg_height();
		g_gray = cvCreateImage(cs, 8, 1);//初始化一张灰度图像
		g_storage = cvCreateMemStorage(0);//初始化内存器大小，为64k
	}
	else
		cvClearMemStorage(g_storage);

	for (int i = 0; i<image.GetImg_height(); i++)
		for (int j = 0; j<image.GetImg_width(); j++)
		{
			cvSetReal2D(g_gray, i, j, cloudMask[i*image.GetImg_width() + j]);

		}
	delete[] cloudMask;

	cvSaveImage("./gray.jpg", g_gray);

	//------------manually modify the image---------------- 手动修改图片
	cout << "If necessary!  modify the gray image!" << endl;

	/*	g_gray=cvLoadImage("gray_14b.jpg",0);*/
	//g_gray = cvLoadImage("gray_test.jpg", 0);

	//上面两行用于手动加载灰度图片
	g_gray = cvLoadImage("./gray.jpg", 0);


	CvSeq* contours = 0;//
	cvThreshold(g_gray, g_gray, g_thresh, 100, CV_THRESH_BINARY);//该函数的典型应用是对灰度图像进行阈值操作得到二值图像。threshold_type=CV_THRESH_BINARY:如果 src(x,y)>threshold ,dst(x,y) = max_value; 否则,dst（x,y）=0;
	cvFindContours(g_gray, g_storage, &contours, sizeof(CvContour), CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE);//检索所有的轮廓，并将其放入list中，压缩水平的、垂直的和斜的部分，也就是，函数只保留他们的终点部分
	cvZero(g_gray);
	if (contours)
	{
		cvDrawContours(g_gray, contours, cvScalarAll(255), cvScalarAll(255), 100);//在图像上绘制外部和内部轮廓，函数cvDrawContours用于在图像上绘制外部和内部轮廓。当thickness >= 0 时，绘制轮廓线；否则填充由轮廓包围的部分，参数100不知道值得是哪个
	}

	cvSaveImage("./initial_boundary.jpg", g_gray);


	vector<vector<POINT>>  ptLists;
	int  maxContourIndex = -1;
	int maxContourCount = -MAXVAL;
	int countourIndex = 0;
	for (; contours != 0; contours = contours->h_next, countourIndex++)
	{
		vector<POINT>  ptList;
		for (int i = 0; i < contours->total; i++)
		{
			CvPoint* cv_pt = (CvPoint*)cvGetSeqElem(contours, i);
			POINT   pt;
			pt.x = cv_pt->x;
			pt.y = cv_pt->y;
			ptList.push_back(pt);

		}

		int cur_count = (ptList.size());
		if (cur_count>maxContourCount)
		{
			maxContourIndex = countourIndex;

		}
		ptLists.push_back(ptList);

	}


	//应该将分离出的轮廓存入一个二维点坐标中
	//delete loop points

	int  pt_count = ptLists[maxContourIndex].size();

	//cout << pt_count << endl;
	int  Search_Number = pt_count / 8;
	bool *   isDeleteList = new bool[pt_count];
	for (int i = 0; i < pt_count; i++)
	{
		isDeleteList[i] = false;

	}

	for (int i = 0; i < pt_count; i++)
	{
		if (isDeleteList[i])
			continue;


		POINT  pt;
		pt.x = (ptLists[maxContourIndex][i]).x;
		pt.y = (ptLists[maxContourIndex][i]).y;

		vector<int> searchIDList;
		for (int j = 1; j <= Search_Number; j++)
		{
			int id = (i + j) % pt_count;

			POINT  cur_pt;
			cur_pt.x = (ptLists[maxContourIndex][id]).x;
			cur_pt.y = (ptLists[maxContourIndex][id]).y;
			if (cur_pt.x == pt.x && cur_pt.y == pt.y)
			{
				for (int k = 0; k<searchIDList.size(); k++)
					isDeleteList[searchIDList[k]] = true;

				isDeleteList[id] = true;



			}
			else
			{
				searchIDList.push_back(id);
			}




		}
		searchIDList.clear();



	}


	boundaryPixel_number = 0;//云边界点的个数
	for (int i = 0; i < pt_count; i++)
	{
		if (!isDeleteList[i])
		{
			boundaryPixel_number++;
		}

	}
	boundaryPixelIndexList = new int[2 * boundaryPixel_number];
	int pixelId = 0;
	for (int i = 0; i < pt_count; i++)
	{
		if (!isDeleteList[i])
		{
			boundaryPixelIndexList[2 * pixelId + 0] = (ptLists[maxContourIndex][i]).x;
			boundaryPixelIndexList[2 * pixelId + 1] = (ptLists[maxContourIndex][i]).y;
			pixelId++;
		}

	}

	//---------------2019.1.25-----liyunfei-------------
	//读取轮廓点
	float *boundaryVertexList = new float[2 * (boundaryPixel_number - 1)];
	for (int i = 0; i<boundaryPixel_number - 1; i++)
	{
		boundaryVertexList[2 * i + 0] = boundaryPixelIndexList[2 * i + 0] * 1.0 / image.img_maxWH;
		boundaryVertexList[2 * i + 1] = (image.img_height - boundaryPixelIndexList[2 * i + 1])*1.0 / image.img_maxWH;
		//cout << image.img_maxWH << " " << endl;
	}

	mesh.edge_point_number = boundaryPixel_number - 1;

	//将轮廓点导入最终点集中
	for (int i = 0; i < boundaryPixel_number - 1; i++){
		//----------------------------------
		Vector2f point(boundaryVertexList[2 * i + 0], boundaryVertexList[2 * i + 1]);

		//Vector2 point;
		//point.x = boundaryVertexList[2 * i + 0];
		//point.y = boundaryVertexList[2 * i + 1];
		mesh.final_points.push_back(point);
	}

	cout << "轮廓点："<<boundaryPixel_number - 1 << endl;

	//将轮廓内部与轮廓边重复点去掉
	for (int i = 0; i < mesh.midPoint.size(); i++){
		int j;
		for (j = 0; j < mesh.final_points.size(); j++){
			if ((mesh.midPoint[i][0] == mesh.final_points[j][0]) && (mesh.midPoint[i][1] == mesh.final_points[j][1])){
				break;
			}
		}
		if (j >= mesh.final_points.size()){
			mesh.final_points.push_back(mesh.midPoint[i]);
		}
	}

	ofstream out("./boudaryVertexList.txt");
	out << mesh.final_points.size() << endl;
	for (int i = 0; i < mesh.final_points.size(); i++){
		out << mesh.final_points[i][0] << "," << mesh.final_points[i][1];
		out << endl;
	}
	//-----------------------------------------------

	//delete[] isDeleteList;
}

void Pixel::SetPixelTypeList(int i, int value)
{
	pixelTypeList[i] = value;
}
