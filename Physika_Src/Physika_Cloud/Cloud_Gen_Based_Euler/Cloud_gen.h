#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<fstream>

//已经默认网格长度为1
#define GRID_LENGTH 1
//判断该网格是否存在物质的阈值
#define F_MIN 1e-6

//------------demo的初始条件(默认构造函数)-----------------
//初始化环境最大/最小温度
#define MAX_AMBIENT_TEMP 25.0
#define MIN_AMBIENT_TEMP 5.0
//初始化地面水气密度
#define INIT_VAPOR_DENSITY 1.0f
//初始化地面水气温度
#define INIT_MATERIAL_TEMP 15.0
//初始化地面水气上升初速度
#define INIT_VELOCITY_W 0.3
//初始化地面水气范围
#define BEGIN_X 28
#define END_X 37
#define BEGIN_Y 28
#define END_Y 37
//---------------------------------------------------------

class CloudGen
{
public:
	//生成模式(0 - 烟, 1 - 云)
	CloudGen(int mode_flag);
	//配置文件路径(二进制文件输入,包含：n_、source_change_times_、Source)
	CloudGen(std::string input_filename, int mode_flag);
	~CloudGen();

	//直接生成对应次数的最终结果(输入运行次数)
	void CloudGenFinalRun(int times);
	//用于demo展示每帧运行()
	void CloudGenFrameRun();

	//---------用于demo显示-----------
	const float GetVelocityU(int i, int j, int k);
	const float GetVelocityV(int i, int j, int k);
	const float GetVelocityW(int i, int j, int k);
	const float GetVaporDensity(int i, int j, int k);
	const float GetCloudDensity(int i, int j, int k);
	const float GetAmbientTemp(int k);
	const float GetMaterialTemp(int i, int j, int k);
	const int GetFrameCount();

private:
	int n_;  //模拟空间大小(不包含边界)
	int size_;
	float dt_;   //时间步长
	int mode_flag;   //生成模式(0-烟,1-云)
	int frame_count_;   //用于记录运行到第几帧

	int source_change_times_;  //用于记录源变化的帧数
	struct Source
	{
		float vapor_density;    //源(地面)水气密度
		float material_temp;    //源(地面)水气温度
		float velocity_w;       //源(地面)水气初始速度
	};
	std::vector<std::vector<Source>>  source_;   //帧数*模拟空间地面（x*y）大小*每个网格的源的情况

	std::vector<float> velocity_u_;   //x方向速度
	std::vector<float> velocity_v_;   //y方向速度
	std::vector<float> velocity_w_;   //z方向速度
	std::vector<float> vapor_density_;  //水气密度
	std::vector<float> cloud_density_;  //云密度
	std::vector<float> ambient_temp_;  //环境温度
	std::vector<float> material_temp_;  //物质温度

	//初始化条件文件输入(其中包含:)
	void InitFileInput(std::string &filename);
	//源(地面)二维坐标位置转化一维向量函数
	const int SourcePosition(int i, int j);
	//空间三维坐标位置转化一维向量函数
	const int Position(int i, int j, int k);
	//计算速度场
	void GetVelocity();
	//计算密度场
	void GetDensity();
	//边界控制(value:控制的对象名  velocity_project_flag:速度投影标记(0-不是，1-x轴，2-y轴，3-z轴)
	void BoundaryCondition(std::vector<float> &value, const int velocity_project_flag);
	//Advect平流(要平流的量，该量上一帧的值，所需要的三个方向的速度)
	void Advect(std::vector<float> &value, const std::vector<float> &pre_value, const std::vector<float> &velocity_u, const std::vector<float> &velocity_v, const std::vector<float> &velocity_w);
	//高斯赛德尔迭代
	void GaussSeidelIteration(std::vector<float> &p, const std::vector<float> &div);
	//Project投影
	void Project();
	//AddBouyancySmoke烟浮力
	void AddBouyancySmoke();
	//AddBouyancyCloud烟浮力
	void AddBouyancyCloud();
	//VorticityConfinement涡旋限制
	void VorticityConfinement();
	//PhaseTransitionCloud云相变
	void PhaseTransitionCloud();
	//CorrectTimestep时间步长纠正(是否合理？)
	void CorrectTimestep();
	//SourceControl源控制
	void SourceControl();
};