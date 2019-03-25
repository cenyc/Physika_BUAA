#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<fstream>

//�Ѿ�Ĭ�����񳤶�Ϊ1
#define GRID_LENGTH 1
//�жϸ������Ƿ�������ʵ���ֵ
#define F_MIN 1e-6

//------------demo�ĳ�ʼ����(Ĭ�Ϲ��캯��)-----------------
//��ʼ���������/��С�¶�
#define MAX_AMBIENT_TEMP 25.0
#define MIN_AMBIENT_TEMP 5.0
//��ʼ������ˮ���ܶ�
#define INIT_VAPOR_DENSITY 1.0f
//��ʼ������ˮ���¶�
#define INIT_MATERIAL_TEMP 15.0
//��ʼ������ˮ���������ٶ�
#define INIT_VELOCITY_W 0.3
//��ʼ������ˮ����Χ
#define BEGIN_X 28
#define END_X 37
#define BEGIN_Y 28
#define END_Y 37
//---------------------------------------------------------

class CloudGen
{
public:
	//����ģʽ(0 - ��, 1 - ��)
	CloudGen(int mode_flag);
	//�����ļ�·��(�������ļ�����,������n_��source_change_times_��Source)
	CloudGen(std::string input_filename, int mode_flag);
	~CloudGen();

	//ֱ�����ɶ�Ӧ���������ս��(�������д���)
	void CloudGenFinalRun(int times);
	//����demoչʾÿ֡����()
	void CloudGenFrameRun();

	//---------����demo��ʾ-----------
	const float GetVelocityU(int i, int j, int k);
	const float GetVelocityV(int i, int j, int k);
	const float GetVelocityW(int i, int j, int k);
	const float GetVaporDensity(int i, int j, int k);
	const float GetCloudDensity(int i, int j, int k);
	const float GetAmbientTemp(int k);
	const float GetMaterialTemp(int i, int j, int k);
	const int GetFrameCount();

private:
	int n_;  //ģ��ռ��С(�������߽�)
	int size_;
	float dt_;   //ʱ�䲽��
	int mode_flag;   //����ģʽ(0-��,1-��)
	int frame_count_;   //���ڼ�¼���е��ڼ�֡

	int source_change_times_;  //���ڼ�¼Դ�仯��֡��
	struct Source
	{
		float vapor_density;    //Դ(����)ˮ���ܶ�
		float material_temp;    //Դ(����)ˮ���¶�
		float velocity_w;       //Դ(����)ˮ����ʼ�ٶ�
	};
	std::vector<std::vector<Source>>  source_;   //֡��*ģ��ռ���棨x*y����С*ÿ�������Դ�����

	std::vector<float> velocity_u_;   //x�����ٶ�
	std::vector<float> velocity_v_;   //y�����ٶ�
	std::vector<float> velocity_w_;   //z�����ٶ�
	std::vector<float> vapor_density_;  //ˮ���ܶ�
	std::vector<float> cloud_density_;  //���ܶ�
	std::vector<float> ambient_temp_;  //�����¶�
	std::vector<float> material_temp_;  //�����¶�

	//��ʼ�������ļ�����(���а���:)
	void InitFileInput(std::string &filename);
	//Դ(����)��ά����λ��ת��һά��������
	const int SourcePosition(int i, int j);
	//�ռ���ά����λ��ת��һά��������
	const int Position(int i, int j, int k);
	//�����ٶȳ�
	void GetVelocity();
	//�����ܶȳ�
	void GetDensity();
	//�߽����(value:���ƵĶ�����  velocity_project_flag:�ٶ�ͶӰ���(0-���ǣ�1-x�ᣬ2-y�ᣬ3-z��)
	void BoundaryCondition(std::vector<float> &value, const int velocity_project_flag);
	//Advectƽ��(Ҫƽ��������������һ֡��ֵ������Ҫ������������ٶ�)
	void Advect(std::vector<float> &value, const std::vector<float> &pre_value, const std::vector<float> &velocity_u, const std::vector<float> &velocity_v, const std::vector<float> &velocity_w);
	//��˹���¶�����
	void GaussSeidelIteration(std::vector<float> &p, const std::vector<float> &div);
	//ProjectͶӰ
	void Project();
	//AddBouyancySmoke�̸���
	void AddBouyancySmoke();
	//AddBouyancyCloud�̸���
	void AddBouyancyCloud();
	//VorticityConfinement��������
	void VorticityConfinement();
	//PhaseTransitionCloud�����
	void PhaseTransitionCloud();
	//CorrectTimestepʱ�䲽������(�Ƿ����)
	void CorrectTimestep();
	//SourceControlԴ����
	void SourceControl();
};