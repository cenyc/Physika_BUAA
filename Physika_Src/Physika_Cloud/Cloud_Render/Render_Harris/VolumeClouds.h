#ifndef _VOLUMETRIC_H
#define _VOLUMETRIC_H

#define SORT_TOWARD 0
#define SORT_AWAY 1

struct CloudPuff
{
	float 	Size;
	int 	ID;
	float 	Angle;
	Vector3 Position;
	float	DistanceToCam;
	Color4 	Color;
	float	Life;	
};

struct VolumetricCloud
{
	vector<CloudPuff> Puffs;
	unsigned 	ImpostorTex;
	Vector3		BoundingBox1, BoundingBox2;
	Vector3 	Center;
	float		Radius;
	Vector3*	VertexBuffer;
	Color4*		ColorBuffer;
	Vector2*	TexCoordBuffer;
	
	Vector3		LastCamera;
	Vector3		LastLight;
	float		DistanceFromCamera;
	int			ImpostorSize;
	Vector3		vx, vy;	//up and right vectors for impostor rendering
};

class VolumetricClouds
{
public:
	VolumetricClouds();
	
	int 	Create(int NumClouds, float PlaneSize, float PlaneHeight);	
	void	Update(Vector3 Sun, Vector3 Camera);
	void	Render(Vector3 Sun, Vector3 Camera);
	void	Destroy();
	void	GetInfo(int* Sprites, int* Impostors)
	{
		*Sprites = NumSprites;
		*Impostors = NumImpostors;
	}

private:

	void	UpdateCloud(VolumetricCloud* Cloud, Vector3 Sun, Vector3 Camera);
	void	RenderCloudImpostor(VolumetricCloud* Cloud, float alpha); 
	void	RenderCloud3D(VolumetricCloud* Cloud, Vector3 Camera, Vector3 Sun, float alpha);

	void	MakeCloudImpostor(VolumetricCloud* Cloud, Vector3 Sun, Vector3 Camera);
	
	void	LightCloud(VolumetricCloud* Cloud, Vector3 Sun);
	void	GrowCloud(VolumetricCloud* Cloud, int level, float radius, Vector3 Position);

	int		GetImpostorSize(float distance2);	
	void	GenerateTexture();
	void	SortParticles(VolumetricCloud* Cloud, int mode);

	class	SortAwayComparison
	{
	public:
		bool operator () (CloudPuff puff1, CloudPuff puff2)
		{
			return puff1.DistanceToCam < puff2.DistanceToCam;
		}
	} SortAway;

	class	SortTowardComparison
	{
	public:
		bool operator () (CloudPuff puff1, CloudPuff puff2)
		{
			return puff1.DistanceToCam > puff2.DistanceToCam;
		}
	} SortToward;

	class	SortCloudsTowardComparison
	{
	public:
		bool operator () (VolumetricCloud cloud1, VolumetricCloud cloud2)
		{
			return cloud1.DistanceFromCamera > cloud2.DistanceFromCamera;
		};
	} SortCloudToward;

	unsigned PuffTexture;
	unsigned PuffImage;

	vector<VolumetricCloud> Clouds;
	
	int	SplatBufferSize;
	int ImpostorSize;
	int NumSprites, NumImpostors;
	float Albedo, Extinction;
};

#endif
