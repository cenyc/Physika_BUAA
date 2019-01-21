#ifndef _CAM_H
#define _CAM_H

class CCamera {

  public:

	CCamera();

	Vector3 UpVector() {	return m_vUpVector;		}

	Vector3 GetCameraView() { return m_vView; }
	Vector3 GetCameraPosition() { return m_vPosition; }
		
	void Create(float fov, int Width, int Height, float maxd);
	
	void PositionCamera(float positionX, float positionY, float positionZ,
			 		    float viewX,     float viewY,     float viewZ,
						float upVectorX, float upVectorY, float upVectorZ);

	void SetViewByMouse(); 

	void RotateAroundPoint(Vector3 vCenter, float angle, float X, float Y, float Z);
	void RotateView(float angle, float x, float y, float z);

	
	
	void CheckForMovement();
	void StrafeCamera(float);
	void MoveCamera(float);
	void Update();
	void Look();
	void LookAt(Vector3 v);
	void SetDist(float v) { fDist = v; }
	void Reset();
	void SetSpeed(float sp) { kSpeed = sp; }

public:

	float	kSpeed;

	Vector3 m_vPosition;					
	Vector3 m_vView;						
	Vector3 m_vUpVector;		
	Vector3 m_lookAt;
	Vector3 m_vStrafe;
	
	float fDist;
	float currentRotX;
	
	int 	ScreenWidth, ScreenHeight;			
};

#endif
