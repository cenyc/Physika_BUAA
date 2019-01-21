#include "Common.h"

CCamera::CCamera()
{
	Reset();
	kSpeed = 2.0f;
}

void CCamera::Create(float fov, int Width, int Height, float maxd)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fov, (float)Width / (float)Height, 1.0f, maxd);
	glMatrixMode(GL_MODELVIEW);
	
	ScreenWidth = Width;
	ScreenHeight = Height;
}

void CCamera::Reset()
{
	Vector3 vZero = Vector3(0.0, 0.0, 0.0);		
	Vector3 vView = Vector3(1.0, 0.0, 0.0);		
	Vector3 vUp   = Vector3(0.0, 1.0, 0.0);		

	m_vPosition	= vZero;					
	m_vView		= vView;					
	m_vUpVector	= vUp;	
	
	currentRotX = 0.0f;
}

void CCamera::PositionCamera(float positionX, float positionY, float positionZ,
				  		     float viewX,     float viewY,     float viewZ,
							 float upVectorX, float upVectorY, float upVectorZ)
{
	Vector3 vPosition	= Vector3(positionX, positionY, positionZ);
	Vector3 vView		= Vector3(viewX, viewY, viewZ);
	Vector3 vUpVector	= Vector3(upVectorX, upVectorY, upVectorZ);

	m_vPosition = vPosition;					
	m_vView     = vView;						
	m_vUpVector = vUpVector;					
}


void CCamera::SetViewByMouse()//�ӵ����������ƶ����ƶ�
{
	unsigned short middleX = ScreenWidth / 2;
	unsigned short middleY = ScreenHeight / 2;
	
	float angleY = 0.0f;
	float angleZ = 0.0f;
	
	int PosX, PosY;
	SDL_GetMouseState(&PosX, &PosY);	
	if( (PosX == middleX) && (PosY == middleY) ) 	return;
	SDL_WarpMouse(middleX, middleY);
	
	angleY = (float)( (middleX - PosX) ) / 500.0f;		
	angleZ = (float)( (middleY - PosY) ) / 500.0f;		

	currentRotX -= angleZ;  

	Vector3 vAxis = Cross(m_vView - m_vPosition, m_vUpVector);//???????????
	vAxis = Normalize(vAxis);

	RotateView(angleY,0,1, 0);

	if(currentRotX > 2.0f)
		currentRotX = 2.0f;	
	else 
		if(currentRotX < -2.0f)
			currentRotX = -2.0f;	
		else
			RotateView(angleZ, vAxis.x, vAxis.y, vAxis.z);//?????????????
	        //RotateView(angleZ, 0, 0, 1);
}

void CCamera::RotateAroundPoint(Vector3 vCenter, float angle, float x, float y, float z)
{
	Vector3 vNewPosition;			

	Vector3 vPos = m_vPosition - vCenter;
	
	float cosTheta = (float)cos(angle);
	float sinTheta = (float)sin(angle);

	vNewPosition.x  = (cosTheta + (1 - cosTheta) * x * x)		* vPos.x;
	vNewPosition.x += ((1 - cosTheta) * x * y - z * sinTheta)	* vPos.y;
	vNewPosition.x += ((1 - cosTheta) * x * z + y * sinTheta)	* vPos.z;
	
	vNewPosition.y  = ((1 - cosTheta) * x * y + z * sinTheta)	* vPos.x;
	vNewPosition.y += (cosTheta + (1 - cosTheta) * y * y)		* vPos.y;
	vNewPosition.y += ((1 - cosTheta) * y * z - x * sinTheta)	* vPos.z;
	
	vNewPosition.z  = ((1 - cosTheta) * x * z - y * sinTheta)	* vPos.x;
	vNewPosition.z += ((1 - cosTheta) * y * z + x * sinTheta)	* vPos.y;
	vNewPosition.z += (cosTheta + (1 - cosTheta) * z * z)		* vPos.z;

	m_vPosition = vCenter + vNewPosition;
}

void CCamera::RotateView(float angle, float x, float y, float z)
{
	Vector3 vNewView;

	// Get the view vector (The direction we are facing)
	Vector3 vView = m_vView - m_vPosition;		

	// Calculate the sine and cosine of the angle once
	float cosTheta = (float)cos(angle);
	float sinTheta = (float)sin(angle);

	// Find the new x position for the new rotated point
	vNewView.x  = (cosTheta + (1 - cosTheta) * x * x)		* vView.x;
	vNewView.x += ((1 - cosTheta) * x * y - z * sinTheta)	* vView.y;
	vNewView.x += ((1 - cosTheta) * x * z + y * sinTheta)	* vView.z;

	// Find the new y position for the new rotated point
	vNewView.y  = ((1 - cosTheta) * x * y + z * sinTheta)	* vView.x;
	vNewView.y += (cosTheta + (1 - cosTheta) * y * y)		* vView.y;
	vNewView.y += ((1 - cosTheta) * y * z - x * sinTheta)	* vView.z;

	// Find the new z position for the new rotated point
	vNewView.z  = ((1 - cosTheta) * x * z - y * sinTheta)	* vView.x;
	vNewView.z += ((1 - cosTheta) * y * z + x * sinTheta)	* vView.y;
	vNewView.z += (cosTheta + (1 - cosTheta) * z * z)		* vView.z;

	// Now we just add the newly rotated vector to our position to set
	// our new rotated view of our camera.
	m_vView = m_vPosition + vNewView;
}

void CCamera::StrafeCamera(float speed)
{	
	// Add the strafe vector to our position
	m_vPosition.x += m_vStrafe.x * speed;
	m_vPosition.z += m_vStrafe.z * speed;
	m_vPosition.y += m_vStrafe.y * speed;

	// Add the strafe vector to our view
	m_vView.x += m_vStrafe.x * speed;
	m_vView.z += m_vStrafe.z * speed;
	m_vView.y += m_vStrafe.y * speed;

}


void CCamera::MoveCamera(float speed)
{
	// Get the current view vector (the direction we are looking)
	Vector3 vVector = m_vView - m_vPosition;
	vVector = Normalize(vVector);

	m_vPosition.x += vVector.x * speed;
	m_vPosition.y += vVector.y * speed;
	m_vPosition.z += vVector.z * speed;
	m_vView.x += vVector.x * speed;
	m_vView.y += vVector.y * speed;
	m_vView.z += vVector.z * speed;
}



void CCamera::CheckForMovement()
{	
	float speed = kSpeed;

	unsigned char *Keys = SDL_GetKeyState(NULL);
	
	if (Keys[SDLK_UP])	MoveCamera(speed);				
	if (Keys[SDLK_DOWN]) MoveCamera(-speed);				
	if (Keys[SDLK_LEFT]) StrafeCamera(-speed);
	if (Keys[SDLK_RIGHT]) StrafeCamera(speed);
}

void CCamera::LookAt(Vector3 v)
{
	m_vView = v;
	if (SqDist(m_vView, m_vPosition) > fDist)
		m_vPosition = m_vView + Normalize(m_vPosition - m_vView) * fDist;
}

void CCamera::Update() 
{
	Vector3 vCross = Cross(m_vView - m_vPosition, m_vUpVector);
	m_vStrafe = Normalize(vCross);

	CheckForMovement();
	SetViewByMouse();
}


void CCamera::Look()
{
	
	gluLookAt(m_vPosition.x, m_vPosition.y, m_vPosition.z,	
			  m_vView.x,	 m_vView.y,     m_vView.z,	
			  m_vUpVector.x, m_vUpVector.y, m_vUpVector.z);
}


