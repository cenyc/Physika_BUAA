

#include "Common.h"

Matrix44::Matrix44()
{
	m[0][0] = 0.0f; m[0][1] = 0.0f; m[0][2] = 0.0f; m[0][3] = 0.0f;
	m[1][0] = 0.0f; m[1][1] = 0.0f; m[1][2] = 0.0f; m[1][3] = 0.0f;
	m[2][0] = 0.0f; m[2][1] = 0.0f; m[2][2] = 0.0f; m[2][3] = 0.0f;
	m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 0.0f;
}

Matrix44::Matrix44(float m00, float m01, float m02, float m03,
					 float m10, float m11, float m12, float m13,
					 float m20, float m21, float m22, float m23,
					 float m30, float m31, float m32, float m33)
{
	m[0][0] = m00; m[0][1] = m01; m[0][2] = m02; m[0][3] = m03;
	m[1][0] = m10; m[1][1] = m11; m[1][2] = m12; m[1][3] = m13;
	m[2][0] = m20; m[2][1] = m21; m[2][2] = m22; m[2][3] = m23;
	m[3][0] = m30; m[3][1] = m31; m[3][2] = m32; m[3][3] = m33;
}

Matrix44 Matrix44::operator +(Matrix44 mat)
{
	return Matrix44(m[0][0] + mat.m[0][0], m[0][1] + mat.m[0][1], m[0][2] + mat.m[0][2], m[0][3] + mat.m[0][3],
					 m[1][0] + mat.m[1][0], m[1][1] + mat.m[1][1], m[1][2] + mat.m[1][2], m[1][3] + mat.m[1][3],
					 m[2][0] + mat.m[2][0], m[2][1] + mat.m[2][1], m[2][2] + mat.m[2][2], m[2][3] + mat.m[2][3],
					 m[3][0] + mat.m[3][0], m[3][1] + mat.m[3][1], m[3][2] + mat.m[3][2], m[3][3] + mat.m[3][3]);
}

Matrix44 Matrix44::operator -(Matrix44 mat)
{
	return Matrix44(m[0][0] - mat.m[0][0], m[0][1] - mat.m[0][1], m[0][2] - mat.m[0][2], m[0][3] - mat.m[0][3],
					 m[1][0] - mat.m[1][0], m[1][1] - mat.m[1][1], m[1][2] - mat.m[1][2], m[1][3] - mat.m[1][3],
					 m[2][0] - mat.m[2][0], m[2][1] - mat.m[2][1], m[2][2] - mat.m[2][2], m[2][3] - mat.m[2][3],
					 m[3][0] - mat.m[3][0], m[3][1] - mat.m[3][1], m[3][2] - mat.m[3][2], m[3][3] - mat.m[3][3]);
}

Matrix44 Matrix44::operator *(Matrix44 mat)
{
	return Matrix44(
		mat.m[0][0] * m[0][0] + mat.m[0][1] * m[1][0] + mat.m[0][2] * m[2][0] + mat.m[0][3] * m[3][0],
		mat.m[0][0] * m[0][1] + mat.m[0][1] * m[1][1] + mat.m[0][2] * m[2][1] + mat.m[0][3] * m[3][1],
		mat.m[0][0] * m[0][2] + mat.m[0][1] * m[1][2] + mat.m[0][2] * m[2][2] + mat.m[0][3] * m[3][2],
		mat.m[0][0] * m[0][3] + mat.m[0][1] * m[1][3] + mat.m[0][2] * m[2][3] + mat.m[0][3] * m[3][3],
			   
		mat.m[1][0] * m[0][0] + mat.m[1][1] * m[1][0] + mat.m[1][2] * m[2][0] + mat.m[1][3] * m[3][0],
		mat.m[1][0] * m[0][1] + mat.m[1][1] * m[1][1] + mat.m[1][2] * m[2][1] + mat.m[1][3] * m[3][1],
		mat.m[1][0] * m[0][2] + mat.m[1][1] * m[1][2] + mat.m[1][2] * m[2][2] + mat.m[1][3] * m[3][2],
		mat.m[1][0] * m[0][3] + mat.m[1][1] * m[1][3] + mat.m[1][2] * m[2][3] + mat.m[1][3] * m[3][3],

		mat.m[2][0] * m[0][0] + mat.m[2][1] * m[1][0] + mat.m[2][2] * m[2][0] + mat.m[2][3] * m[3][0],
		mat.m[2][0] * m[0][1] + mat.m[2][1] * m[1][1] + mat.m[2][2] * m[2][1] + mat.m[2][3] * m[3][1],
		mat.m[2][0] * m[0][2] + mat.m[2][1] * m[1][2] + mat.m[2][2] * m[2][2] + mat.m[2][3] * m[3][2],
		mat.m[2][0] * m[0][3] + mat.m[2][1] * m[1][3] + mat.m[2][2] * m[2][3] + mat.m[2][3] * m[3][3],

		mat.m[3][0] * m[0][0] + mat.m[3][1] * m[1][0] + mat.m[3][2] * m[2][0] + mat.m[3][3] * m[3][0],
		mat.m[3][0] * m[0][1] + mat.m[3][1] * m[1][1] + mat.m[3][2] * m[2][1] + mat.m[3][3] * m[3][1],
		mat.m[3][0] * m[0][2] + mat.m[3][1] * m[1][2] + mat.m[3][2] * m[2][2] + mat.m[3][3] * m[3][2],
		mat.m[3][0] * m[0][3] + mat.m[3][1] * m[1][3] + mat.m[3][2] * m[2][3] + mat.m[3][3] * m[3][3]);
}


void Matrix44::operator =(Matrix44 mat)
{
	memcpy(m, mat.m, sizeof(float) * 16);
}

void Matrix44::operator =(Vector3 v)
{
	LoadIdentity();
	m[0][3] = v.x;
	m[1][3] = v.y;
	m[2][3] = v.z;
}

void Matrix44::RotateAxis(Vector3 Axis, float angle)
{
	float u = Axis.x, v = Axis.y, w = Axis.z;
	
	float L2 = Dot(Axis, Axis);
	float L = sqrtf(L2);
	float cosA = cos(angle);
	float sinA = cos(angle);
	
	Matrix44 rot;
	
	rot.m[0][0] = (u * u + (v * v + w * w) * cosA) / L2;
	rot.m[0][1] = (u * v * (1 - cosA) - w * L * sinA) / L2;
	rot.m[0][2] = (u * w * (1 - cosA) + v * L * sinA) / L2;
	rot.m[0][3] = 0;
	rot.m[1][0] = (u * v * (1 - cosA) + w * L * sinA) / L2;
	rot.m[1][1] = (v * v + (u * u + w * w) * cosA) / L2;
	rot.m[1][2] = (v * w * (1 - cosA) - u * L * sinA) / L2;
	rot.m[1][3] = 0;
	rot.m[2][0] = (u * w * (1 - cosA) - v * L * sinA) / L2;
	rot.m[2][1] = (v * w * (1 - cosA) + u * L * sinA) / L2;
	rot.m[2][2] = (w * w + (u * u + v * v) * cosA) / L2;
	rot.m[2][3] = 0;
	rot.m[3][0] = 0;
	rot.m[3][1] = 0;
	rot.m[3][2] = 0;
	rot.m[3][3] = 1;
	
	(*this) = (*this) * rot;
}

void Matrix44::RotateX(float angle)
{
	Matrix44 rot(1,          0,           0, 0,
		          0, cos(angle), -sin(angle), 0,
				  0, sin(angle),  cos(angle), 0,
				  0,          0,           0, 1);

	(*this) = (*this) * rot;
}

void Matrix44::RotateY(float angle)
{
	Matrix44 rot(cos(angle), 0, -sin(angle), 0,
		                 0, 1,           0, 0,
				sin(angle), 0,  cos(angle), 0,
				         0, 0,           0, 1);

	(*this) = (*this) * rot;
}

void Matrix44::RotateZ(float angle)
{
	Matrix44 rot(cos(angle), -sin(angle), 0, 0,
		          sin(angle),  cos(angle), 0, 0,
				           0,           0, 1, 0,
				           0,           0, 0, 1);

	(*this) = (*this) * rot;
}

void Matrix44::Translate(float x, float y, float z)
{
	Matrix44 tran(1, 0, 0, x,
		           0, 1, 0, y,
				   0, 0, 1, z,
				   0, 0, 0, 1);
	(*this) = (*this) * tran;
}

void Matrix44::Scale(float x, float y, float z)
{
	Matrix44 sc(x, 0, 0, 0,
		         0, y, 0, 0,
				 0, 0, z, 0,
				 0, 0, 0, 1);

	(*this) = (*this) * sc;
}

void Matrix44::LoadIdentity()
{
	(*this) = Matrix44(1, 0, 0, 0,
		                0, 1, 0, 0,
						0, 0, 1, 0,
						0, 0, 0, 1);
}
void Matrix44::GetGLMatrix(float* v)
{
	v[ 0] = m[0][0];
	v[ 1] = m[1][0];
	v[ 2] = m[2][0];
	v[ 3] = m[3][0];
	v[ 4] = m[0][1];
	v[ 5] = m[1][1];
	v[ 6] = m[2][1];
	v[ 7] = m[3][1];
	v[ 8] = m[0][2];
	v[ 9] = m[1][2];
	v[10] = m[2][2];
	v[11] = m[3][2];
	v[12] = m[0][3];
	v[13] = m[1][3];
	v[14] = m[2][3];
	v[15] = m[3][3];
}

void Matrix44::LoadGLMatrix(float *v)
{
	m[0][0] = v[ 0];
	m[1][0] = v[ 1];
	m[2][0] = v[ 2];
	m[3][0] = v[ 3];
	m[0][1] = v[ 4];
	m[1][1] = v[ 5];
	m[2][1] = v[ 6];
	m[3][1] = v[ 7];
	m[0][2] = v[ 8];
	m[1][2] = v[ 9];
	m[2][2] = v[10];
	m[3][2] = v[11];
	m[0][3] = v[12];
	m[1][3] = v[13];
	m[2][3] = v[14];
	m[3][3] = v[15];
}

Vector3 Matrix44::TransformVertex(Vector3 v)
{
	Vector3 result;

	result.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3];
	result.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3];
	result.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3];

	return result;
}

void Matrix44::Transpose()
{
	m[0][0] = m[0][0];
	m[0][1] = m[1][0];
	m[0][2] = m[2][0];
	m[0][3] = m[3][0];
	m[1][0] = m[0][1];
	m[1][1] = m[1][1];
	m[1][2] = m[2][1];
	m[1][3] = m[3][1];
	m[2][0] = m[0][2];
	m[2][1] = m[1][2];
	m[2][2] = m[2][2];
	m[2][3] = m[3][2];
	m[3][0] = m[0][3];
	m[3][1] = m[1][3];
	m[3][2] = m[2][3];
	m[3][3] = m[3][3];
}
