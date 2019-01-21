

#ifndef _MATRIX_H

#define _MATRIX_H

class Matrix44
{
public:

	Matrix44();

	Matrix44  (float m00, float m01, float m02, float m03,
			  float m10, float m11, float m12, float m13,
			  float m20, float m21, float m22, float m23,
			  float m30, float m31, float m32, float m33);

	Matrix44 operator * (Matrix44 mat);
	Matrix44 operator + (Matrix44 mat);
	Matrix44 operator - (Matrix44 mat);
	void operator = (Matrix44 mat);
	void operator = (Vector3 v);

	void RotateX(float angle);
	void RotateY(float angle);
	void RotateZ(float angle);
	void RotateAxis(Vector3 Axis, float angle);
	void Translate(float x, float y, float z);
	void Scale(float x, float y, float z);
	void LoadIdentity();

	Vector3 TransformVertex(Vector3 v);

	void GetGLMatrix(float*);
	void LoadGLMatrix(float *);

	void Transpose();

	float m[4][4];
};

#endif
