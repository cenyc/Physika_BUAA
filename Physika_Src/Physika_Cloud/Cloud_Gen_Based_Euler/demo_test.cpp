#include"Cloud_gen.h"

int main(int argc, char **argv){
	CloudGen cloud("C:/Users/李云飞/Desktop/test/1.dat", 0);
	cloud.CloudGenFinalRun(20);
	return 0;
}

//---------------------------------------------------------
/*

#define GLEW_STATIC

#include<GL\glut.h>
#include"Cloud_gen.h"


//CloudGen cloud(0);
CloudGen cloud("C:/Users/李云飞/Desktop/test/1.dat", 0);

const int windowWidth = 1000;
const int windowHeight = 800;

GLfloat translationX;
GLfloat translationY;
GLfloat translationZ;

GLfloat rotationX;
GLfloat rotationY;
GLfloat rotationZ;

GLfloat alpha = 0.05;

int N = 64;

bool drawVelocity = false;
//输出控制
bool output = false;

void Init(){
	//------------------------------
	rotationX = -90.0f;
	rotationY = 0;
	rotationZ = 0;

	translationX = -1.0;
	translationY = -1.0;
	translationZ = -1.0;
}

void DrawGrid(){
	glLineWidth(1.0f);

	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(1.3f, 0.0f, 0.0f);

	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 1.3f, 0.0f);

	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 1.3f);

	glVertex3f(1.3f, 0.0f, 0.0f);
	glVertex3f(1.3f, 1.3f, 0.0f);

	glVertex3f(1.3f, 1.3f, 0.0f);
	glVertex3f(0.0f, 1.3f, 0.0f);

	glVertex3f(0.0f, 1.3f, 1.3f);
	glVertex3f(0.0f, 0.0f, 1.3f);

	glVertex3f(0.0f, 1.3f, 1.3f);
	glVertex3f(0.0f, 1.3f, 0.0f);

	glVertex3f(1.3f, 0.0f, 0.0f);
	glVertex3f(1.3f, 0.0f, 1.3f);

	glVertex3f(0.0f, 0.0f, 1.3f);
	glVertex3f(1.3f, 0.0f, 1.3f);

	glVertex3f(1.3f, 1.3f, 0.0f);
	glVertex3f(1.3f, 1.3f, 1.3f);

	glVertex3f(1.3f, 1.3f, 1.3f);
	glVertex3f(1.3f, 0.0f, 1.3f);

	glVertex3f(0.0f, 1.3f, 1.3f);
	glVertex3f(1.3f, 1.3f, 1.3f);

	glEnd();
}

void DrawDensity(){
	GLfloat positionX;
	GLfloat positionY;
	GLfloat positionZ;

	GLfloat density000;
	GLfloat density010;
	GLfloat density100;
	GLfloat density110;
	GLfloat density001;
	GLfloat density011;
	GLfloat density101;
	GLfloat density111;

	GLfloat h = 1.3f / N;

	glBegin(GL_QUADS);
	for (int x = 1; x <= N; x++)
	{
		positionX = (x - 0.5f) * h;

		for (int y = 1; y <= N; y++)
		{
			positionY = (y - 0.5f) * h;

			for (int z = 1; z <= N; z++)
			{
				positionZ = (z - 0.5f) * h;

				if (cloud.GetVaporDensity(x, y, z)>0.0f){

					//水
					density000 = cloud.GetVaporDensity(x, y, z) / INIT_VAPOR_DENSITY;
					density010 = cloud.GetVaporDensity(x, y + 1, z) / INIT_VAPOR_DENSITY;
					density100 = cloud.GetVaporDensity(x + 1, y, z) / INIT_VAPOR_DENSITY;
					density110 = cloud.GetVaporDensity(x + 1, y + 1, z) / INIT_VAPOR_DENSITY;

					density001 = cloud.GetVaporDensity(x, y, z + 1) / INIT_VAPOR_DENSITY;
					density011 = cloud.GetVaporDensity(x, y + 1, z + 1) / INIT_VAPOR_DENSITY;
					density101 = cloud.GetVaporDensity(x + 1, y, z + 1) / INIT_VAPOR_DENSITY;
					density111 = cloud.GetVaporDensity(x + 1, y + 1, z + 1) / INIT_VAPOR_DENSITY;

					////云
					//density000 = cloud.GetCloudDensity(x, y, z);
					//density010 = cloud.GetCloudDensity(x, y + 1, z);
					//density100 = cloud.GetCloudDensity(x + 1, y, z);
					//density110 = cloud.GetCloudDensity(x + 1, y + 1, z);

					//density001 = cloud.GetCloudDensity(x, y, z + 1);
					//density011 = cloud.GetCloudDensity(x, y + 1, z + 1);
					//density101 = cloud.GetCloudDensity(x + 1, y, z + 1);
					//density111 = cloud.GetCloudDensity(x + 1, y + 1, z + 1);

					glColor4f(density111, density111, density111, alpha);
					glVertex3f(positionX + h, positionY + h, positionZ + h);

					glColor4f(density011, density011, density011, alpha);
					glVertex3f(positionX, positionY + h, positionZ + h);

					glColor4f(density001, density001, density001, alpha);
					glVertex3f(positionX, positionY, positionZ + h);

					glColor4f(density101, density101, density101, alpha);
					glVertex3f(positionX + h, positionY, positionZ + h);

					glColor4f(density110, density110, density110, alpha);
					glVertex3f(positionX + h, positionY + h, positionZ);

					glColor4f(density111, density111, density111, alpha);
					glVertex3f(positionX + h, positionY + h, positionZ + h);

					glColor4f(density101, density101, density101, alpha);
					glVertex3f(positionX + h, positionY, positionZ + h);

					glColor4f(density100, density100, density100, alpha);
					glVertex3f(positionX + h, positionY, positionZ);

					glColor4f(density010, density010, density010, alpha);
					glVertex3f(positionX, positionY + h, positionZ);

					glColor4f(density110, density110, density110, alpha);
					glVertex3f(positionX + h, positionY + h, positionZ);

					glColor4f(density100, density100, density100, alpha);
					glVertex3f(positionX + h, positionY, positionZ);

					glColor4f(density000, density000, density000, alpha);
					glVertex3f(positionX, positionY, positionZ);

					glColor4f(density011, density011, density011, alpha);
					glVertex3f(positionX, positionY + h, positionZ + h);

					glColor4f(density010, density010, density010, alpha);
					glVertex3f(positionX, positionY + h, positionZ);

					glColor4f(density000, density000, density000, alpha);
					glVertex3f(positionX, positionY, positionZ);

					glColor4f(density001, density001, density001, alpha);
					glVertex3f(positionX, positionY, positionZ + h);

					glColor4f(density100, density100, density100, alpha);
					glVertex3f(positionX + h, positionY, positionZ);

					glColor4f(density000, density000, density000, alpha);
					glVertex3f(positionX, positionY, positionZ);

					glColor4f(density001, density001, density001, alpha);
					glVertex3f(positionX, positionY, positionZ + h);

					glColor4f(density101, density101, density101, alpha);
					glVertex3f(positionX + h, positionY, positionZ + h);

					glColor4f(density110, density110, density110, alpha);
					glVertex3f(positionX + h, positionY + h, positionZ);

					glColor4f(density010, density010, density010, alpha);
					glVertex3f(positionX, positionY + h, positionZ);

					glColor4f(density011, density011, density011, alpha);
					glVertex3f(positionX, positionY + h, positionZ + h);

					glColor4f(density111, density111, density111, alpha);
					glVertex3f(positionX + h, positionY + h, positionZ + h);
				}
			}
		}
	}
	glEnd();
}

void DrawVelocity(){
	GLfloat positionX;
	GLfloat positionY;
	GLfloat positionZ;

	GLfloat h = 1.3f / N;
	glColor3f(1.0, 1.0, 1.0);

	for (int x = 1; x <= N; x++)
	{
		positionX = (x - 0.5f) * h;
		for (int y = 1; y <= N; y++)
		{
			positionY = (y - 0.5f) * h;
			for (int z = 1; z <= N; z++)
			{
				positionZ = (z - 0.5f) * h;
				glBegin(GL_LINES);
				glVertex3f(positionX, positionY, positionZ);
				glVertex3f(positionX + cloud.GetVelocityU(x, y, z) / 2,
					positionY + cloud.GetVelocityV(x, y, z) / 2,
					positionZ + cloud.GetVelocityW(x, y, z) / 2);
				glEnd();
			}
		}
	}
}

void Display(){
	//若使用RGBA模式指定颜色，则需要激活OpenGL颜色调和操作
	glEnable(GL_BLEND);
	//描述颜色调和方法（源调和因子：将四个调和因子设为源alpha；目标调和因子：( 1.0, 1.0, 1.0, 1.0)，即完全使用颜色参与混合运算）
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);

	//启用alpha测试
	glEnable(GL_ALPHA_TEST);
	//设置Alpha测试条件为大于0则通过
	glAlphaFunc(GL_GREATER, 0);

	//********
	glClear(GL_COLOR_BUFFER_BIT);


	//复制活动栈顶的当前状态矩阵并将其存入第二个栈位置
	glPushMatrix();
	glTranslatef(translationX, translationY, translationZ);
	glRotatef(rotationX, 1.0f, 0, 0);
	glRotatef(rotationZ, 0, 0, 1.0f);

	if (drawVelocity) DrawVelocity();
	else DrawDensity();

	DrawGrid();

	//出栈恢复原矩阵？？
	glPopMatrix();
	glutSwapBuffers();
}

void SpecialKeys(int key, int x, int y)
{
	if (key == GLUT_KEY_UP)
	{
		rotationX -= 5.0f;

		//rotationX = (GLfloat)((const int)xRot60);
		//rotationY=(GLfloat)((const int)yRot60);
	}

	if (key == GLUT_KEY_DOWN)
	{
		rotationX += 5.0f;

		//rotationX=(GLfloat)((const int)rotationX * 60);
		//rotationY=(GLfloat)((const int)ro);
	}

	if (key == GLUT_KEY_LEFT)
	{
		rotationZ += 5.0f;

		//rotationX=(GLfloat)((const int)xRot60);
		//rotationY=(GLfloat)((const int)yRot60);
	}

	if (key == GLUT_KEY_RIGHT)
	{
		rotationZ -= 5.0f;

		//rotationX=(GLfloat)((const int)xRot60);
		//rotationY=(GLfloat)((const int)yRot60);
	}

	glutPostRedisplay();
}

void Reshape(GLint width, GLint height){
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(45.0, (float)width / height, 0.001, 100.0);
	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -3.0f);
}

void KeyBoard(unsigned char key, int x, int y)
{
	if (key == 'w' || key == 'W')
	{
		translationZ = translationZ + 0.1;
	}

	if (key == 's' || key == 'S')
	{
		translationZ = translationZ - 0.1;
	}

	if (key == 'a' || key == 'A')
	{
		translationX = translationX - 0.1;
	}

	if (key == 'd' || key == 'D')
	{
		translationX = translationX + 0.1;
	}

	if (key == 'q' || key == 'Q')
	{
		translationY = translationY + 0.1;
	}

	if (key == 'e' || key == 'E')
	{
		translationY = translationY - 0.1;
	}

	if (key == 'v' || key == 'V')
	{
		drawVelocity = !drawVelocity;
	}
}

void Idle()
{
	cloud.CloudGenFrameRun();
	std::cout << "第" << cloud.GetFrameCount() << "帧" << std::endl;
	glutPostRedisplay();
}

int main(int argc, char **argv){
	glutInit(&argc, argv);
	Init();
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(windowWidth, windowHeight);
	glutCreateWindow("Stable Fluids 3D");

	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutSpecialFunc(SpecialKeys);
	glutKeyboardFunc(KeyBoard);
	glutIdleFunc(Idle);

	glutMainLoop();

}

*/
