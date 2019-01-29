
#include "Common.h"
#include <unistd.h>
#include "utils_.h"
#include "utils_.h"

#include <string.h>
CSkyDome 			Sky;
CCamera 			Camera;
CFrustum 			Frustum;
GLText				Text;
VolumetricClouds 	VolumeClouds;

SDL_Surface*	MainWindow;

int 	ScreenWidth = 800;
int 	ScreenHeight = 600;
int 	ScreenBpp = 32;
int 	ScreenDepth = 24;

float 	CloudPlaneSize = 400.0f;
float 	CloudPlaneHeight = 100.0f;
float 	SkyRadius = 200.0f;

bool	fullscreen = false;

#ifdef WIN32
double winGetFrequency()
{
	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);
	return (double)(freq.QuadPart);
}
#endif

double sysGetTime()
{
#ifdef WIN32
	static double frequency = winGetFrequency();
	
	LARGE_INTEGER current;
	QueryPerformanceCounter(&current);
	return (double)(current.QuadPart) / frequency;
#else
//	struct timeval t;
//	gettimeofday(&t, NULL);
	
//	return (double)(t.tv_sec) + (VTime)(t.tv_usec) / 1000000.0;
    return 0;
#endif
}

int init_game()
{
    char tga_path[100];
    strcpy(tga_path,(get_cwd()+"/../../data/Cloud_Render/Harris/font.tga").c_str());
	if (Text.LoadFont(
            tga_path
	        , ScreenWidth, ScreenHeight)) return 1;

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	Text.TextOut(ScreenWidth / 2 - 100, ScreenHeight / 2, "PLEASE WAIT, LOADING...");

	SDL_GL_SwapBuffers();

	Camera.Create(45.0f, ScreenWidth, ScreenHeight, 4096.0f);
	if (Sky.Initialize(SkyRadius, 24, 32, false)) return 1;	
    if (VolumeClouds.Create(1, CloudPlaneSize / 2, CloudPlaneHeight)) return 1;
	//hack * 10 :)
	//Camera.PositionCamera(0, 0, 0, Sky.GetSunVector().x * 10, Sky.GetSunVector().y * 10, Sky.GetSunVector().z * 10, 0, 1, 0);
	Camera.PositionCamera(16, CloudPlaneHeight + 16, 100, 16, CloudPlaneHeight + 16, 0, 0, 1, 0);

	SDL_ShowCursor(SDL_DISABLE);

	unsigned e;		
	if ( (e = glGetError()) )
	{
		cout << "OpenGL Error: " << gluErrorString(e) << endl;
		return 1;
	}
			
	return 0;
}

int update_frame()
{		
	static float FPS = 0;
	static double LastTime;

	double CurrentTime = sysGetTime();	
	FPS = 1 / (CurrentTime - LastTime);
	LastTime = CurrentTime;

	//the update part
	glLoadIdentity();				
	//Camera.Update();
	Camera.Look();
	Frustum.CalculateFrustum();	

	VolumeClouds.Update(Sky.GetSunVector() * SkyRadius, Camera.GetCameraPosition());
	
	//the rendering part
	glClearColor(0.1f, 0.44f, 0.62f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glColor3f(1.0f, 1.0f, 1.0f);


		
	//Sky.Render(Camera.GetCameraPosition());
	VolumeClouds.Render(Camera.GetCameraPosition(), Sky.GetSunVector()* SkyRadius );//
	
	int spr, imp;
	VolumeClouds.GetInfo(&spr, &imp);
	Text.TextOut(0, 0, "%3.1f PARTICLES: %d IMPOSTORS %d", FPS, spr, imp);

	unsigned e;		
	if ( (e = glGetError()) )
	{
		cout << "OpenGL Error: " << gluErrorString(e) << endl;
		return 1;
	}
	
	SDL_GL_SwapBuffers();
	
	return 0;
}

void terminate_game()
{		
	Sky.Destroy();
	Text.Destroy();
}

int init_opengl()
{	
	glRenderMode(GL_RENDER);
	glClearColor(0.0f, 0.8f, 1.0f, 0.0f);
	glDepthFunc(GL_LEQUAL);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	return 0;
}

int set_display_mode(int width, int height, int bpp, bool fullscreen)
{
	unsigned int flags = SDL_OPENGL | SDL_HWSURFACE;  // | SDL_RESIZABLE;

	if (fullscreen) 	flags |= SDL_FULLSCREEN;
		    	
	SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );		// tell SDL that the GL drawing is going to be double buffered
    SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, ScreenDepth);		// size of depth buffer
    SDL_GL_SetAttribute( SDL_GL_STENCIL_SIZE, 0);		// we aren't going to use the stencil buffer
    SDL_GL_SetAttribute( SDL_GL_ACCUM_RED_SIZE, 0);		// this and the next three lines set the bits allocated per pixel -
    SDL_GL_SetAttribute( SDL_GL_ACCUM_GREEN_SIZE, 0);	// - for the accumulation buffer to 0
    SDL_GL_SetAttribute( SDL_GL_ACCUM_BLUE_SIZE, 0);
    SDL_GL_SetAttribute( SDL_GL_ACCUM_ALPHA_SIZE, 0);
	SDL_GL_SetAttribute( SDL_GL_RED_SIZE, 8);
	SDL_GL_SetAttribute( SDL_GL_BLUE_SIZE, 8);
	SDL_GL_SetAttribute( SDL_GL_GREEN_SIZE, 8);
	SDL_GL_SetAttribute( SDL_GL_ALPHA_SIZE, 8);

	if (MainWindow)
		SDL_FreeSurface(MainWindow);

	MainWindow = SDL_SetVideoMode(width, height, bpp, flags);

	if (!MainWindow)
	{
		cout << "INIT: Cannot initialise video mode or create window in " << width << "x" << height << "x" << bpp << " DepthSize: " << ScreenDepth << " " << endl;
		return 1;
	}
	
	cout << "INIT: Video mode set to " << width << "x" << height << "x" << bpp << " DepthSize: " << ScreenDepth << " " << endl;
	
	return 0;
}

#ifdef WIN32
int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
#else
int main(int argc, char *argv[])
#endif
{
//	char cwnd[100];
//	cout << get_cwd() <<"\n"
//	<< get_current_dir_name()<<"\n"
//
//	;


	cout << "Initializing Skydome demo" << endl;
	if (SDL_Init(SDL_INIT_VIDEO) <  0) 
	{
	    cout << "Could not initialize Sdl\n";
		cout << "Could not initialize SDL:" << SDL_GetError() << endl;
		SDL_Quit();
		return 1;
	} 
	
	cout << "SDL Initialization successful" << endl;

	if (set_display_mode(ScreenWidth, ScreenHeight, ScreenBpp, fullscreen))
	{
		SDL_Quit();
		return 1;
	}
		
	init_opengl();
	
	ilInit();
	
	if (init_game())
	{
		ilShutDown();
		SDL_Quit();
		return 1;
	}
	
	bool running = true;
	
	SDL_Event event;
	
	while (running)
	{
		while (SDL_PollEvent(&event))
		{
			switch (event.type)
			{
				case SDL_VIDEORESIZE:
					ScreenWidth = event.resize.w;
					ScreenHeight = event.resize.h;					
					set_display_mode(ScreenWidth, ScreenHeight, ScreenBpp, false);		
					break;
					
				case SDL_QUIT:
					running = false;
					break;
				case SDL_KEYDOWN:
					switch (event.key.keysym.sym)
					{
						case SDLK_ESCAPE:
							running = false;
							break;
						case SDLK_u:
							glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);//�������ģʽ
							break;
						case SDLK_i:
							glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
							break;
						case SDLK_w:
							Camera.MoveCamera(10);
							break;
						case SDLK_s:
							Camera.MoveCamera(-10);
							break;
						case SDLK_a:
							Camera.RotateAroundPoint(Camera.m_vView,0.2,0,1,0);
							break;
						case SDLK_d:
							Camera.RotateAroundPoint(Camera.m_vView,0.2,0,-1,0);
							break;
						case SDLK_LEFT:
							Camera.m_vStrafe.x=1;
							Camera.StrafeCamera(1);
							Camera.m_vStrafe = Vector3();
							break;
						case SDLK_RIGHT:
							Camera.m_vStrafe.x=1;
							Camera.StrafeCamera(-1);
							Camera.m_vStrafe = Vector3();
							break;
						case SDLK_UP:
							Camera.m_vStrafe.y=1;
							Camera.StrafeCamera(1);
							Camera.m_vStrafe = Vector3();
							break;
						case SDLK_DOWN:
							Camera.m_vStrafe.y=1;
							Camera.StrafeCamera(-1);
							Camera.m_vStrafe = Vector3();
							break;

						default:
							break;

					}
					break;

			}			
		}
		//Camera.MoveCamera(-0.1);
		//Camera.RotateView(0.01,0,0,1);
		//Camera.RotateAroundPoint(Camera.m_vView,0.01,0,0,0);
//		Camera.m_vStrafe.y=1;
//		Camera.StrafeCamera(0.01);
		//Camera.Update();
		int r = update_frame();
		if (r){
            break;
		}
	}	
	
	terminate_game();
	
	ilShutDown();
	
	SDL_Quit();
	
	return 0;
}
