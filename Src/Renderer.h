#ifndef RENDERER_H
#define RENDERER_H

#include "Headers.h"

class GUI;
class CameraManager;
class FrameBufferObjectMultisample;
class FrameBufferObject;
class Shader;
class Scene;
class VertexBufferObject;
class Texture;
class Scanner;

class Renderer
{
public:
    Renderer(Scene *scene, CameraManager *camManager, GUI *gui);
    ~Renderer();

    void init();
    void render(Transform &trans);
    
    void resize(int width, int height);
    void toggleBGColor();
    void togglePolygonMode();

private:
    void renderScene(const Transform &trans);

public:
	Scanner *m_scanner;

private:
    GUI *m_gui;
    CameraManager *m_cameraManager;
	Scene *m_scene;

    int m_width;
    int m_height;

    vec4 m_bgColor; 

    const GLuint m_samples;
    GLuint m_bgMode;

	int m_bufferWidth;
	int m_bufferHeight;
	
};

#endif

 