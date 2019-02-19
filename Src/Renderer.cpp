#include "Renderer.h"
#include "GUI.h"
#include "CameraManager.h"
#include "FrameBufferObject.h"
#include "Shader.h"
#include "VertexBufferObject.h"
#include "Scene.h"
#include "Mesh.h"
#include "Light.h"
#include "Texture.h"
#include "Scanner.h"

Renderer::Renderer(Scene *scene, CameraManager *camManager, GUI *gui)
: m_scene(scene),
  m_gui(gui),
  m_cameraManager(camManager),
  m_bgColor(0.1f, 0.1f, 0.1f, 1.0f),
  m_width(0),
  m_height(0),
  m_samples(16),
  m_bgMode(0),
  m_bufferWidth(0),
  m_bufferHeight(0)
{
    init();

	m_scanner = new Scanner(m_scene);
}

Renderer::~Renderer()
{
}

void Renderer::init()
{
}

void Renderer::render(Transform &trans)
{
    if(params::inst()->applyShadow)
    {
        trans.lightViews.clear();
        trans.lightViews.resize(m_scene->m_lights.size());

        for(int i=0; i<m_scene->m_lights.size(); ++i)
        {
            m_scene->m_lights[i]->setIntensity(params::inst()->lightIntensity);
            m_scene->m_lights[i]->setDirection(m_scene->m_lights[i]->position());
            m_scene->m_lights[i]->renderLightView(trans.lightViews[i]); 
        }
    }    

	m_scanner->renderScan(trans, true);
    renderScene(trans);	

    if(params::inst()->renderTextures)
    {
        for(int i=0; i<m_scene->m_lights.size(); ++i)
        {
            renderTexture(m_scene->m_lights[i]->m_fboLight->texAttachment(GL_COLOR_ATTACHMENT0), 220+i*210, m_height-200, 200, 200);
        }
    }

    m_gui->render();
}

void Renderer::renderScene(const Transform &trans)
{
    glViewport(0, 0, m_width, m_height);
    glClearColor(m_bgColor.x, m_bgColor.y, m_bgColor.z, m_bgColor.w);    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);     

    glEnable(GL_MULTISAMPLE);        

    if (params::inst()->polygonMode == 1)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    //glEnable(GL_POLYGON_OFFSET_FILL);
    //glEnable(GL_POLYGON_OFFSET_POINT);
    //glEnable(GL_POLYGON_OFFSET_LINE);
    //glDepthFunc(GL_LEQUAL);
    //glDepthRange(0.1, 1.0);
    //glPolygonOffset(-1.0, -10.5);
		    
    m_scene->renderObjects(trans);
    m_scene->renderWorld(trans);  

	if (params::inst()->renderPointCloud)
		m_scanner->renderPointCloud(trans);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void Renderer::resize(int width, int height)
{
    m_width = width;
    m_height = height;  

	m_bufferWidth = m_width / 8;
	m_bufferHeight = m_height / 8;

	m_scanner->resize(m_bufferWidth, m_bufferHeight);
}

void Renderer::toggleBGColor()
{
    m_bgMode ++;
    if(m_bgMode > 2)
        m_bgMode = 0;

    if(m_bgMode == 0)
	{
        m_bgColor = vec4(0.1f, 0.1f, 0.1f, 1.0f);
		m_gui->setFontColor(vec4(0.9f, 0.9f, 0.9f, 1.0f));
	}
   
    if(m_bgMode == 1)
	{
        m_bgColor = vec4(0.5f, 0.5f, 0.5f, 1.0f);
		m_gui->setFontColor(vec4(0.0f, 0.0f, 0.0f, 1.0f));
	}

    if(m_bgMode == 2)
	{
        m_bgColor = vec4(1.0f, 1.0f, 1.0f, 1.0f);
		m_gui->setFontColor(vec4(0.0f, 0.0f, 0.0f, 1.0f));
	}
}

