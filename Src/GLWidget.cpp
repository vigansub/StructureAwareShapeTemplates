#include "GLWidget.h"
#include "Shader.h"
#include "Camera.h"
#include "Renderer.h"
#include "Scene.h"
#include "GUI.h"
#include "CameraManager.h"
#include "Light.h"
#include "Grammar.h"
#include "Scanner.h"
#include "postoffice.h"


GLWidget::GLWidget(QGLContext *context, int width, int height)
: QGLWidget(context),
  m_leftButton(false),
  m_rightButton(false),
  m_width(width),
  m_height(height),
  m_ctrlPressed(false),
  m_altPressed(false),
  m_shiftPressed(false),
  m_noOtherKey(true),
  m_renderOffline(false),
  m_oldTime(0),
  m_frameNr(0)
{
    this->resize(m_width, m_height);
    this->setMouseTracking(true);
}

GLWidget::~GLWidget()
{
	delete m_scene;
	delete m_cameraManager;
	delete m_gui;
	delete m_renderer;
}

void GLWidget::initializeGL() 
{
    glewInit();    

	initParams();
    initShaders();

	m_cameraManager = new CameraManager();
    m_scene = new Scene(this, m_cameraManager);
	m_gui = new GUI(m_cameraManager, m_scene);
    m_renderer = new Renderer(m_scene, m_cameraManager, m_gui);

    m_renderer->toggleBGColor();
    m_renderer->toggleBGColor();
    loop(params::inst()->gridRenderMode, 0, 4);
    loop(params::inst()->gridRenderMode, 0, 4);
    loop(params::inst()->gridRenderMode, 0, 4);

	m_gui->toggleMode();

    glEnable(GL_DEPTH_TEST);     
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glCullFace(GL_BACK);
    glFrontFace(GL_CCW);    

    glEnable(GL_SAMPLE_SHADING);
    glMinSampleShading(params::inst()->sampleShading);


    Postoffice::setGLWidget(this);
    Postoffice::setRenderer(m_renderer);
}

void GLWidget::initParams()
{
    GlobalObjectParams *p = params::inst();

    p->camPos              = vec3(0.0f, 0.0f, 0.0f);
	p->blur                = vec2(8.0f, 0.64f);
//    p->shadowMapSize       = vec2(2048, 2048);
    p->shadowMapSize       = vec2(2048, 2048);
    p->applyShadow         = true;
    p->gridRenderMode      = 0;
    p->shadowIntensity     = 0.6f;
    p->sampleShading       = 1.0f;
    p->polygonMode         = 0;
    p->activeLight         = 0;
    
    p->windowSize          = vec2(m_width, m_height);
    
    p->renderMesh          = false;
    p->renderObjects       = true;
    p->renderTextures      = false;
    p->renderWireframe     = false;
    p->renderNormals       = false;
    p->renderMisc          = false;
	p->renderTarget        = false;
	p->renderCurrent       = false;
    p->renderPointCloud    = true;
    p->renderGenerations   = false;
    p->isPartial           = true;
    p->renderBoxes         = true;
    p->renderConnectors    = true;
    p->boxTransparency     = 1.0f;

    p->clipPlaneGround     = vec4(0.0f, -1.0f, 0.0f, 4.0f);
    p->ncp                 = 0.0f;
    p->fcp                 = 0.0f;
    p->fov                 = 0.0f;
    p->lightIntensity      = 1.2f;
    p->renderHeightOffset = 1.0f;
    
    p->polygonOffsetUnits  = 1.0f;
    p->polygonOffsetFactor = 0.5f;
    p->depthRangeMax       = 1.0f;
    p->depthRangeMin       = 0.0f;
    
    p->nrVertices          = 0;
    p->nrActiveVertices    = 0;

    p->weightProjection = 0.3f;
    p->weightBoundingBox = 0.8f;
    p->weightNegative = 1.0f;
	p->weightDistanceField = 0.0f;
    p->weightGrammar = 0.0f;
	p->weightDisentangle = 0.4f;
    p->weightSkin = 0.8f;

    p->optIterations = 700;
    p->optParticles = 100;
    p->optProbability = 0.8f;
    p->optScale = 0.7f;
    p->optStrategy = 0;
    p->optTechnique = 0;
    p->optMean = 0.2f;
    p->optStd = 0.005f;
    p->optW = 0.7f;
    p->optphip = 2.05f;
    p->optphig = 2.05f;
    p->optPopulation = 200;
    p->optCrossover = 0.7f;
    p->optMutation = 0.001f;
    p->optBits = 10;

	p->error = 0.0;
}

void GLWidget::initShaders()
{
#ifdef WIN32
    shaders::inst()->default = new Shader("Shader/Default.vert.glsl", "Shader/Default.frag.glsl");
    shaders::inst()->default->bindAttribLocations();
#else
    shaders::inst()->default_ = new Shader("./../_main/Shader/Default.vert.glsl", "./../_main/Shader/Default.frag.glsl");
    shaders::inst()->default_->bindAttribLocations();
#endif

#ifdef WIN32
    shaders::inst()->defaultLight = new Shader("Shader/DefaultLight.vert.glsl", "Shader/DefaultLight.frag.glsl");
    shaders::inst()->defaultLight->bindAttribLocations();

    shaders::inst()->defaultDepth = new Shader("Shader/Default.vert.glsl", "Shader/DefaultDepth.frag.glsl");
    shaders::inst()->defaultDepth->bindAttribLocations();

    shaders::inst()->blur = new Shader("Shader/Blur.vert.glsl", "Shader/Blur.frag.glsl");
    shaders::inst()->blur->bindAttribLocations();

    shaders::inst()->grid = new Shader("Shader/NiceGrid.vert.glsl", "Shader/NiceGrid.frag.glsl");
    shaders::inst()->grid->bindAttribLocations();

    shaders::inst()->object = new Shader("Shader/Object.vert.glsl", "Shader/Object.frag.glsl");
    shaders::inst()->object->bindAttribLocations();

    shaders::inst()->objectDepth = new Shader("Shader/ObjectDepth.vert.glsl", "Shader/ObjectDepth.frag.glsl");
    shaders::inst()->objectDepth->bindAttribLocations();

    shaders::inst()->objectLines = new Shader("Shader/ObjectLines.vert.glsl", "Shader/ObjectLines.frag.glsl");
    shaders::inst()->objectLines->bindAttribLocations();

	shaders::inst()->gui = new Shader("Shader/GUI.vert.glsl", "Shader/GUI.frag.glsl");
	shaders::inst()->gui->bindAttribLocations();

	shaders::inst()->noise = new Shader("Shader/Noise.vert.glsl", "Shader/Noise.frag.glsl");
	shaders::inst()->noise->bindAttribLocations();

	shaders::inst()->cookTorrance = new Shader("Shader/CookTorrance.vert.glsl", "Shader/CookTorrance.frag.glsl");
	shaders::inst()->cookTorrance->bindAttribLocations();

	shaders::inst()->sphericalHarmonic = new Shader("Shader/SphericalHarmonic.vert.glsl", "Shader/SphericalHarmonic.frag.glsl");
	shaders::inst()->sphericalHarmonic->bindAttribLocations();

    shaders::inst()->tessellation = new Shader();
    shaders::inst()->tessellation->attachVertexShader("Shader/TessInterp.vert.glsl");
    shaders::inst()->tessellation->attachControlShader("Shader/TessInterp.cont.glsl");
    shaders::inst()->tessellation->attachEvaluationShader("Shader/TessInterp.eval.glsl");
    shaders::inst()->tessellation->attachGeometryShader("Shader/TessInterp.geom.glsl");
    shaders::inst()->tessellation->attachFragmentShader("Shader/TessInterp.frag.glsl");
    shaders::inst()->tessellation->bindAttribLocations();

	shaders::inst()->box = new Shader("Shader/Box.vert.glsl", "Shader/Box.frag.glsl");
	shaders::inst()->box->bindAttribLocations();

	shaders::inst()->scan = new Shader("Shader/Scan.vert.glsl", "Shader/Scan.frag.glsl");
	shaders::inst()->scan->bindAttribLocations();

	shaders::inst()->points = new Shader("Shader/Points.vert.glsl", "Shader/Points.frag.glsl");
	shaders::inst()->points->bindAttribLocations();
#else
    shaders::inst()->defaultLight = new Shader("./../_main/Shader/DefaultLight.vert.glsl", "./../_main/Shader/DefaultLight.frag.glsl");
    shaders::inst()->defaultLight->bindAttribLocations();

    shaders::inst()->defaultDepth = new Shader("./../_main/Shader/Default.vert.glsl", "./../_main/Shader/DefaultDepth.frag.glsl");
    shaders::inst()->defaultDepth->bindAttribLocations();

    shaders::inst()->blur = new Shader("./../_main/Shader/Blur.vert.glsl", "./../_main/Shader/Blur.frag.glsl");
    shaders::inst()->blur->bindAttribLocations();

    shaders::inst()->grid = new Shader("./../_main/Shader/NiceGrid.vert.glsl", "./../_main/Shader/NiceGrid.frag.glsl");
    shaders::inst()->grid->bindAttribLocations();

    shaders::inst()->object = new Shader("./../_main/Shader/Object.vert.glsl", "./../_main/Shader/Object.frag.glsl");
    shaders::inst()->object->bindAttribLocations();

    shaders::inst()->objectDepth = new Shader("./../_main/Shader/ObjectDepth.vert.glsl", "./../_main/Shader/ObjectDepth.frag.glsl");
    shaders::inst()->objectDepth->bindAttribLocations();

    shaders::inst()->objectLines = new Shader("./../_main/Shader/ObjectLines.vert.glsl", "./../_main/Shader/ObjectLines.frag.glsl");
    shaders::inst()->objectLines->bindAttribLocations();

    shaders::inst()->gui = new Shader("./../_main/Shader/GUI.vert.glsl", "./../_main/Shader/GUI.frag.glsl");
    shaders::inst()->gui->bindAttribLocations();

    shaders::inst()->noise = new Shader("./../_main/Shader/Noise.vert.glsl", "./../_main/Shader/Noise.frag.glsl");
    shaders::inst()->noise->bindAttribLocations();

    shaders::inst()->cookTorrance = new Shader("./../_main/Shader/CookTorrance.vert.glsl", "./../_main/Shader/CookTorrance.frag.glsl");
    shaders::inst()->cookTorrance->bindAttribLocations();

    shaders::inst()->sphericalHarmonic = new Shader("./../_main/Shader/SphericalHarmonic.vert.glsl", "./../_main/Shader/SphericalHarmonic.frag.glsl");
    shaders::inst()->sphericalHarmonic->bindAttribLocations();

    shaders::inst()->tessellation = new Shader();
    shaders::inst()->tessellation->attachVertexShader("./../_main/Shader/TessInterp.vert.glsl");
    shaders::inst()->tessellation->attachControlShader("./../_main/Shader/TessInterp.cont.glsl");
    shaders::inst()->tessellation->attachEvaluationShader("./../_main/Shader/TessInterp.eval.glsl");
    shaders::inst()->tessellation->attachGeometryShader("./../_main/Shader/TessInterp.geom.glsl");
    shaders::inst()->tessellation->attachFragmentShader("./../_main/Shader/TessInterp.frag.glsl");
    shaders::inst()->tessellation->bindAttribLocations();

    shaders::inst()->box = new Shader("./../_main/Shader/Box.vert.glsl", "./../_main/Shader/Box.frag.glsl");
    shaders::inst()->box->bindAttribLocations();

    shaders::inst()->scan = new Shader("./../_main/Shader/Scan.vert.glsl", "./../_main/Shader/Scan.frag.glsl");
    shaders::inst()->scan->bindAttribLocations();

    shaders::inst()->points = new Shader("./../_main/Shader/Points.vert.glsl", "./../_main/Shader/Points.frag.glsl");
    shaders::inst()->points->bindAttribLocations();
#endif
}

void GLWidget::paintGL()
{   
	params::inst()->windowSize.x = this->width();
	params::inst()->windowSize.y = this->height();
    params::inst()->lights = m_scene->m_lights;
    params::inst()->nrActiveVertices = 0;

    m_renderer->render(m_trans);

    if (m_renderOffline)
    {
        saveFrameBuffer(this, m_frameNr);
        m_frameNr++;
    }

    DWORD time = GetTickCount();
    float delta = time - m_oldTime;
    if(delta > 25)
    {
        m_scene->update(1.0f/delta);
        m_oldTime = time;

        m_cameraManager->currentPerspective(m_trans);
        m_cameraManager->currentCamParams();
    }

    this->update();
}

void GLWidget::resizeGL(int w, int h)
{
	m_width = w;
	m_height = h;

	m_renderer->resize(w, h);
	m_cameraManager->resize(w, h);
    m_gui->resize(w, h);
}

void GLWidget::wheelEvent(QWheelEvent *event)
{
    if(!m_altPressed && !m_ctrlPressed && !m_rightButton)
    {
	    m_cameraManager->onMouseWheel(event->delta());
    }

    if (m_altPressed && !m_ctrlPressed)
    {
        vec4 lPos = m_scene->m_lights[params::inst()->activeLight]->position();
	    float delta = lPos.y * 0.1;

        if (event->delta() > 0) 
			m_scene->m_lights[params::inst()->activeLight]->setPosition(vec3(lPos.x, lPos.y+delta, lPos.z));
        else 
            m_scene->m_lights[params::inst()->activeLight]->setPosition(vec3(lPos.x, lPos.y-delta, lPos.z));
    }

	if (m_rightButton) 
	{
		if (event->delta() > 0)
			m_cameraManager->changeRotHeight(-0.2f);
		else
			m_cameraManager->changeRotHeight(0.2f);
	}

    if (m_altPressed && m_ctrlPressed)
    {
        float lPos = m_scene->m_ballEraserPos.y;
        float delta = lPos * 0.1;

        if (event->delta() > 0)
            m_scene->erasePoints(0.0, 0.0, delta);
        else
            m_scene->erasePoints(0.0, 0.0, -delta);
    }

    event->accept();    
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    QPoint newMouse(event->x(), event->y());
	QPoint diff = newMouse - m_mouse;

    m_gui->onMouseMove(event->x(), event->x());

    if(!m_altPressed && !m_ctrlPressed)
    {
	    if(m_leftButton)
	    {            		
		    m_cameraManager->onMouseMove(diff.x(), diff.y(), 0);
	    }
	     else if(m_rightButton)
	    {
		    m_cameraManager->onMouseMove(diff.x(), diff.y(), 1);
	    }
    }

    if(m_leftButton && m_altPressed && !m_ctrlPressed)
    {
        if(m_ctrlPressed)
            m_scene->m_lights[params::inst()->activeLight]->recordPath(true);
        else
            m_scene->m_lights[params::inst()->activeLight]->recordPath(false);

		m_scene->m_lights[params::inst()->activeLight]->move(m_cameraManager, diff.x()*0.1, diff.y()*0.1);
    }

    if(m_leftButton && m_ctrlPressed)
    {
        m_scene->move(m_trans, diff.x(), diff.y());
    }

    if (m_leftButton && m_altPressed && m_ctrlPressed)
    {
        m_scene->erasePoints(diff.x()*0.01, diff.y()*0.01, 0.0f);
    }

    m_mouse.setX(event->x());
    m_mouse.setY(event->y());    
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    if(!m_gui->onMouseClick(event->x(), event->y()))
    {
        if(event->button() == Qt::LeftButton)
            m_leftButton = true; 

        if(event->button() == Qt::RightButton)
            m_rightButton = true;
    }

    if(m_ctrlPressed)
    {
        m_scene->resetSelection();
        m_scene->select(m_trans, this->width(), this->height(), event->x(), event->y());
    }
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    m_gui->onMouseRelease();

    if(event->button() == Qt::LeftButton)
        m_leftButton = false;   

    if(event->button() == Qt::RightButton)
        m_rightButton = false;
}

void GLWidget::keyPressEvent(QKeyEvent *event)
{   
	m_cameraManager->onKeyPress(event->key());

    switch (event->key()) 
    {
        case Qt::Key_Tab:
            if (m_noOtherKey)
                m_gui->toggleMode();
            break;
		case Qt::Key_PageUp: 
			m_cameraManager->increaseSpeed();
            break;
        case Qt::Key_PageDown:
			m_cameraManager->decreaseSpeed();
            break;
        case Qt::Key_Left:
            break;
        case Qt::Key_Right:
            break;
		case Qt::Key_Space:
			loop(params::inst()->polygonMode, 0, 1, 1);
			break;
        case Qt::Key_Plus:
            break;
        case Qt::Key_Minus:
            break;
        case Qt::Key_Enter:
        case Qt::Key_Return:
            break;

        case Qt::Key_A:
            //Camera Strafe Left
            break;
        case Qt::Key_B:
            m_scene->m_grammar->softCleanUp();
            m_scene->m_grammar->derive();
            break;
        case Qt::Key_C:
            m_cameraManager->lockCurCamera();
            break;
        case Qt::Key_D:
            //Camera Strafe Right
            break;
        case Qt::Key_E:
            m_scene->newGenerate();
            break;
        case Qt::Key_F:
            break;
        case Qt::Key_G:
            params::inst()->renderMisc = !params::inst()->renderMisc;
            break;
        case Qt::Key_H:
            m_renderer->m_scanner->updateBuffer();
            break;
        case Qt::Key_I:   
            m_cameraManager->toggleInterpolation();
            break;
        case Qt::Key_J:
			m_scene->testIndependentVariables();
            break;
        case Qt::Key_K:
            loop(params::inst()->activeLight, 0, (int)params::inst()->lights.size()-1);
            break;
        case Qt::Key_L:
            m_scene->m_lights[params::inst()->activeLight]->toggleMode();
            break;
        case Qt::Key_M:		
			m_scene->m_grammar->printParameters();
            break;
        case Qt::Key_N:
            stats.print();
            break;
        case Qt::Key_O:
			m_scene->optimize();
            break;
        case Qt::Key_P:
            loop(params::inst()->gridRenderMode, 0, 4);
            break;
        case Qt::Key_Q:
            //Shape Merging
            m_scene->mergeShapes("chair_01", "chair_04");
            break;
        case Qt::Key_U:
            params::inst()->applyShadow = !params::inst()->applyShadow;
            break;
        case Qt::Key_R:
            //Generating from Obtained mean and Std
            m_scene->renderOptimal();
            //m_scene->newShapeResize();
            break;
        case Qt::Key_S:
            //Generate Error Matrix Parallel
            //m_scene->generateErrorMatrix_Parallel();
            m_scene->compareMinhyuk();
            m_renderer->m_scanner->updateBuffer();
            break;
        case Qt::Key_T:
            //Generate Error Matrix
            m_scene->generateErrorMatrix();
            break;
        case Qt::Key_V:
            m_scene->m_grammar->update();
            break;
        case Qt::Key_W:
            //Camera Forward
            m_scene->obtain_gram();
            break;
        case Qt::Key_X:
			m_renderer->m_scanner->storeView();
            //m_renderer->m_scanner->m_scanning = true;
            //m_renderer->m_scanner->autoScan();
            break;
        case Qt::Key_Y:
			m_renderer->m_scanner->clearBuffer();
            break;
        case Qt::Key_Z:
            //m_scene->newGenerate();
            m_scene->savePartialResults();
            //m_scene->skinGenerate();
            //m_scene->generateNewDatabase();
            //m_scene->test_MultipleShapes();
            //m_scene->generateContinuousSkinning();
            //m_scene->generateContinuousHeatMap();
            break;

        case Qt::Key_F1:
            this->setWindowState(this->windowState() ^ Qt::WindowFullScreen); 
            break;
        case Qt::Key_F2:
			m_cameraManager->toggleCam();
            break;        
        case Qt::Key_F3:            
            saveFrameBuffer(this);
            break;        
        case Qt::Key_F4:
			m_renderer->toggleBGColor();
            break;
        case Qt::Key_F5:
            if(m_noOtherKey)
			    m_gui->toggleMode();

            if(m_ctrlPressed)
                m_cameraManager->toggleFrameset();
            break;
       case Qt::Key_F6:			
            if(m_ctrlPressed)
                m_cameraManager->addFrame();
            break;
       case Qt::Key_F7:			
            if(m_ctrlPressed)
                m_cameraManager->clearFrameset();
            break;
       case Qt::Key_F8:			
            if(m_ctrlPressed)
                m_cameraManager->saveFrameset();
            break;
       case Qt::Key_F9:
            m_renderOffline = !m_renderOffline;
            m_frameNr = 0;
            break;
       case Qt::Key_F10:			
            break;
       case Qt::Key_F11:	
           //m_scene->writeBoxes_grammars();
           m_scene->writeBoxes();
          //m_scene->saveObjboxes("test/test_boxes.obj");
            break;
       case Qt::Key_F12:	
           //m_scene->saveAllImages();
           //m_scene->obtainBestFitShape();
           //m_scene->sceneCompletion();
           m_scene->findClosestShape();
            break;

		case Qt::Key_Control:
			m_ctrlPressed = true;
            m_noOtherKey = false;
			break;
		case Qt::Key_Alt:
			m_altPressed = true;
            m_noOtherKey = false;
			break;
		case Qt::Key_Shift:
			m_shiftPressed = true;
            m_noOtherKey = false;
			break;
        case Qt::Key_Escape:            
            exit(0);
            break;
    }
}

void GLWidget::keyReleaseEvent(QKeyEvent *event)
{
	m_cameraManager->onKeyRelease(event->key());

    switch (event->key()) 
    {
		case Qt::Key_Control:
			m_ctrlPressed = false;
            m_noOtherKey = true;
			break;
		case Qt::Key_Alt:
			m_altPressed = false;
            m_noOtherKey = true;
			break;
		case Qt::Key_Shift:
			m_shiftPressed = false;
            m_noOtherKey = true;
			break;
		default:
            break;
    }
}
