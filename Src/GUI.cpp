#include "GUI.h"
#include "CameraManager.h"
#include "Mesh.h"
#include "Shader.h"
#include "VertexBufferObject.h"
#include "Scene.h"
#include "Light.h"

GUI::GUI(CameraManager *cameraManager, Scene *m_scene)
: m_width(0),
  m_height(0),
  m_fps(0),
  m_fpsCount(0),
  m_mode(0),
  m_oldTime(0),
  m_cameraManager(cameraManager),
  m_scene(m_scene),
  m_fontColor(0.9f, 0.9f, 0.9f, 1.0), 
  m_guiAll(5),
  m_tabNames(5),
  m_maxMode(0), 
  m_curTab(-1)
{
	initShortcuts();
	initUIElements();

	connect(&m_timer, SIGNAL(timeout()), this, SLOT(onTimer()));
	m_timer.start(200);

    m_hpTimer.reset();

	m_vboQuad = Mesh::quad(0, 0, 1, 1, vec4(0.0f, 0.0f, 0.0f, 0.7f));
}

GUI::~GUI()
{
	delete m_vboQuad;
}

void GUI::initShortcuts()
{
	m_shortcuts.push_back(QString("Toggle Fullscreen (F1)"));
	m_shortcuts.push_back(QString("Toggle Cameras (F2)"));
	m_shortcuts.push_back(QString("Screen Shot (F3)"));
	m_shortcuts.push_back(QString("Toggle Background Color (F4)"));
    m_shortcuts.push_back(QString("Toggle Shadow (U)"));
    m_shortcuts.push_back(QString("Cam Speed (PgUp/PgDown)"));
    m_shortcuts.push_back(QString("Interpolate Cam (I)"));
}

void GUI::initUIElements()
{
    int nrTabs = 5;
    m_maxMode = 1 + (nrTabs-1);

    initTab1(TAB1);
    initTab2(TAB2);
    initTab3(TAB3);
    initTab4(TAB4);
    initTab5(TAB5);
}

void GUI::initTab1(Slots slot)
{
    m_tabNames[slot] = "Grammar";

    GlobalObjectParams *p = params::inst();

    int y = 20;
    int x = 10;

	setupFloatSlider(slot, x, y += 45, 0.0f, 1.0f, QString("Weight Projection: "), &p->weightProjection);
	setupFloatSlider(slot, x, y += 45, 0.0f, 1.0f, QString("Weight BB: "), &p->weightBoundingBox);
	setupFloatSlider(slot, x, y += 45, 0.0f, 1.0f, QString("Weight Negative: "), &p->weightNegative);
	setupFloatSlider(slot, x, y += 45, 0.0f, 1.0f, QString("Weight Skinning: "), &p->weightSkin);
	setupFloatSlider(slot, x, y += 45, 0.0f, 1.0f, QString("Weight Disentanglement: "), &p->weightDisentangle);
	setupFloatSlider(slot, x, y += 45, 0.0f, 1.0f, QString("Weight DF: "), &p->weightDistanceField);


	setupIntSlider(slot, x, y += 65, 0, 10000, QString("Iterations: "), &p->optIterations);
	setupFloatSlider(slot, x, y += 45, 0.0f, 1.0f, QString("Scale: "), &p->optScale);
	setupFloatSlider(slot, x, y += 45, 0.0f, 1.0f, QString("Probability: "), &p->optProbability);
    setupFloatSlider(slot, x, y += 45, 0.0f, 2.0f, QString("CMA-ES Mean: "), &p->optMean);
    setupFloatSlider(slot, x, y += 45, 0.0f, 0.5f, QString("CMA-ES Standard Deviation: "), &p->optStd);
    setupIntSlider(slot, x, y += 65, 0, 1000, QString("Swarm Particles: "), &p->optParticles);
    setupFloatSlider(slot, x, y += 45, -5.0f, 5.0f, QString("Swarm Omega: "), &p->optW);
    setupFloatSlider(slot, x, y += 45, -5.0f, 5.0f, QString("Swarm Local Variant: "), &p->optphip);
    setupFloatSlider(slot, x, y += 45, -5.0f, 5.0f, QString("Swarm Global Variant: "), &p->optphig);
    setupIntSlider(slot, x, y += 45, 0, 1000, QString("Population Size: "), &p->optPopulation);
    setupFloatSlider(slot, x, y += 45, 0.0f, 1.0f, QString("Genetic Crossover: "), &p->optCrossover);
    setupFloatSlider(slot, x, y += 45, 0.0f, 0.2f, QString("Genetic Mutation: "), &p->optMutation);
    setupIntSlider(slot, x, y += 45, 4, 25, QString("Gene Bits: "), &p->optBits);

	setupCheckBox(slot, x, y += 45, QString("Target"), &p->renderTarget);
	setupCheckBox(slot, x + 90, y, QString("Current"), &p->renderCurrent);

	setupCheckBox(slot, x, y += 45, QString("Objects"), &p->renderObjects);
	setupCheckBox(slot, x+90, y, QString("Point Cloud"), &p->renderPointCloud);

    setupCheckBox(slot, x, y+=45, QString("Generations"), &p->renderGenerations);
    setupCheckBox(slot, x + 90, y, QString("Partial"), &p->isPartial);

    ComboBox *cb = new ComboBox(x, y += 45, 200, 20, "Strategy");
    cb->setVariable(&params::inst()->optStrategy);
    cb->addItem(0, "Best1Exp");
    cb->addItem(1, "Rand1Exp");
    cb->addItem(2, "RandToBest1Exp");
    cb->addItem(3, "Best2Exp");
    cb->addItem(4, "Rand2Exp");
    cb->addItem(5, "Best1Bin");
    cb->addItem(6, "Rand1Bin");
    cb->addItem(7, "RandToBest1Bin");
    cb->addItem(8, "Best2Bin");
    cb->addItem(9, "Rand2Bin");
    cb->setActiveIdx(0);

    ComboBox *cb2 = new ComboBox(x, y += 45, 200, 20, "Optimizer");
    cb2->setVariable(&params::inst()->optTechnique);
    cb2->addItem(0, "CMAES");
    cb2->addItem(1, "DE Solver");
    cb2->addItem(2, "Simulated Annealing");
    cb2->addItem(3, "Swarm Optimization");
    cb2->addItem(4, "Genetic Algorithm");
    cb2->setActiveIdx(0);

    m_guiAll[slot].push_back(cb);
    m_guiAll[slot].push_back(cb2);
}

void GUI::initTab2(Slots slot)
{
    m_tabNames[slot] = "Rendering";

    GlobalObjectParams *p = params::inst();

    int y = 20;
    int x = 10;

    setupFloatSlider(slot, x, y += 45, 0.0f, 8.0f, QString("Shadow A: "), &p->blur.x);
    setupFloatSlider(slot, x, y += 45, 0.0f, 8.0f, QString("Shadow B: "), &p->blur.y);
    setupFloatSlider(slot, x, y += 45, 0.0f, 1.0f, QString("Shadow I: "), &p->shadowIntensity);
    setupFloatSlider(slot, x, y += 45, -5.0f, 5.0f, QString("Shadow I: "), &p->renderHeightOffset);
    setupFloatSlider(slot, x, y += 45, 0.0f, 1.0f, QString("Box Transp: "), &p->boxTransparency);
}

void GUI::initTab3(Slots slot)
{
    m_tabNames[slot] = "Tab3";
}

void GUI::initTab4(Slots slot)
{
    m_tabNames[slot] = "Tab4";
}

void GUI::initTab5(Slots slot)
{
    m_tabNames[slot] = "Tab5";
}

void GUI::setupCheckBox(int slot, int x, int y, const QString &name, bool *var)
{
    CheckBox *cb = new CheckBox(x, y, 15, 15, name);
    cb->setVariable(var);
    cb->setState(*var);
    m_guiAll[slot].push_back(cb);	
}

void GUI::setupFloatSlider(int slot, int x, int y, float rangeMin, float rangeMax, const QString &name, float *var)
{
    Slider<float> *s = new Slider<float>(x, y, 200, 15, name);
    s->setRange(rangeMin, rangeMax);
    s->setVariable(var);
    s->setValue(*var);
    m_guiAll[slot].push_back(s);
}

void GUI::setupIntSlider(int slot, int x, int y, int rangeMin, int rangeMax, const QString &name, int *var)
{
    Slider<int> *s = new Slider<int>(x, y, 200, 15, name);
    s->setRange(rangeMin, rangeMax);
    s->setVariable(var);
    s->setValue(*var);
    m_guiAll[slot].push_back(s);
}

void GUI::render()
{
	m_width = params::inst()->windowSize.x;
	m_height = params::inst()->windowSize.y;

    double frameTime = m_hpTimer.time();

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);

	glDisable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);  
	glDisable(GL_MULTISAMPLE);
    glDisable(GL_TEXTURE_2D);

    //QString helpStr = QString("(F5) Toggle GUI ... ");

    switch (m_mode)
    {
    case 0:         
        //renderString(helpStr.toStdString().c_str(), 10, 15, vec4(m_fontColor.x, m_fontColor.y, m_fontColor.z, 1.0f), GLUT_BITMAP_HELVETICA_12);
        break;
    //case 1:
    //    renderBars()
    //    break;
    default:
        renderTabs();
        break;
    }

    QString fpsStr = QString::number(m_fps);
    renderString(fpsStr.toStdString().c_str(), m_width - 33, 15, vec4(m_fontColor.x, m_fontColor.y, m_fontColor.z, 1.0f), GLUT_BITMAP_HELVETICA_10);

    glPopClientAttrib();
    glPopAttrib();

	m_fpsCount++;

    m_oldTime = frameTime;
}

void GUI::renderTabs()
{
    vec2 size(220, m_height);

    mat4 model = mat4::identitiy();
    mat4 view = mat4::translate(0, 0, -1);
    mat4 projection = mat4::orthographic(0, m_width, m_height, 0, -1, 1);

    Transform trans;
    trans.view = view;
    trans.projection = projection;

    Shader *shader = shaders::inst()->gui;

    shader->bind();

        shader->setMatrix("matView", view, GL_TRUE);
        shader->setMatrix("matProjection", projection, GL_TRUE);

        model = mat4::scale(size.x, size.y, 1.0f);
        shader->setMatrix("matModel", model, GL_TRUE);

        m_vboQuad->render();

    shader->release();

    for (int j = 0; j < m_guiAll.size(); ++j)
    {
        if (m_curTab == j)
        {
            renderString(m_tabNames[j].toStdString().c_str(), 10, 30, vec4(m_fontColor.x, m_fontColor.y, m_fontColor.z, 1.0f), GLUT_BITMAP_HELVETICA_18);

            for (uint i = 0; i < m_guiAll[j].size(); ++i)
            {
                m_guiAll[j][i]->render();
            }
        }
    }

    QString strCamPos = ("LOD CamPos: " + QString::number(params::inst()->camPos.x, 'f', 2) + " " + QString::number(params::inst()->camPos.y, 'f', 2) + " " + QString::number(params::inst()->camPos.z, 'f', 2));
    renderString(strCamPos.toStdString().c_str(), 230, 20, vec4(1.0f, 1.0f, 1.0f, 1.0f), GLUT_BITMAP_HELVETICA_10);

#ifdef WIN32
    vec3 &lightPos = params::inst()->lights[params::inst()->activeLight]->position();
#else
    vec3 lightPos = params::inst()->lights[params::inst()->activeLight]->position();
#endif
    QString strLightPos = ("LightPos: " + QString::number(lightPos.x, 'f', 2) + " " + QString::number(lightPos.y, 'f', 2) + " " + QString::number(lightPos.z, 'f', 2));
    renderString(strLightPos.toStdString().c_str(), 400, 20, vec4(1.0f, 1.0f, 1.0f, 1.0f), GLUT_BITMAP_HELVETICA_10);

    QString strFrameSet = ("Current Frameset: " + m_cameraManager->currentFramesetName());
    renderString(strFrameSet.toStdString().c_str(), 545, 20, vec4(1.0f, 1.0f, 1.0f, 1.0f), GLUT_BITMAP_HELVETICA_10);

    //QString strLightDir = ("LightDir: " + QString::number(params::inst()->lightDir.x, 'f', 2) + " " + QString::number(params::inst()->lightDir.y, 'f', 2) + " " + QString::number(params::inst()->lightDir.z, 'f', 2));
    //renderString(strLightDir.toStdString().c_str(), 10, 60, vec4(1.0f, 1.0f, 1.0f, 1.0f), GLUT_BITMAP_HELVETICA_10);

    QString strVertices = ("Nr Vertices: " + QString::number(params::inst()->nrVertices));
    renderString(strVertices.toStdString().c_str(), 230, 45, vec4(1.0f, 1.0f, 1.0f, 1.0f), GLUT_BITMAP_HELVETICA_10);

    QString strActiveVertices = ("Nr Active Vertices: " + QString::number(params::inst()->nrActiveVertices));
    renderString(strActiveVertices.toStdString().c_str(), 340, 45, vec4(1.0f, 1.0f, 1.0f, 1.0f), GLUT_BITMAP_HELVETICA_10);
}

void GUI::renderBars()
{
    mat4 model = mat4::identitiy();
    mat4 view = mat4::translate(0, 0, -1);
    mat4 projection = mat4::orthographic(0, m_width, m_height, 0, -1, 1);

    Transform trans;
    trans.view = view;
    trans.projection = projection;

    float heightY = 80;
    uint stepX = 170, stepY = 16;

    uint nrX = max<int>(2, m_width / stepX);
    uint nrY = max<int>(2, (heightY-stepY) / stepY);

    uint startX = 10;
    uint startY = m_height - heightY + stepY + 3;
    uint posX = startX, posY = startY;


    Shader *shader = shaders::inst()->gui;
    shader->bind();

        shader->setMatrix("matView", view, GL_TRUE);
        shader->setMatrix("matProjection", projection, GL_TRUE);

        model = mat4::scale(m_width, heightY, 1.0f);
        shader->setMatrix("matModel", model, GL_TRUE);

        m_vboQuad->render();
        
        model = mat4::translate(0.0f, m_height - heightY, 0.0f) * mat4::scale(m_width, heightY, 1.0f);
        shader->setMatrix("matModel", model, GL_TRUE);

        m_vboQuad->render();

    shader->release();

    //top
    QString strCamPos = ("LOD CamPos: " + QString::number(params::inst()->camPos.x, 'f', 2) + " " + QString::number(params::inst()->camPos.y, 'f', 2) + " " + QString::number(params::inst()->camPos.z, 'f', 2));
    renderString(strCamPos.toStdString().c_str(), 10, 20, vec4(1.0f, 1.0f, 1.0f, 1.0f), GLUT_BITMAP_HELVETICA_10);

    //QString strLightPos = ("LightPos: " + QString::number(params::inst()->lightPos.x, 'f', 2) + " " + QString::number(params::inst()->lightPos.y, 'f', 2) + " " + QString::number(params::inst()->lightPos.z, 'f', 2));
    //renderString(strLightPos.toStdString().c_str(), 10, 40, vec4(1.0f, 1.0f, 1.0f, 1.0f), GLUT_BITMAP_HELVETICA_10);

    //QString strLightDir = ("LightDir: " + QString::number(params::inst()->lightDir.x, 'f', 2) + " " + QString::number(params::inst()->lightDir.y, 'f', 2) + " " + QString::number(params::inst()->lightDir.z, 'f', 2));
    //renderString(strLightDir.toStdString().c_str(), 10, 60, vec4(1.0f, 1.0f, 1.0f, 1.0f), GLUT_BITMAP_HELVETICA_10);

    QString strFrameSet = ("Current Frameset: " + m_cameraManager->currentFramesetName());
    renderString(strFrameSet.toStdString().c_str(), 200, 20, vec4(1.0f, 1.0f, 1.0f, 1.0f), GLUT_BITMAP_HELVETICA_10);


    //bottom
    for (uint i = 1; i<m_shortcuts.size() + 1; ++i)
    {
        QString text = m_shortcuts[i - 1];

        renderString(text.toStdString().c_str(), posX, posY, vec4(1.0f, 1.0f, 1.0f, 1.0f), GLUT_BITMAP_HELVETICA_10);

        posY += stepY;

        if ((i%nrY) == 0)
        {
            posY = startY;
            posX += stepX;
        }
    }

    QString fpsStr = QString("FPS: ") + QString::number(m_fps);
    renderString(fpsStr.toStdString().c_str(), m_width - 65, 15, vec4(1.0f, 1.0f, 1.0f, 1.0f), GLUT_BITMAP_HELVETICA_12);
}

void GUI::onTimer()
{
	m_fps = m_fpsCount * 5;
	m_fpsCount = 0;
}

void GUI::toggleMode()
{
    m_mode++;

    if (m_mode >= 1)
    {
        m_curTab++;

        if (m_curTab < m_guiAll.size())
        {
            while (m_guiAll[m_curTab].size() == 0) 
            {
                m_mode++;
                m_curTab++;

                if (m_curTab >= m_guiAll.size())
                {
                    break;
                }
            }                
        }
        else 
        {
            m_curTab = 0;
            m_mode++;
        }
    }

     if (m_mode > m_maxMode)
    {
        m_mode = 0;
        m_curTab = -1;
    }
}

int GUI::currentMode()
{
    return m_mode;
}

void GUI::setFontColor(const vec4 &color)
{
	m_fontColor = color;
}

void GUI::resize(int width, int height)
{
	m_width = width;
	m_height = height;
}

bool GUI::onMouseClick(uint mx, uint my)
{
    bool clicked = false;
 
    for (uint j = 0; j < m_guiAll.size(); ++j)
    {
        if (m_curTab == j)
        {
            for (uint i = 0; i < m_guiAll[j].size(); ++i)
            {
                clicked |= m_guiAll[j][i]->onMouseClick(mx, my);
            }
        }   
    }

	return clicked;
}

void GUI::onMouseMove(uint mx, uint my)
{    
    for (uint j = 0; j < m_guiAll.size(); ++j)
    {
        if (m_curTab == j)
        {
            for (uint i = 0; i < m_guiAll[j].size(); ++i)
            {
                m_guiAll[j][i]->onMouseMove(mx, my);
            }
        }
    }
    
}

void GUI::onMouseRelease()
{
    for (uint j = 0; j < m_guiAll.size(); ++j)
    {
        if (m_curTab == j)
        {
            for (uint i = 0; i < m_guiAll[j].size(); ++i)
            {
                m_guiAll[j][i]->onMouseRelease();
            }
        }
    }
}
