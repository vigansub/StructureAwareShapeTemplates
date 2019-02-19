#include "Scanner.h"
#include "Shader.h"
#include "VertexBufferObject.h"
#include "FrameBufferObject.h"
#include "Scene.h"


Scanner::Scanner(Scene *scene) 
: m_scene(scene), 
  m_bufferWidth(0), 
  m_bufferHeight(0), 
  m_nrBufferVertices(100000), 
  m_moveAutomatic(true), 
  m_fcp(0.01),
  m_ncp(100),
  m_aspect(1.0),
  m_fov(45.0f), 
  m_activePos(0, 0), 
  m_rotHeight(0.0f), 
  m_zoom(-20.0), 
  m_fbo(nullptr), 
  m_scanning(false)
{
    buildVBO();
}

Scanner::~Scanner()
{
    delete m_vbo;
    delete m_fbo;
}

void Scanner::scan()
{
    Transform t;
    renderScan(t, false);

    storeView();
}

void Scanner::autoScan()
{
    if(m_scanning)
    {
        scan();
        int result = setToNextPos();

        if(result == -1)
        {
            m_scanning = false;
        }
    }
}

// TODO or TOREMOVE
vector<vec3> Scanner::scan_data()
{
    vector<vec3> results;
    m_scanning = true;
    scan();
    int result = setToNextPos();

    if(result == -1)
    {
        m_scanning = false;
    }
    return results;
}

void Scanner::renderScan(const Transform &trans, bool useView)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);

    Transform t;

    t.projection = trans.projection;
    t.view = trans.view;

    if(!useView)
    {
        mat4 projection = mat4::perspective(m_fov, m_aspect, m_ncp, m_fcp);
	    mat4 view = mat4::identitiy();

	    view = mat4::translate(vec3(0, m_rotHeight, 0));
	    view *= mat4::translate(vec3(0, 0, -5));
	    view *= mat4::rotate(m_activePos.x, 1, 0, 0);
	    view *= mat4::rotate(m_activePos.y, 0, 1, 0); 

        t.projection = projection;
        t.view = view;
    }

    m_fbo->bind();

        glViewport(0, 0, m_bufferWidth, m_bufferHeight);
        glClearColor(0, 0, 0, 1);    
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);     

        glEnable(GL_MULTISAMPLE);        

        if (params::inst()->polygonMode == 1)
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        else
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		
        m_scene->renderObjectsScan(t);

    m_fbo->release();          

	glPopClientAttrib();
    glPopAttrib();
}

void Scanner::storeView()
{
    qDebug() << "Scanner::storeView()" << endl;
    int w = m_bufferWidth;
    int h = m_bufferHeight;
    int s = w * h * 4;

    float *dataPos = new float[w * h * 4];
    float *dataNor = new float[w * h * 4];

    fill_n(dataPos, s, 0.0f);
    fill_n(dataNor, s, 0.0f);

    m_fbo->bind();
       	    
        glReadBuffer(GL_COLOR_ATTACHMENT0);
	    glReadPixels(0, 0, w, h, GL_RGBA, GL_FLOAT, dataPos);	 

        glReadBuffer(GL_COLOR_ATTACHMENT1);
	    glReadPixels(0, 0, w, h, GL_RGBA, GL_FLOAT, dataNor);	 

    m_fbo->release();

    int sum = 0;
    for(int i=0; i<s; i+=4)
    {        
        vec3 p = vec3(dataPos[i], dataPos[i+1], dataPos[i+2]);        

        if(length(p) > 0.0f) 
        {                
			int prec = 6;// params::inst()->scanPrecision;
            
            Vertex v;
            vec3 n = vec3(dataNor[i], dataNor[i+1], dataNor[i+2]);

            float rnd_range = 0.001f;
            v.position = p + vec3(rand(-rnd_range, rnd_range), rand(-rnd_range, rnd_range), rand(-rnd_range, rnd_range));
            v.normal = n;

            QString x = QString("%1").arg(p.x, 0, 'f', prec);
            QString y = QString("%1").arg(p.y, 0, 'f', prec);
            QString z = QString("%1").arg(p.z, 0, 'f', prec);

            QString str = x + y + z;
            //m_vertices.insert(str, v);
            m_vertices[str] = v;
			m_scene->m_verticesScan.insert(str, v);
            sum ++;
        }
    }

    updateBuffer();      

    delete [] dataPos;
    delete [] dataNor;
}

void Scanner::updateBuffer()
{
    m_vbo->bind();

        int size = 0;
	    glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);
	    size /= sizeof(VertexBufferObject::DATA);

	    VertexBufferObject::DATA *attrData = (VertexBufferObject::DATA *)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);   

        memset(attrData, 0, size * sizeof(VertexBufferObject::DATA));
        
        QList<Vertex> vertices = m_scene->m_verticesScan.values();//m_vertices.values();

        QList<Vertex>::iterator iter = vertices.begin();

        for(int i=0, j=0; i<size && j<vertices.size(); ++i, ++j)
        {
                Vertex *v = &(*iter);

                attrData[i].vx = v->position.x;
                attrData[i].vy = v->position.y;
                attrData[i].vz = v->position.z;
                attrData[i].vw = 1.0f;

                attrData[i].nx = v->normal.x;
                attrData[i].ny = v->normal.y;
                attrData[i].nz = v->normal.z;
                attrData[i].nw = 0.0f;

                attrData[i].cx = 0.0f;
                attrData[i].cy = 0.0f;
                attrData[i].cz = 0.0f;
                attrData[i].cw = 1.0f;

                attrData[i].tx = 0.0f;
                attrData[i].ty = 0.0f;
                attrData[i].tz = 0.0f;
		        attrData[i].tw = 0.0f;        

                iter++;
        }

        glUnmapBuffer(GL_ARRAY_BUFFER);

    m_vbo->release();
}

void Scanner::resize(int w, int h)
{
    m_bufferWidth = w;
    m_bufferHeight = h;

    m_fbo = new FrameBufferObject(m_bufferWidth, m_bufferHeight, 2, 0, GL_TRUE, GL_NEAREST);
}

void Scanner::buildVBO()
{
   uint nrVertices = m_nrBufferVertices;
   VertexBufferObject::DATA *attrData = new VertexBufferObject::DATA[nrVertices];

    for(uint i=0; i<nrVertices; ++i)
    {    
        attrData[i].vx = 0.0f;
        attrData[i].vy = 0.0f;
        attrData[i].vz = 0.0f;
        attrData[i].vw = 1.0f;

        attrData[i].nx = 0.0f;
        attrData[i].ny = 0.0f;
        attrData[i].nz = 0.0f;
        attrData[i].nw = 0.0f;

        attrData[i].cx = 1.0f;
        attrData[i].cy = 1.0f;
        attrData[i].cz = 0.0f;
        attrData[i].cw = 1.0f;

        attrData[i].tx = 0.0f;
        attrData[i].ty = 0.0f;
        attrData[i].tz = 0.0f;
		attrData[i].tw = 0.0f;
    }

    m_vbo = new VertexBufferObject();
    m_vbo->setData(attrData, GL_STATIC_DRAW, nrVertices, GL_POINTS);
	m_vbo->bindDefaultAttribs();

    delete[] attrData;    
}

int Scanner::setToNextPos()
{
    int result = 0;

    m_activePos.x += 60;
    if(m_activePos.x >= 360)
    {
        m_activePos.x = 0;
        m_activePos.y += 60;
    }

    if(m_activePos.y >= 360)
    {
        m_activePos.y = 0;
        result = -1;
    }

    return result;
}

void Scanner::renderPointCloud(const Transform &trans)
{
    mat4 projection = trans.projection;
    mat4 view = trans.view;
    mat4 model = mat4::identitiy();

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);

    glPointSize(2.0f);

#ifndef WIN32
    Shader *shader = shaders::inst()->default_;
#else
    Shader *shader = shaders::inst()->default;
#endif

    shader->bind();
        shader->setMatrix("matProjection", projection, GL_TRUE);
        shader->setMatrix("matView", view, GL_TRUE);
        shader->setMatrix("matModel", model, GL_TRUE); 
        
        if(m_vbo)// && params::inst()->renderScan)
            m_vbo->render();

    shader->release();

	glPopClientAttrib();
    glPopAttrib();
}

void Scanner::renderPointCloudDepth(const Transform &trans)
{
    mat4 projection = trans.projection;
    mat4 view = trans.view;
    mat4 model = mat4::identitiy();

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);

    glPointSize(1.0f);

#ifndef WIN32
    Shader *shader = shaders::inst()->defaultDepth;
#else
    Shader *shader = shaders::inst()->defaultDepth;
#endif

    shader->bind();
    shader->setMatrix("matProjection", projection, GL_TRUE);
    shader->setMatrix("matView", view, GL_TRUE);
    shader->setMatrix("matModel", model, GL_TRUE);

    if (m_vbo)// && params::inst()->renderScan)
        m_vbo->render();

    shader->release();

    glPopClientAttrib();
    glPopAttrib();
}

void Scanner::clearBuffer()
{
    m_vertices.clear();

    m_vbo->bind();

    int size = 0;
	glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);
	size /= sizeof(VertexBufferObject::DATA);

	VertexBufferObject::DATA *attrData = (VertexBufferObject::DATA *)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);   

    for(int i=0; i<size ; ++i)
    {
            attrData[i].vx = 0.0f;
            attrData[i].vy = 0.0f;
            attrData[i].vz = 0.0f;
            attrData[i].vw = 0.0f;

            attrData[i].nx = 0.0f;
            attrData[i].ny = 0.0f;
            attrData[i].nz = 0.0f;
            attrData[i].nw = 0.0f;

            attrData[i].cx = 0.0f;
            attrData[i].cy = 0.0f;
            attrData[i].cz = 0.0f;
            attrData[i].cw = 0.0f;

            attrData[i].tx = 0.0f;
            attrData[i].ty = 0.0f;
            attrData[i].tz = 0.0f;
		    attrData[i].tw = 0.0f;        
    }

    glUnmapBuffer(GL_ARRAY_BUFFER);
    m_vbo->release();
}
