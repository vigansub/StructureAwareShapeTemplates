#include "Headers.h"
#include "DistanceField.h"
#include "VertexBufferObjectAttribs.h"
//#include "VertexBufferObject.h"

#include "Shader.h"
#include "Texture.h"

#include <flann\flann.hpp>
#include <ppl.h>
#include <windows.h> 

DistanceField::DistanceField(int w, int h, int d)  
: m_w(w), m_h(h), m_d(d), m_vbo(nullptr), m_vboLines(nullptr), m_vboPlanes1(nullptr), m_tex1(nullptr), m_vboPlanes2(nullptr), m_tex2(nullptr), m_gridOffset(0)
{
    m_data = new float[m_w*m_h*m_d];
}

DistanceField::~DistanceField()
{
    delete [] m_data;
}

void DistanceField::computeDistance(const vector<vec3> &points)
{
    if(points.size() == 0)
        return;

    int nn = 3;

    float *sourceData = new float[points.size()*4];

    int idxData = 0 ;
    for(int i=0; i<points.size(); ++i)
    {
        vec3 p = points[i];

        sourceData[idxData]   = p.x * m_w;
        sourceData[idxData+1] = p.y * m_h;
        sourceData[idxData+2] = p.z * m_d;

        idxData+=4;
    }

    float *queryData = new float[nr() * 4];

    int w = m_w, h=m_h, d=m_d;    
    int idxQuery = 0;
    for(int z=0; z<d; ++z)
    {
        for(int y=0; y<h; ++y)
        {
            for(int x=0; x<w; ++x)
            {
                queryData[idxQuery]   = x;
                queryData[idxQuery+1] = y;
                queryData[idxQuery+2] = z;

                idxQuery += 4;
            }
        }

    }

    flann::Matrix<float> dataset(sourceData, points.size(), 3, sizeof(float) * 4);
    flann::Matrix<float> query(queryData, nr(), 3, sizeof(float) * 4);

    flann::Matrix<int> indices(new int[query.rows*nn], query.rows, nn);
    flann::Matrix<float> dists(new float[query.rows*nn], query.rows, nn);

    //flann::Index<flann::L2<float> > index(dataset, flann::KDTreeIndexParams(1));
    flann::KDTreeSingleIndex<flann::L2<float> > index(dataset, flann::KDTreeSingleIndexParams(16));
    index.buildIndex();
    
    index.knnSearch(query, indices, dists, nn, flann::SearchParams(128));

    for (int i = 0; i < nr(); ++i)
    {
        m_data[i] = dists[i][0];        
    }    

    m_maxDist = computeMax();
    
    buildVBOLines();
    buildVBOPlanes1();
    buildVBOPlanes2();
    buildVBO();
    buildTextures();
    

    delete[] dataset.ptr();
    delete[] query.ptr();
    delete[] indices.ptr();
    delete[] dists.ptr();    
}

float DistanceField::dist(int x, int y, int z)
{
    int idx = x + y * m_w + z * m_w * m_h;   

    if(idx < 0 || idx > nr())
        return max();
    else
        return m_data[idx];
}

int DistanceField::nr()
{
    return m_w * m_h * m_d;
}

float DistanceField::computeMax()
{
    float ma = math_minfloat;
    for(int i=0; i<nr(); ++i)
    {
        if(m_data[i] > ma)
           ma =  m_data[i];
    }

    return ma;
}

float DistanceField::max()
{
    return m_maxDist;
}

float DistanceField::width()
{
    return m_w;
}

float DistanceField::height()
{
    return m_h;
}

float DistanceField::depth()
{
    return m_d;
}

void DistanceField::print()
{
    for(int i=0; i<nr(); ++i)
    {
        printf("%.1f ", m_data[i]);
        if((i+1)%m_w == 0)
            printf("\n");
    }
}

vector<float> DistanceField::distances()
{
    vector<float> distances;

    for(int i=0; i<nr(); ++i)
    {
        distances.push_back(m_data[i]);
    }

    return distances;
}

void DistanceField::buildVBO()//
{
    delete m_vbo;

    int w = m_w;
    int h = m_h;
    int d = m_d;

    vector<vec3> vertices;
    vector<vec3> normals;

    int count = 0;
    for(int z=0; z<d; z++)
    {
        for(int y=0; y<h; y++)
        {
            for(int x=0; x<w; x++)
            {
                int idx = x + y * w + z * w * h;
                float d = m_data[idx];

                if(d < 1.0)
                {
                    vec3 mi = vec3(0.0, 0.0, 0.0) + vec3(x, y, z);
                    vec3 ma = vec3(1.0, 1.0, 1.0) + vec3(x, y, z);

                    vertices.push_back(vec3(mi.x, ma.y, mi.z));
                    vertices.push_back(vec3(mi.x, ma.y, ma.z));
                    vertices.push_back(vec3(ma.x, ma.y, ma.z));
                    vertices.push_back(vec3(ma.x, ma.y, mi.z));

                    normals.push_back(vec3(0.0f, 1.0f, 0.0f));
                    normals.push_back(vec3(0.0f, 1.0f, 0.0f));
                    normals.push_back(vec3(0.0f, 1.0f, 0.0f));
                    normals.push_back(vec3(0.0f, 1.0f, 0.0f));


                    vertices.push_back(vec3(mi.x, mi.y, mi.z));
                    vertices.push_back(vec3(ma.x, mi.y, mi.z));
                    vertices.push_back(vec3(ma.x, mi.y, ma.z));
                    vertices.push_back(vec3(mi.x, mi.y, ma.z));

                    normals.push_back(vec3(0.0f, -1.0f, 0.0f));
                    normals.push_back(vec3(0.0f, -1.0f, 0.0f));
                    normals.push_back(vec3(0.0f, -1.0f, 0.0f));
                    normals.push_back(vec3(0.0f, -1.0f, 0.0f));    


                    vertices.push_back(vec3(mi.x, mi.y, mi.z));
                    vertices.push_back(vec3(mi.x, ma.y, mi.z));
                    vertices.push_back(vec3(ma.x, ma.y, mi.z));
                    vertices.push_back(vec3(ma.x, mi.y, mi.z));

                    normals.push_back(vec3(0.0f, 0.0f, -1.0f));
                    normals.push_back(vec3(0.0f, 0.0f, -1.0f));
                    normals.push_back(vec3(0.0f, 0.0f, -1.0f));
                    normals.push_back(vec3(0.0f, 0.0f, -1.0f));


                    vertices.push_back(vec3(mi.x, mi.y, ma.z));
                    vertices.push_back(vec3(ma.x, mi.y, ma.z));
                    vertices.push_back(vec3(ma.x, ma.y, ma.z));
                    vertices.push_back(vec3(mi.x, ma.y, ma.z));

                    normals.push_back(vec3(0.0f, 0.0f, 1.0f));
                    normals.push_back(vec3(0.0f, 0.0f, 1.0f));
                    normals.push_back(vec3(0.0f, 0.0f, 1.0f));
                    normals.push_back(vec3(0.0f, 0.0f, 1.0f));


                    vertices.push_back(vec3(mi.x, mi.y, mi.z));
                    vertices.push_back(vec3(mi.x, mi.y, ma.z));
                    vertices.push_back(vec3(mi.x, ma.y, ma.z));
                    vertices.push_back(vec3(mi.x, ma.y, mi.z));

                    normals.push_back(vec3(-1.0f, 0.0f, 0.0f));
                    normals.push_back(vec3(-1.0f, 0.0f, 0.0f));
                    normals.push_back(vec3(-1.0f, 0.0f, 0.0f));
                    normals.push_back(vec3(-1.0f, 0.0f, 0.0f));


                    vertices.push_back(vec3(ma.x, mi.y, mi.z));
                    vertices.push_back(vec3(ma.x, ma.y, mi.z));
                    vertices.push_back(vec3(ma.x, ma.y, ma.z));
                    vertices.push_back(vec3(ma.x, mi.y, ma.z));

                    normals.push_back(vec3(1.0f, 0.0f, 0.0f));
                    normals.push_back(vec3(1.0f, 0.0f, 0.0f));
                    normals.push_back(vec3(1.0f, 0.0f, 0.0f));
                    normals.push_back(vec3(1.0f, 0.0f, 0.0f));
                    
                    count ++;
                }
            }
        }
    }

    writeObj(vertices, normals);

    int s = m_w * m_h * m_d;
    //qDebug() << s << count << (float)count / (float)s;

   vec4 color = vec4(1, 0, 0, 1);

   uint nrVertices = vertices.size();
   VertexBufferObjectAttribs::DATA *attrData = new VertexBufferObjectAttribs::DATA[nrVertices];

    for(uint i=0; i<nrVertices; ++i)
    {    
        vec3 v = vertices[i];
        vec3 n = normals[i];

        attrData[i].vx = v.x;
        attrData[i].vy = v.y;
        attrData[i].vz = v.z;
        attrData[i].vw = 1.0f;

        attrData[i].nx = n.x;
        attrData[i].ny = n.y;
        attrData[i].nz = n.z;
        attrData[i].nw = 0.0f;

        attrData[i].cx = color.x;
        attrData[i].cy = color.y;
        attrData[i].cz = color.z;
        attrData[i].cw = color.w;

        attrData[i].tx = 0.0f;
        attrData[i].ty = 0.0f;
        attrData[i].tz = 0.0f;
        attrData[i].tw = 0.0f;
    }

    m_vbo = new VertexBufferObjectAttribs();
    m_vbo->setData(attrData, GL_STATIC_DRAW, nrVertices, GL_QUADS);

    m_vbo->addAttrib(VERTEX_POSITION);
    m_vbo->addAttrib(VERTEX_NORMAL);
    m_vbo->addAttrib(VERTEX_COLOR);
    m_vbo->addAttrib(VERTEX_TEXTURE);
    m_vbo->bindAttribs();

    delete[] attrData;    
}

void DistanceField::buildVBOPlanes1()
{
    delete m_vboPlanes1;

    int w = m_w;
    int h = m_h;
    int d = m_d;

    vector<vec3> vertices, texCoord;

    vec3 mi = vec3(0, 0, d/2);
    vec3 ma = vec3(w, h, d/2);

    vertices.push_back(vec3(mi.x, mi.y, mi.z));
    vertices.push_back(vec3(ma.x, mi.y, mi.z));
    vertices.push_back(vec3(ma.x, ma.y, ma.z));
    vertices.push_back(vec3(mi.x, ma.y, ma.z));

    texCoord.push_back(vec3(0, 0, 0));
    texCoord.push_back(vec3(1, 0, 0));
    texCoord.push_back(vec3(1, 1, 0));
    texCoord.push_back(vec3(0, 1, 0));

   vec4 color = vec4(0, 0, 1, 1);

   uint nrVertices = vertices.size();
   VertexBufferObjectAttribs::DATA *attrData = new VertexBufferObjectAttribs::DATA[nrVertices];

    for(uint i=0; i<nrVertices; ++i)
    {    
        vec3 v = vertices[i];
        vec3 n = vec3();
        vec3 t = texCoord[i];

        attrData[i].vx = v.x;
        attrData[i].vy = v.y;
        attrData[i].vz = v.z;
        attrData[i].vw = 1.0f;

        attrData[i].nx = n.x;
        attrData[i].ny = n.y;
        attrData[i].nz = n.z;
        attrData[i].nw = 0.0f;

        attrData[i].cx = color.x;
        attrData[i].cy = color.y;
        attrData[i].cz = color.z;
        attrData[i].cw = color.w;

        attrData[i].tx = t.x;
        attrData[i].ty = t.y;
        attrData[i].tz = 0.0f;
        attrData[i].tw = 0.0f;
    }

    m_vboPlanes1 = new VertexBufferObjectAttribs();
    m_vboPlanes1->setData(attrData, GL_STATIC_DRAW, nrVertices, GL_QUADS);

    m_vboPlanes1->addAttrib(VERTEX_POSITION);
    m_vboPlanes1->addAttrib(VERTEX_NORMAL);
    m_vboPlanes1->addAttrib(VERTEX_COLOR);
    m_vboPlanes1->addAttrib(VERTEX_TEXTURE);
    m_vboPlanes1->bindAttribs();

    delete[] attrData;    
}

void DistanceField::buildVBOPlanes2()
{
    delete m_vboPlanes2;

    int w = m_w;
    int h = m_h;
    int d = m_d;

    vector<vec3> vertices, texCoord;

    vec3 mi = vec3(0, h/2, 0);
    vec3 ma = vec3(w, h/2, d);

    vertices.push_back(vec3(mi.x, mi.y, mi.z));
    vertices.push_back(vec3(ma.x, mi.y, mi.z));
    vertices.push_back(vec3(ma.x, ma.y, ma.z));
    vertices.push_back(vec3(mi.x, ma.y, ma.z));

    texCoord.push_back(vec3(0, 1, 0));
    texCoord.push_back(vec3(1, 1, 0));
    texCoord.push_back(vec3(1, 0, 0));
    texCoord.push_back(vec3(0, 0, 0));

   vec4 color = vec4(0, 0, 1, 1);

   uint nrVertices = vertices.size();
   VertexBufferObjectAttribs::DATA *attrData = new VertexBufferObjectAttribs::DATA[nrVertices];

    for(uint i=0; i<nrVertices; ++i)
    {    
        vec3 v = vertices[i];
        vec3 n = vec3();
        vec3 t = texCoord[i];

        attrData[i].vx = v.x;
        attrData[i].vy = v.y;
        attrData[i].vz = v.z;
        attrData[i].vw = 1.0f;

        attrData[i].nx = n.x;
        attrData[i].ny = n.y;
        attrData[i].nz = n.z;
        attrData[i].nw = 0.0f;

        attrData[i].cx = color.x;
        attrData[i].cy = color.y;
        attrData[i].cz = color.z;
        attrData[i].cw = color.w;

        attrData[i].tx = t.x;
        attrData[i].ty = t.y;
        attrData[i].tz = 0.0f;
        attrData[i].tw = 0.0f;
    }

    m_vboPlanes2 = new VertexBufferObjectAttribs();
    m_vboPlanes2->setData(attrData, GL_STATIC_DRAW, nrVertices, GL_QUADS);

    m_vboPlanes2->addAttrib(VERTEX_POSITION);
    m_vboPlanes2->addAttrib(VERTEX_NORMAL);
    m_vboPlanes2->addAttrib(VERTEX_COLOR);
    m_vboPlanes2->addAttrib(VERTEX_TEXTURE);
    m_vboPlanes2->bindAttribs();

    delete[] attrData;    
}

void DistanceField::buildVBOLines()
{
    delete m_vboLines;

    int w = m_w;
    int h = m_h;
    int d = m_d;

    vector<vec3> vertices;
    vector<vec3> normals;

    int count = 0;
    for(int z=0; z<d; z++)
    {
        for(int y=0; y<h; y++)
        {
            for(int x=0; x<w; x++)
            {
                int idx = x + y * w + z * w * h;
                float d = m_data[idx];

                if(d < 1.0)
                {
                    vec3 mi = vec3(0.0, 0.0, 0.0) + vec3(x, y, z);
                    vec3 ma = vec3(1.0, 1.0, 1.0) + vec3(x, y, z);

                    vertices.push_back(vec3(mi.x, ma.y, mi.z));
                    vertices.push_back(vec3(mi.x, ma.y, ma.z));

                    vertices.push_back(vec3(mi.x, ma.y, ma.z));
                    vertices.push_back(vec3(ma.x, ma.y, ma.z));

                    vertices.push_back(vec3(ma.x, ma.y, ma.z));
                    vertices.push_back(vec3(ma.x, ma.y, mi.z));

                    vertices.push_back(vec3(ma.x, ma.y, mi.z));
                    vertices.push_back(vec3(mi.x, ma.y, mi.z));


                    vertices.push_back(vec3(mi.x, mi.y, mi.z));
                    vertices.push_back(vec3(ma.x, mi.y, mi.z));

                    vertices.push_back(vec3(ma.x, mi.y, mi.z));
                    vertices.push_back(vec3(ma.x, mi.y, ma.z));

                    vertices.push_back(vec3(ma.x, mi.y, ma.z));
                    vertices.push_back(vec3(mi.x, mi.y, ma.z));

                    vertices.push_back(vec3(mi.x, mi.y, ma.z));
                    vertices.push_back(vec3(mi.x, mi.y, mi.z));


                    vertices.push_back(vec3(mi.x, mi.y, mi.z));
                    vertices.push_back(vec3(mi.x, ma.y, mi.z));

                    vertices.push_back(vec3(mi.x, ma.y, mi.z));
                    vertices.push_back(vec3(ma.x, ma.y, mi.z));

                    vertices.push_back(vec3(ma.x, ma.y, mi.z));
                    vertices.push_back(vec3(ma.x, mi.y, mi.z));

                    vertices.push_back(vec3(ma.x, mi.y, mi.z));
                    vertices.push_back(vec3(mi.x, mi.y, mi.z));


                    vertices.push_back(vec3(mi.x, mi.y, ma.z));
                    vertices.push_back(vec3(ma.x, mi.y, ma.z));

                    vertices.push_back(vec3(ma.x, mi.y, ma.z));
                    vertices.push_back(vec3(ma.x, ma.y, ma.z));

                    vertices.push_back(vec3(ma.x, ma.y, ma.z));
                    vertices.push_back(vec3(mi.x, ma.y, ma.z));

                    vertices.push_back(vec3(mi.x, ma.y, ma.z));
                    vertices.push_back(vec3(mi.x, mi.y, ma.z));


                    vertices.push_back(vec3(mi.x, mi.y, mi.z));
                    vertices.push_back(vec3(mi.x, mi.y, ma.z));

                    vertices.push_back(vec3(mi.x, mi.y, ma.z));
                    vertices.push_back(vec3(mi.x, ma.y, ma.z));

                    vertices.push_back(vec3(mi.x, ma.y, ma.z));
                    vertices.push_back(vec3(mi.x, ma.y, mi.z));

                    vertices.push_back(vec3(mi.x, ma.y, mi.z));
                    vertices.push_back(vec3(mi.x, mi.y, mi.z));


                    vertices.push_back(vec3(ma.x, mi.y, mi.z));
                    vertices.push_back(vec3(ma.x, ma.y, mi.z));

                    vertices.push_back(vec3(ma.x, ma.y, mi.z));
                    vertices.push_back(vec3(ma.x, ma.y, ma.z));

                    vertices.push_back(vec3(ma.x, ma.y, ma.z));
                    vertices.push_back(vec3(ma.x, mi.y, ma.z));

                    vertices.push_back(vec3(ma.x, mi.y, ma.z));
                    vertices.push_back(vec3(ma.x, mi.y, mi.z));

                    count ++;
                }
            }
        }
    }

    int s = m_w * m_h * m_d;
    ///qDebug() << s << count << (float)count / (float)s;

   vec4 color = vec4(0, 0, 0, 1);

   uint nrVertices = vertices.size();
   VertexBufferObjectAttribs::DATA *attrData = new VertexBufferObjectAttribs::DATA[nrVertices];

    for(uint i=0; i<nrVertices; ++i)
    {    
        vec3 v = vertices[i];
        vec3 n = vec3();//normals[i];

        attrData[i].vx = v.x;
        attrData[i].vy = v.y;
        attrData[i].vz = v.z;
        attrData[i].vw = 1.0f;

        attrData[i].nx = n.x;
        attrData[i].ny = n.y;
        attrData[i].nz = n.z;
        attrData[i].nw = 0.0f;

        attrData[i].cx = color.x;
        attrData[i].cy = color.y;
        attrData[i].cz = color.z;
        attrData[i].cw = color.w;

        attrData[i].tx = 0.0f;
        attrData[i].ty = 0.0f;
        attrData[i].tz = 0.0f;
        attrData[i].tw = 0.0f;
    }

    m_vboLines = new VertexBufferObjectAttribs();
    m_vboLines->setData(attrData, GL_STATIC_DRAW, nrVertices, GL_LINES);

    m_vboLines->addAttrib(VERTEX_POSITION);
    m_vboLines->addAttrib(VERTEX_NORMAL);
    m_vboLines->addAttrib(VERTEX_COLOR);
    m_vboLines->addAttrib(VERTEX_TEXTURE);
    m_vboLines->bindAttribs();

    delete[] attrData;    
}

void DistanceField::render(const Transform &trans)
{
    mat4 projection = trans.projection;
    mat4 view = trans.view;

    float s = 1.0f / g_max<float>(m_w, g_max<float>(m_h, m_d));

    float t = 1.0f;
    mat4 model = mat4::translate(-t, -t, -t) * mat4::scale(vec3(2, 2, 2) * s);
    //mat4 model = mat4::translate(-t, -t, -t) * s;

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);

    if(params::inst()->renderOccupancy)
    {
	    Shader *shader = shaders::inst()->normalLight;	
	    shader->bind();
        
            shader->set3f("lightPos", params::inst()->lightPos);
     	    shader->setMatrix("matProjection", projection, GL_TRUE);
		    shader->setMatrix("matView", view, GL_TRUE);			
		    shader->setMatrix("matModel", model, GL_TRUE);
		
            if(m_vbo)
                m_vbo->render();

            if(m_vboLines)
                m_vboLines->render();
		
	    shader->release();    
    }
    
    if(params::inst()->renderDistTexture)
    {
	    Shader *shaderP1 = shaders::inst()->planes1;	
	    shaderP1 ->bind();

     	    shaderP1->setMatrix("matProjection", projection, GL_TRUE);
		    shaderP1->setMatrix("matView", view, GL_TRUE);			
		    shaderP1->setMatrix("matModel", model, GL_TRUE);

            if(m_tex1)
            {
                glActiveTexture(GL_TEXTURE0);
                glBindTexture(GL_TEXTURE_2D, m_tex1->id());    
                shaderP1->seti("tex", 0);    
            }      
		
            if(m_vboPlanes1)
                m_vboPlanes1->render();

        shaderP1->release();


	    Shader *shaderP2 = shaders::inst()->planes2;	
	    shaderP2 ->bind();

     	    shaderP2->setMatrix("matProjection", projection, GL_TRUE);
		    shaderP2->setMatrix("matView", view, GL_TRUE);			
		    shaderP2->setMatrix("matModel", model, GL_TRUE);

            if(m_tex2)
            {
                glActiveTexture(GL_TEXTURE0);
                glBindTexture(GL_TEXTURE_2D, m_tex2->id());    
                shaderP2->seti("tex", 0);    
            }

            //float uintY = 2.0f / params::inst()->DFDimension;
            //model = mat4::translate(0, uintY*m_gridOffset, 0) * model;
            shaderP2 ->setMatrix("matModel", model, GL_TRUE);

            if(m_vboPlanes2)
                m_vboPlanes2->render();
		
	    shaderP2->release();    
    }


	glPopClientAttrib();
    glPopAttrib();
}

void DistanceField::save(const QString &path)
{
    QFile file(path);
    if (!file.open(QIODevice::WriteOnly))
        return;

    QDataStream out(&file); 

    for(int i=0; i<nr(); ++i)
    {
        out << m_data[i];
    }

    file.close();
}

void DistanceField::load(const QString &path)
{
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly))
        return;

    memset(m_data, 0, nr() * sizeof(float));

    QDataStream in(&file); 

    for(int i=0; i<nr(); ++i)
    {
        in >> m_data[i];
    }

    m_maxDist = computeMax();
    buildVBO();

    file.close();
}

void DistanceField::writeObj(vector<vec3> &positions, vector<vec3> &normals)
{
	vector<Face> faces;
	for (int i = 0; i<positions.size() - 3; i += 4)
	{
		Face f;
		f.a = i + 0 + 1;
		f.b = i + 1 + 1;
		f.c = i + 2 + 1;
        f.d = i + 3 + 1;

		faces.push_back(f);
	}

	QFile file("df_" + QString::number(params::inst()->DFDimension) + ".obj");
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
		return;

	QTextStream out(&file);

	for (int i = 0; i<positions.size(); ++i)
	{
		vec3 p = positions[i];
		out << "v " << p.x << " " << p.y << " " << p.z << endl;
	}

	for (int i = 0; i<normals.size(); ++i)
	{
		vec3 n = normals[i];
		out << "vn " << n.x << " " << n.y << " " << n.z << endl;
	}

	for (int i = 0; i<faces.size(); ++i)
	{
		Face f = faces[i];

		out << "f " << f.a << "/" << f.a << "/" << f.a << " " << f.b << "/" << f.b << "/" << f.b
			<< " " << f.c << "/" << f.c << "/" << f.c << " " << f.d << "/" << f.d << "/" << f.d << endl;
	}

	file.close();
}

void DistanceField::buildTextures()
{
    delete m_tex1;

    int w = m_w;
    int h = m_h;
    int d = m_d;

    float *data1 = new float[w * h * 4];

    //for(int z=0; z<d; z++)

    int ref = 0;
    int z = d/2;
    {
        for(int y=0; y<h; y++)
        {
            for(int x=0; x<w; x++)
            {
                int idx = x + y * w + z * w * h;
                
                float d  = sqrt(m_data[idx] / max());

                data1[ref]   = d;
                data1[ref+1] = d;
                data1[ref+2] = d;
                data1[ref+3] = 1;

                ref+=4;
            }
        }
    }

    m_tex1 = new Texture(w, h, GL_RGBA, GL_RGBA, GL_FLOAT, data1);
    m_tex1->setFilter(GL_LINEAR, GL_LINEAR);
    m_tex1->setWrapMode(GL_MIRRORED_REPEAT);




    float *data2 = new float[w * h * 4];

    //for(int z=0; z<d; z++)

    ref = 0;
    for(int z = 0; z<d; z++)
    {
        //for(int y=0; y<h; y++)
        int y = h/2 + m_gridOffset;
        {
            for(int x=0; x<w; x++)
            {
                int idx = x + y * w + z * w * h;
                
                float d  = sqrt(m_data[idx] / max());

                data2[ref]   = d;
                data2[ref+1] = d;
                data2[ref+2] = d;
                data2[ref+3] = 1;

                ref+=4;
            }
        }
    }

    m_tex2 = new Texture(w, h, GL_RGBA, GL_RGBA, GL_FLOAT, data2);
    m_tex2->setFilter(GL_LINEAR, GL_LINEAR);
    m_tex2->setWrapMode(GL_MIRRORED_REPEAT);



    delete data1;
    delete data2;
}