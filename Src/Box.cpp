#include "Box.h"
#include "Shader.h"
#include "VertexBufferObject.h"
#include "Light.h"

Box::Box() 
 : m_size(0.5f), 
 m_vbo(nullptr)
{
    init();
}

Box::~Box()
{
    delete m_vbo;
}

void Box::init()
{   
    float s = m_size;

    m_vertices.push_back(vec3(-s, -s,  s)); // a
    m_vertices.push_back(vec3( s, -s,  s)); // b
    m_vertices.push_back(vec3( s, -s, -s)); // c
    m_vertices.push_back(vec3(-s, -s, -s)); // d
    m_vertices.push_back(vec3(-s,  s,  s)); // e
    m_vertices.push_back(vec3( s,  s,  s)); // f
    m_vertices.push_back(vec3( s,  s, -s)); // g
    m_vertices.push_back(vec3(-s,  s, -s)); // h

    Face front, back, right, left, top, bottom;

    front.a = 0;
    front.b = 1;
    front.c = 5;
    front.d = 4;
    front.type = FRONT;
    front.color = vec4(1, 0, 0, 1); 
    front.normal = vec3(0, 0, 1);

    back.a = 2;
    back.b = 3;
    back.c = 7;
    back.d = 6;
    back.type = BACK;
    back.color = vec4(0.5, 0, 0, 1);
    back.normal = vec3(0, 0, -1);

    right.a = 1;
    right.b = 2;
    right.c = 6;
    right.d = 5;
    right.type = RIGHT;
    right.color = vec4(0, 1.0, 0, 1);
    right.normal = vec3(1, 0, 0);

    left.a = 3;
    left.b = 0;
    left.c = 4;
    left.d = 7;
    left.type = LEFT;
    left.color = vec4(0, 0.5, 0, 1);
    left.normal = vec3(-1, 0, 0);

    top.a = 4;
    top.b = 5;
    top.c = 6;
    top.d = 7;
    top.type = TOP;
    top.color = vec4(0, 0, 1, 1);
    top.normal = vec3(0, 1, 0);

    bottom.a = 3;
    bottom.b = 2;
    bottom.c = 1;
    bottom.d = 0;
    bottom.type = BOTTOM;
    bottom.color = vec4(0, 0, 0.5, 1);
    bottom.normal = vec3(0, -1, 0);

    m_faces.resize(6);
    m_faces[FRONT] = front;
    m_faces[BACK] = back;
    m_faces[LEFT] = right;
    m_faces[RIGHT] = left;
    m_faces[TOP] = top;
    m_faces[BOTTOM] = bottom;

    buildVBO();
    buildVBOLines();
}

void Box::buildVBO()
{
    vector<vec3> vertices;
    vector<vec4> colors;
    vector<vec3> normals;

    for (int i = 0; i < m_faces.size(); ++i)
    {
        Face &f = m_faces[i];

        vec4 color = f.color;
        vec3 normal = f.normal;

        vec3 a = m_vertices[f.a];
        vec3 b = m_vertices[f.b];
        vec3 c = m_vertices[f.c];
        vec3 d = m_vertices[f.d];

        vertices.push_back(a);
        vertices.push_back(b);
        vertices.push_back(c);
        vertices.push_back(c);
        vertices.push_back(d);
        vertices.push_back(a);

        normals.push_back(normal);
        normals.push_back(normal);
        normals.push_back(normal);
        normals.push_back(normal);
        normals.push_back(normal);
        normals.push_back(normal);

        colors.push_back(color);
        colors.push_back(color);
        colors.push_back(color);
        colors.push_back(color);
        colors.push_back(color);
        colors.push_back(color);
    }

    uint nrVertices = vertices.size();
    VertexBufferObject::DATA *attrData = new VertexBufferObject::DATA[nrVertices];

    for (uint i = 0; i<nrVertices; ++i)
    {
        vec3 v = vertices[i];
        vec4 c = colors[i];
        vec3 n = normals[i];

        attrData[i].vx = v.x;
        attrData[i].vy = v.y;
        attrData[i].vz = v.z;
        attrData[i].vw = 1.0f;

        attrData[i].nx = n.x;
        attrData[i].ny = n.y;
        attrData[i].nz = n.z;
        attrData[i].nw = 1.0;

        attrData[i].cx = c.x;
        attrData[i].cy = c.y;
        attrData[i].cz = c.z;
        attrData[i].cw = c.w;

        attrData[i].tx = 0.0f;
        attrData[i].ty = 0.0f;
        attrData[i].tz = 0.0f;
        attrData[i].tw = 0.0f;
    }

    m_vbo = new VertexBufferObject();
    m_vbo->setData(attrData, GL_STATIC_DRAW, nrVertices, GL_TRIANGLES);

    m_vbo->addAttrib(VERTEX_POSITION);
    m_vbo->addAttrib(VERTEX_NORMAL);
    m_vbo->addAttrib(VERTEX_COLOR);
    m_vbo->addAttrib(VERTEX_TEXTURE);
    m_vbo->bindAttribs();

    delete[] attrData;
}

void Box::buildVBOLines()
{
    vector<vec3> vertices;
    vector<vec4> colors;
    vector<vec3> normals;

    for (int i = 0; i < m_faces.size(); ++i)
    {
        Face &f = m_faces[i];

        vec4 color = f.color;
        vec3 normal = f.normal;

        vec3 a = m_vertices[f.a];
        vec3 b = m_vertices[f.b];
        vec3 c = m_vertices[f.c];
        vec3 d = m_vertices[f.d];

        vertices.push_back(a);
        vertices.push_back(b);

        vertices.push_back(b);
        vertices.push_back(c);

        vertices.push_back(c);
        vertices.push_back(d);

        vertices.push_back(d);
        vertices.push_back(a);
    }

    uint nrVertices = vertices.size();
    VertexBufferObject::DATA *attrData = new VertexBufferObject::DATA[nrVertices];

    for (uint i = 0; i<nrVertices; ++i)
    {
        vec3 v = vertices[i];

        attrData[i].vx = v.x;
        attrData[i].vy = v.y;
        attrData[i].vz = v.z;
        attrData[i].vw = 1.0f;

        attrData[i].nx = 0.0f;
        attrData[i].ny = 0.0f;
        attrData[i].nz = 0.0f;
        attrData[i].nw = 0.0f;

        attrData[i].cx = 0.0f;
        attrData[i].cy = 0.0f;
        attrData[i].cz = 0.0f;
        attrData[i].cw = 1.0f;

        attrData[i].tx = 0.0f;
        attrData[i].ty = 0.0f;
        attrData[i].tz = 0.0f;
        attrData[i].tw = 0.0f;
    }

    m_vboLines = new VertexBufferObject();
    m_vboLines->setData(attrData, GL_STATIC_DRAW, nrVertices, GL_LINES);

    m_vboLines->addAttrib(VERTEX_POSITION);
    m_vboLines->addAttrib(VERTEX_NORMAL);
    m_vboLines->addAttrib(VERTEX_COLOR);
    m_vboLines->addAttrib(VERTEX_TEXTURE);
    m_vboLines->bindAttribs();

    delete[] attrData;
}

void Box::render(const Transform &trans, const mat4 &model, int colorMode, int shaderSelector, bool isConnector)
{
    mat4 projection = trans.projection;
    mat4 view = trans.view;

    Shader *shader = shaders::inst()->default;

    //Render Lines
    if(!isConnector)
    {    
        shader->bind();
            shader->setMatrices(trans, model, true, true, true, true);
        
            if(isConnector)
                shader->setf("transparency", 1.0f);
            else
                shader->setf("transparency", params::inst()->boxTransparency);

            m_vboLines->render();

        shader->release();
    }

    //Render Box
    shader = shaders::inst()->box;

    if(shaderSelector == 1)
        shader = shaders::inst()->scan;

    vec4 color;
    switch(colorMode) 
    {
    // Current default
    //case 0:
    //    color = vec4(0.5f, 1.0f, 1.0f, 1.0f);
    //    break;
    //case 1:
    //    color = vec4(0.0f, 0.0f, 1.0f, 1.0f);
    //    break;
    //case 2:
    //    color = vec4(0.0f, 1.0f, 0.0f, 1.0f);
    //    break;
    //case 3:
    //    color = vec4(1.0f, 0.5f, 0.5f, 1.0f);
    //    break;
    //case 4:
    //    color = vec4(1.0f, 0.5f, 0.0f, 1.0f);
    //    break;
    //case 5:
    //    color = vec4(0.0f, 0.5f, 1.0f, 1.0f);
    //    break;
    //case 6:
    //    color = vec4(0.5f, 0.5f, 1.0f, 1.0f);
    //    break;
    //default:
    //    break;

    //Option 1
    case 0:
        color = vec4(1.0f, 0.5f, 0.0f, 1.0f);
        break;
    case 1:
        color = vec4(0.0f, 0.0f, 1.0f, 1.0f);
        break;
    case 2:
        color = vec4(1.0f, 1.0f, 0.0f, 1.0f);
        break;
    case 3:
        color = vec4(1.0f, 0.3f, 0.3f, 1.0f);
        break;
    case 4:
        color = vec4(1.0f, 0.5f, 0.0f, 1.0f);
        break;
    case 5:
        color = vec4(0.0f, 0.5f, 1.0f, 1.0f);
        break;
    case 6:
        color = vec4(0.5f, 0.5f, 1.0f, 1.0f);
        break;
    default:
        break;

    //Jet Cargo
    //case 0:
    //    color = vec4(0.34f, 0.56f, 0.62f, 1.0f);
    //    break;
    //case 1:
    //    color = vec4(0.01f, 0.11f, 0.25f, 1.0f);
    //    break;
    //case 2:
    //    color = vec4(0.03f, 0.21f, 0.44f, 1.0f);
    //    break;
    //case 3:
    //    color = vec4(0.88f, 0.47f, 0.13f, 1.0f);
    //    break;
    //case 4:
    //    color = vec4(1.0f, 0.5f, 0.0f, 1.0f);
    //    break;
    //case 5:
    //    color = vec4(0.64f, 0.32f, 0.13f, 1.0f);
    //    break;
    //case 6:
    //    color = vec4(0.5f, 0.5f, 1.0f, 1.0f);
    //    break;
    //default:
        break;
    }

    if(isConnector)
        color = vec4(1.0f, 0.0f, 0.0f, 1.0f);

    glEnable(GL_CULL_FACE);

    shader->bind();
        shader->setMatrices(trans, model, true, true, true, true);
        shader->set3f("lightPos", params::inst()->lights[0]->position());
        shader->set4f("color", color);
        shader->seti("colorMode", colorMode);
        
        if(isConnector)
            shader->setf("transparency", 1.0f);
        else
            shader->setf("transparency", params::inst()->boxTransparency);

        m_vbo->render();

    shader->release();

    glDisable(GL_CULL_FACE);


}

void Box::renderDepth(const Transform &trans, const mat4 &model, int colorMode, int shaderSelector)
{
    mat4 projection = trans.projection;
    mat4 view = trans.view;

    Shader *shader = shaders::inst()->defaultDepth;

    shader->bind();
        shader->setMatrices(trans, model, true, true, true, true);
        shader->set3f("lightPos", params::inst()->lights[0]->position());
        shader->set4f("color", vec4(1));

        m_vbo->render();

    shader->release();
}

vector<Vertex> Box::vertices()
{
    vector<Vertex> verts;

    for (int i = 0; i < m_faces.size(); ++i)
    {
        Face &f = m_faces[i];

        vec4 color = f.color;
        vec3 normal = f.normal;

        vec3 a = m_vertices[f.a];
        vec3 b = m_vertices[f.b];
        vec3 c = m_vertices[f.c];
        vec3 d = m_vertices[f.d];

        Vertex v0, v1, v2, v3, v4, v5;

        v0.position = a;
        v0.normal = normal;

        v1.position = b;
        v1.normal = normal;

        v2.position = c;
        v2.normal = normal;

        v3.position = c;
        v3.normal = normal;

        v4.position = d;
        v4.normal = normal;

        v5.position = a;
        v5.normal = normal;

        verts.push_back(v0);
        verts.push_back(v1);
        verts.push_back(v2);
        verts.push_back(v3);
        verts.push_back(v4);
        verts.push_back(v5);
    }

    return verts;
}
