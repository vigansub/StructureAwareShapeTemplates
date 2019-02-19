#include "Object.h"
#include "Shader.h"
#include "VertexBufferObject.h"
#include "Mesh.h"
#include "Light.h"
#include <string.h>

Object::Object(const QString &fileName, bool normalize, bool buildLineVBO, bool buildNormalVBO, const vec3 &pos, const vec3 &scale, const vec4 &rot, const vec4 &color)
: m_fileName(fileName),
  m_position(pos),
  m_scale(scale),
  m_rotation(rot),
  m_color(color),
  m_isSelected(false),
  m_nrTriangles(0),
  m_nrVertices(0), 
  m_normalize(normalize), 
  m_lines(buildLineVBO), 
  m_normals(buildNormalVBO), 
  m_vboComplete(nullptr)
{
    prepareData(m_fileName);
    buildVBOMeshComplete(m_vertices);
}

Object::~Object()
{
    for(int i=m_vbosTriangles.size()-1; i>=0; --i)
    {
        VertexBufferObject *vbo = m_vbosTriangles[i];
        delete vbo;
    }

    for(int i=m_vbosLines.size()-1; i>=0; --i)
    {
        VertexBufferObject *vbo = m_vbosLines[i];
        delete vbo;
    }
}

void Object::prepareData(const QString &fileName)
{ 
    QFileInfo fi(fileName);
    QString fn = fi.fileName();
    QString baseName = fileName;
    baseName.replace(fn, "");

    qDebug() << fileName << baseName;
    vector<tinyobj::shape_t> shapes;
    vector<tinyobj::material_t> materials;

    string err;
    //bool ret = tinyobj::LoadObj(shapes, materials, err, m_fileName.toStdString().c_str(), baseName.toStdString().c_str(), true);
    bool ret = tinyobj::LoadObj(shapes, materials, err, fileName.toStdString().c_str(), baseName.toStdString().c_str(), true);

    if(err.length() > 0)
    {
        qDebug() << "Error: " << err.c_str();
    }

    vector<Material> mats;
    if(materials.size() > 0)
    {
        for (size_t i = 0; i < materials.size(); i++) 
        {
            materials[i].name;

            vec3 Ka = vec3(materials[i].ambient[0], materials[i].ambient[1], materials[i].ambient[2]);
            vec3 Kd = vec3(materials[i].diffuse[0], materials[i].diffuse[1], materials[i].diffuse[2]);
            vec3 Ks = vec3(materials[i].specular[0], materials[i].specular[1], materials[i].specular[2]);
            float Ns = materials[i].shininess;
            QString dTexName = baseName + QString(materials[i].diffuse_texname.c_str());
            mats.push_back(Material(Ka, Kd, Ks, Ns, dTexName));
        }
    }

    vector<vector<Vertex>> allVertices;
    vector<vector<uint>> allIndices;

    for (size_t i = 0; i < shapes.size(); i++) 
    {
        vector<uint> &indices = shapes[i].mesh.indices;

        vec3 p, n;
        vec2 t;
        vector<Vertex> vertices;

        for(size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++) 
        {
            p = vec3(shapes[i].mesh.positions[3*v+0], shapes[i].mesh.positions[3*v+1], shapes[i].mesh.positions[3*v+2]);

            int n1 = 3*v+0, n2 = 3*v+1, n3 = 3*v+2;
            int t1 = 2*v+0, t2 = 2*v+1;

            if(n1 < shapes[i].mesh.normals.size() && n2 < shapes[i].mesh.normals.size() && n3 < shapes[i].mesh.normals.size())
            {
                n = vec3(shapes[i].mesh.normals[n1], shapes[i].mesh.normals[n2], shapes[i].mesh.normals[n3]);
            }

            if(t1 < shapes[i].mesh.texcoords.size() && t2 < shapes[i].mesh.texcoords.size())
            {
                t = vec2(shapes[i].mesh.texcoords[t1], shapes[i].mesh.texcoords[t2]);
            }

            vertices.push_back(Vertex(p, n, vec4(), t));            
        }

        Material defMat(vec3(0.0f), vec3(0.5, 0.5, 0.5), vec3(0.2), 10);

        if(shapes[i].mesh.material_ids.size() > 0 && shapes[i].mesh.material_ids[0] != -1)
        {
            int matId = shapes[i].mesh.material_ids[0];
            m_materials.push_back(mats[matId]);
        }
        else
        {
            m_materials.push_back(defMat);
        }

        allIndices.push_back(indices);
        allVertices.push_back(vertices);
    }

    m_allIndices = allIndices;
    m_allVertices = allVertices;

    if(m_normalize)
    {
        normalizeGeometry(allVertices, m_position, m_scale, m_rotation);
        //normalizeGeometry(allVertices, vec3(0.0f), vec3(1.0f), vec4(0.0f));
    }

    for(int i=0; i<allVertices.size(); ++i)
    {
        buildVBOMesh(allVertices[i], allIndices[i]);

        if(m_lines)
        {
            buildVBOLines(allVertices[i], allIndices[i]);
        }

        if(m_normals)
        {
            buildVBONormals(allVertices[i], allIndices[i]);
        }
    }

    for (int i = 0; i < allIndices.size(); ++i)
    {
        vector<uint> indices = allIndices[i];
        vector<Vertex> vertices = allVertices[i];

        for (int j = 0; j < indices.size(); ++j)
        {
            Vertex v = vertices[indices[j]];
            m_vertices.push_back(v);
        }
    }

}

void Object::render(const Transform &trans)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);

    mat4 model = mat4::identitiy();
    if(!m_normalize)
    {
        model = mat4::translate(m_position) * mat4::rotateY(m_rotation.y) * mat4::scale(m_scale); 
    }

	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);
    glEnable(GL_CLIP_DISTANCE0);    

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(params::inst()->polygonOffsetFactor, params::inst()->polygonOffsetUnits);

    //glDepthFunc(GL_LEQUAL);
    //glDepthRange(params::inst()->depthRangeMin, params::inst()->depthRangeMax);   

    if (params::inst()->renderObjects)
    {
        Shader *shader = shaders::inst()->object;
        shader->bind();

        shader->setMatrices(trans, model, true, true, true, true);
        shader->set3f("camPos", params::inst()->camPos);
        shader->seti("applyShadow", params::inst()->applyShadow);
        shader->setf("shadowIntensity", params::inst()->shadowIntensity);
        shader->seti("isSelected", m_isSelected);
        shader->set4f("clipPlane", params::inst()->clipPlaneGround);
        shader->seti("colorMode", 0);

        shader->setLights(params::inst()->lights);

        for (uint i = 0; i < m_vbosTriangles.size(); ++i)
        {
            shader->setMaterial(m_materials[i]);
            m_vbosTriangles[i]->render();
        }

        shader->release();

//        if (m_lines && params::inst()->renderWireframe)
//        {
//            shader = shaders::inst()->objectLines;
//            shader->bind();
//
//            shader->setMatrices(trans, model, true, true, true, false);
//
//            for (uint i = 0; i < m_vbosLines.size(); ++i)
//            {
//                m_vbosLines[i]->render();
//            }
//
//            shader->release();
//        }
//
//        if (m_normals && params::inst()->renderNormals)
//        {
//#ifdef WIN32
//            shader = shaders::inst()->default;
//#else
//            shader = shaders::inst()->default_;
//#endif
//            shader->bind();
//
//            shader->setMatrices(trans, model, true, true, true, false);
//
//            for (uint i = 0; i < m_vbosNormals.size(); ++i)
//            {
//                m_vbosNormals[i]->render();
//            }
//
//            shader->release();
//        }
    }

    if (params::inst()->renderGenerations)
    {
        glPointSize(4.0f);
        Shader *shader = shaders::inst()->object;
        shader->bind();

        shader->setMatrices(trans, model, true, true, true, true);
        shader->set3f("camPos", params::inst()->camPos);
        shader->seti("applyShadow", params::inst()->applyShadow);
        shader->setf("shadowIntensity", params::inst()->shadowIntensity);
        shader->seti("isSelected", m_isSelected);
        shader->set4f("clipPlane", params::inst()->clipPlaneGround);
        shader->seti("colorMode", 0);

        shader->setLights(params::inst()->lights);

        shader->setMatrices(trans, model, true, true, true, false);
        shader->set3f("lightPos", params::inst()->lights[0]->position());

        for (uint i = 0; i < m_vbosTriangles_skin.size(); ++i)
        {
            shader->setMaterial(m_materials[i]);
            m_vbosTriangles_skin[i]->render();
        }

        //m_vboComplete->render();
        shader->release();
    }

	glDisable(GL_CULL_FACE);    

	glPopClientAttrib();
	glPopAttrib();
}

void Object::renderScan(const Transform &trans)
{
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);

    mat4 model = mat4::identitiy();
    if (!m_normalize)
    {
        model = mat4::translate(m_position) * mat4::rotateY(m_rotation.y) * mat4::scale(m_scale);
    }

	mat4 view = trans.view;
	mat4 projection = trans.projection;

	Shader *shaderScan = shaders::inst()->scan;
	shaderScan->bind();

	shaderScan->setMatrix("matModel", model, GL_TRUE);
	shaderScan->setMatrix("matView", view, GL_TRUE);
	shaderScan->setMatrix("matProjection", projection, GL_TRUE);


	for (uint i = 0; i<m_vbosTriangles.size(); ++i)
	{
		m_vbosTriangles[i]->render();
	}

	shaderScan->release();

	glPopClientAttrib();
	glPopAttrib();
}

void Object::renderDepth(const Transform &trans)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);

    mat4 model = mat4::identitiy();

    if (params::inst()->renderObjects)
    {
        glEnable(GL_CLIP_DISTANCE0);

        
        if (!m_normalize)
        {
            model = mat4::translate(m_position) * mat4::rotateY(m_rotation.y) * mat4::scale(m_scale);
        }

        glDisable(GL_CULL_FACE);
        glCullFace(GL_FRONT);

        glClearDepth(1.0);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glDepthFunc(GL_LEQUAL);

        glDepthRange(params::inst()->depthRangeMin, params::inst()->depthRangeMax);
        glPolygonOffset(params::inst()->polygonOffsetFactor, params::inst()->polygonOffsetUnits);

        Shader *shader = shaders::inst()->objectDepth;
        shader->bind();

        shader->setMatrices(trans, model, true, true, true, true);
        shader->set4f("clipPlane", params::inst()->clipPlaneGround);

        for (uint i = 0; i < m_vbosTriangles.size(); ++i)
        {
            m_vbosTriangles[i]->render();
        }

        shader->release();
    }

    if (params::inst()->renderGenerations)
    {
        Shader *shader = shaders::inst()->defaultDepth;
        shader->bind();
        shader->setMatrices(trans, model, true, true, true, false);
        //m_vboComplete->render();

        for (uint i = 0; i < m_vbosTriangles_skin.size(); ++i)
        {
            m_vbosTriangles_skin[i]->render();
        }

        shader->release();
    }

	glPopClientAttrib();
	glPopAttrib();
}

void Object::buildVBOMesh(vector<Vertex> &vertices, vector<uint> &indices)
{    
    vector<Vertex> tmpVertices;
    for (int i = 0; i < indices.size(); ++i)
    {
        Vertex v = vertices[indices[i]];
        tmpVertices.push_back(v);
    }

    vector<vec3> normals;
    for (uint i = 0; i < tmpVertices.size(); i += 3)
    {
        vec3 a = tmpVertices[i].position;
        vec3 b = tmpVertices[i + 1].position;
        vec3 c = tmpVertices[i + 2].position;

        vec3 s = a - b;
        vec3 t = a - c;

        vec3 n = normalize(cross(s, t));

        normals.push_back(n);
        normals.push_back(n);
        normals.push_back(n);
    }

    VertexBufferObject::DATA *data = new VertexBufferObject::DATA[tmpVertices.size()];

	for (uint i = 0; i<tmpVertices.size(); ++i)
	{
        vec3 v = tmpVertices[i].position;
        vec3 n = normals[i];// normalize(vertices[i].normal);
        vec2 t = tmpVertices[i].texture;

		data[i].vx = v.x;
		data[i].vy = v.y;
		data[i].vz = v.z;
		data[i].vw = 0.0f;

        data[i].cx = m_color.x;
        data[i].cy = m_color.y;
        data[i].cz = m_color.z;
		data[i].cw = m_color.w;
			
        data[i].nx = n.x;
		data[i].ny = n.y;
		data[i].nz = n.z;
		data[i].nw = 0.0f;
            
		data[i].tx = t.x;
		data[i].ty = t.y;
        data[i].tz = 0.0f;
        data[i].tw = 0.0f;
	}

	VertexBufferObject* vbo = new VertexBufferObject();
    vbo->setData(data, GL_STATIC_DRAW, tmpVertices.size(), GL_TRIANGLES);
    //vbo->setIndexData(indices.data(), GL_STATIC_DRAW, indices.size());
    vbo->bindDefaultAttribs();

    m_vbosTriangles.push_back(vbo); 

	delete[] data;    
}

void Object::buildVBOMesh_skin(vector<Vertex> &vertices, vector<uint> &indices)
{
    vector<Vertex> tmpVertices;
    for (int i = 0; i < indices.size(); ++i)
    {
        Vertex v = vertices[indices[i]];
        tmpVertices.push_back(v);
    }

    vector<vec3> normals;
    for (uint i = 0; i < tmpVertices.size(); i += 3)
    {
        vec3 a = tmpVertices[i].position;
        vec3 b = tmpVertices[i + 1].position;
        vec3 c = tmpVertices[i + 2].position;

        vec3 s = a - b;
        vec3 t = a - c;

        vec3 n = normalize(cross(s, t));

        normals.push_back(n);
        normals.push_back(n);
        normals.push_back(n);
    }

    VertexBufferObject::DATA *data = new VertexBufferObject::DATA[tmpVertices.size()];

    for (uint i = 0; i<tmpVertices.size(); ++i)
    {
        vec3 v = tmpVertices[i].position;
        vec3 n = normals[i];// normalize(vertices[i].normal);
        vec2 t = tmpVertices[i].texture;
        vec4 c = tmpVertices[i].color;

        data[i].vx = v.x;
        data[i].vy = v.y;
        data[i].vz = v.z;
        data[i].vw = 0.0f;

        data[i].cx = c.x;
        data[i].cy = c.y;
        data[i].cz = c.z;
        data[i].cw = c.w;

        data[i].nx = n.x;
        data[i].ny = n.y;
        data[i].nz = n.z;
        data[i].nw = 0.0f;

        data[i].tx = t.x;
        data[i].ty = t.y;
        data[i].tz = 0.0f;
        data[i].tw = 0.0f;
    }

    VertexBufferObject* vbo = new VertexBufferObject();
    vbo->setData(data, GL_STATIC_DRAW, tmpVertices.size(), GL_TRIANGLES);
    //vbo->setIndexData(indices.data(), GL_STATIC_DRAW, indices.size());
    vbo->bindDefaultAttribs();

    m_vbosTriangles_skin.push_back(vbo);

    delete[] data;
}

void Object::buildVBOLines(vector<Vertex> &vertices, vector<uint> &indices)
{
    vector<Vertex> tmp;
	for (uint i = 0; i<indices.size()-3; i+=3)
	{
        Vertex &a = vertices[indices[i]];
        Vertex &b = vertices[indices[i+1]];
        Vertex &c = vertices[indices[i+2]];

        tmp.push_back(a);
        tmp.push_back(b);

        tmp.push_back(a);
        tmp.push_back(c);

        tmp.push_back(b);
        tmp.push_back(c);
    }    

    VertexBufferObject::DATA *data = new VertexBufferObject::DATA[tmp.size()];

	for (uint i = 0; i<tmp.size(); ++i)
	{
		vec3 v = tmp[i].position;
		vec3 n = tmp[i].normal;
		vec2 t = tmp[i].texture;

		data[i].vx = v.x;
		data[i].vy = v.y;
		data[i].vz = v.z;
		data[i].vw = 0.0f;

        data[i].cx = m_color.x;
        data[i].cy = m_color.y;
        data[i].cz = m_color.z;
		data[i].cw = m_color.w;
			
        data[i].nx = n.x;
		data[i].ny = n.y;
		data[i].nz = n.z;
		data[i].nw = 0.0f;
            
		data[i].tx = t.x;
		data[i].ty = t.y;
        data[i].tz = 0.0f;
        data[i].tw = 0.0f;
	}

	VertexBufferObject* vbo = new VertexBufferObject();
	vbo->setData(data, GL_STATIC_DRAW, tmp.size(), GL_LINES);
    vbo->bindDefaultAttribs();
    m_vbosLines.push_back(vbo);

	delete[] data;    
}

void Object::buildVBOLines_skin(vector<Vertex> &vertices, vector<uint> &indices)
{
    vector<Vertex> tmp;
    for (uint i = 0; i<indices.size() - 3; i += 3)
    {
        Vertex &a = vertices[indices[i]];
        Vertex &b = vertices[indices[i + 1]];
        Vertex &c = vertices[indices[i + 2]];

        tmp.push_back(a);
        tmp.push_back(b);

        tmp.push_back(a);
        tmp.push_back(c);

        tmp.push_back(b);
        tmp.push_back(c);
    }

    VertexBufferObject::DATA *data = new VertexBufferObject::DATA[tmp.size()];

    for (uint i = 0; i<tmp.size(); ++i)
    {
        vec3 v = tmp[i].position;
        vec3 n = tmp[i].normal;
        vec2 t = tmp[i].texture;

        data[i].vx = v.x;
        data[i].vy = v.y;
        data[i].vz = v.z;
        data[i].vw = 0.0f;

        data[i].cx = m_color.x;
        data[i].cy = m_color.y;
        data[i].cz = m_color.z;
        data[i].cw = m_color.w;

        data[i].nx = n.x;
        data[i].ny = n.y;
        data[i].nz = n.z;
        data[i].nw = 0.0f;

        data[i].tx = t.x;
        data[i].ty = t.y;
        data[i].tz = 0.0f;
        data[i].tw = 0.0f;
    }

    VertexBufferObject* vbo = new VertexBufferObject();
    vbo->setData(data, GL_STATIC_DRAW, tmp.size(), GL_LINES);
    vbo->bindDefaultAttribs();
    m_vbosLines_skin.push_back(vbo);

    delete[] data;
}

void Object::buildVBONormals(vector<Vertex> &vertices, vector<uint> &indices)
{
    float s = 0.1;
    vector<vec3> positions;
	for (uint i = 0; i<indices.size()-3; i+=3)
	{
        Vertex &a = vertices[indices[i]];
        Vertex &b = vertices[indices[i+1]];
        Vertex &c = vertices[indices[i+2]];

        positions.push_back(a.position);
        positions.push_back(a.position + a.normal * s);

        positions.push_back(b.position);
        positions.push_back(b.position + b.normal * s);

        positions.push_back(c.position);
        positions.push_back(c.position + c.normal * s);
    }    

    VertexBufferObject::DATA *data = new VertexBufferObject::DATA[positions.size()];

	for (uint i = 0; i<positions.size(); ++i)
	{
		vec3 v = positions[i];

		data[i].vx = v.x;
		data[i].vy = v.y;
		data[i].vz = v.z;
		data[i].vw = 0.0f;

        data[i].cx = v.x;
        data[i].cy = v.y;
        data[i].cz = v.z;
		data[i].cw = 1.0f;
			
        data[i].nx = 0.0f;
		data[i].ny = 0.0f;
		data[i].nz = 0.0f;
		data[i].nw = 0.0f;
            
		data[i].tx = 0.0f;
		data[i].ty = 0.0f;
        data[i].tz = 0.0f;
        data[i].tw = 0.0f;
	}

	VertexBufferObject* vbo = new VertexBufferObject();
	vbo->setData(data, GL_STATIC_DRAW, positions.size(), GL_LINES);
    vbo->bindDefaultAttribs();
    m_vbosNormals.push_back(vbo);

	delete[] data;    
}

void Object::buildVBONormals_skin(vector<Vertex> &vertices, vector<uint> &indices)
{
    float s = 0.1;
    vector<vec3> positions;
    for (uint i = 0; i<indices.size() - 3; i += 3)
    {
        Vertex &a = vertices[indices[i]];
        Vertex &b = vertices[indices[i + 1]];
        Vertex &c = vertices[indices[i + 2]];

        positions.push_back(a.position);
        positions.push_back(a.position + a.normal * s);

        positions.push_back(b.position);
        positions.push_back(b.position + b.normal * s);

        positions.push_back(c.position);
        positions.push_back(c.position + c.normal * s);
    }

    VertexBufferObject::DATA *data = new VertexBufferObject::DATA[positions.size()];

    for (uint i = 0; i<positions.size(); ++i)
    {
        vec3 v = positions[i];

        data[i].vx = v.x;
        data[i].vy = v.y;
        data[i].vz = v.z;
        data[i].vw = 0.0f;

        data[i].cx = v.x;
        data[i].cy = v.y;
        data[i].cz = v.z;
        data[i].cw = 1.0f;

        data[i].nx = 0.0f;
        data[i].ny = 0.0f;
        data[i].nz = 0.0f;
        data[i].nw = 0.0f;

        data[i].tx = 0.0f;
        data[i].ty = 0.0f;
        data[i].tz = 0.0f;
        data[i].tw = 0.0f;
    }

    VertexBufferObject* vbo = new VertexBufferObject();
    vbo->setData(data, GL_STATIC_DRAW, positions.size(), GL_LINES);
    vbo->bindDefaultAttribs();
    m_vbosNormals_skin.push_back(vbo);

    delete[] data;
}

void Object::normalizeGeometry(vector<vector<Vertex>> &vertices, const vec3 &translate, const vec3 &scale, const vec4 &rotate)
{
	vec3 mi = vec3(math_maxfloat, math_maxfloat, math_maxfloat);
	vec3 ma = vec3(math_minfloat, math_minfloat, math_minfloat);

	for (int i = 0; i<vertices.size(); ++i)
	{
        for(int j=0; j<vertices[i].size(); ++j)
        {
		    vec3 &a = vertices[i][j].position;

		    if (a.x > ma.x) ma.x = a.x;
		    if (a.y > ma.y) ma.y = a.y;
		    if (a.z > ma.z) ma.z = a.z;

		    if (a.x < mi.x) mi.x = a.x;
		    if (a.y < mi.y) mi.y = a.y;
		    if (a.z < mi.z) mi.z = a.z;
        }
	}    

	vec3 d = ma - mi;
	float s = max(d.x, max(d.y, d.z));

	vec3 shift = d / s /2;
	for (int i = 0; i<vertices.size(); ++i)
	{
        for(int j=0; j<vertices[i].size(); ++j)
        {
		    vec3 &a = vertices[i][j].position;

		    a -= mi;
		    a /= s;
		    a -= vec3(shift.x, 0.0f, shift.z);

		    mat4 m = mat4::identitiy();
		    m *= mat4::translate(translate);
		    m *= mat4::rotate(rotate.x, vec3(rotate.y, rotate.z, rotate.w));
		    m *= mat4::scale(scale);

		    vec4 ta = m * vec4(a);
		    vertices[i][j].position = vec3(ta.x, ta.y, ta.z);
        }
	}

	mi = vec3(math_maxfloat, math_maxfloat, math_maxfloat);
	ma = vec3(math_minfloat, math_minfloat, math_minfloat);

	for (int i = 0; i<vertices.size(); ++i)
	{
        for(int j=0; j<vertices[i].size(); ++j)
        {
		    vec3 &a = vertices[i][j].position;

		    if (a.x > ma.x) ma.x = a.x;
		    if (a.y > ma.y) ma.y = a.y;
		    if (a.z > ma.z) ma.z = a.z;

		    if (a.x < mi.x) mi.x = a.x;
		    if (a.y < mi.y) mi.y = a.y;
		    if (a.z < mi.z) mi.z = a.z;
        }
	} 

    m_bb = BoundingBox(mi, ma);
}

float Object::selected(Picking &pick, const Transform &trans, int sw, int sh, int mx, int my)
{     
    mat4 model = mat4::identitiy();        

    if(!m_normalize)
    {
        model = mat4::translate(m_position) * mat4::rotateY(m_rotation.y) * mat4::scale(m_scale);
    }
        
    return pick.select(trans, model, m_bb.mi(), m_bb.ma(), sw, sh, mx, my);  
}

void Object::move(int x, int y, const vec3 &dir, const vec3 &up, const vec3 &right, const vec3 &pos)
{
	vec3 movZ, movX, movY;	
	vec3 onPlane = cross(right, vec3(0, 1, 0));
    vec3 p = m_position;

    movY = vec3(0, p.y, 0);
    movX = vec3(p.x, 0, 0) + right * x * 0.1;
    movZ = vec3(0, 0, p.z) + onPlane * y * 0.1;
        	
    vec3 result = movX + movY + movZ ;                    
    m_position = result;    
}


void Object::buildVBOMeshComplete(vector<Vertex> vertices)
{
    //delete m_vboComplete;

    for (int i = 0; i<m_allVertices_skin.size(); ++i)
    {
        buildVBOMesh_skin(m_allVertices_skin[i], m_allIndices[i]);

        if (m_lines)
        {
            buildVBOLines_skin(m_allVertices_skin[i], m_allIndices[i]);
        }

        if (m_normals)
        {
            buildVBONormals_skin(m_allVertices_skin[i], m_allIndices[i]);
        }
    }

    //uint nrVertices = vertices.size();
    //VertexBufferObject::DATA *attrData = new VertexBufferObject::DATA[nrVertices];
    //
    //vector<vec3> normals;
    //for (uint i = 0; i < nrVertices; i+=3)
    //{
    //    vec3 a = vertices[i].position;
    //    vec3 b = vertices[i+1].position;
    //    vec3 c = vertices[i+2].position;
    //
    //    vec3 s = a - b;
    //    vec3 t = a - c;
    //
    //    vec3 n = normalize(cross(s, t));
    //
    //    normals.push_back(n);
    //    normals.push_back(n);
    //    normals.push_back(n);
    //}
    //
    //for (uint i = 0; i<nrVertices; ++i)
    //{
    //    vec3 v = vertices[i].position;
    //    vec3 n = normals[i]; //vertices[i].normal;
    //    vec4 c = vertices[i].color;
    //    c.w = 1.0f;
    //
    //    attrData[i].vx = v.x;
    //    attrData[i].vy = v.y;
    //    attrData[i].vz = v.z;
    //    attrData[i].vw = 1.0f;
    //
    //    attrData[i].nx = n.x;
    //    attrData[i].ny = n.y;
    //    attrData[i].nz = n.z;
    //    attrData[i].nw = 0.0f;
    //
    //    attrData[i].cx = c.x;
    //    attrData[i].cy = c.y;
    //    attrData[i].cz = c.z;
    //    attrData[i].cw = 1.0f;
    //
    //    attrData[i].tx = 0.0f;
    //    attrData[i].ty = 0.0f;
    //    attrData[i].tz = 0.0f;
    //    attrData[i].tw = 0.0f;
    //
    //}
    //
    //m_vboComplete = new VertexBufferObject();
    //m_vboComplete->setData(attrData, GL_STATIC_DRAW, nrVertices, GL_TRIANGLES);
    //m_vboComplete->bindDefaultAttribs();
    //
    //delete[] attrData;    
}

vector<Vertex> Object::transform(Proc_Model p, Proc_Model q)
{
    vector<int> box_list;
    vector<vec3> pos1 = p.position(), pos2 = q.position(), size1 = p.size(), size2 = q.size();
    vector<vector<double>> box_distances;
    box_list.resize(m_vertices.size());
    box_distances.resize(m_vertices.size());
    int mode = 3;
    if (mode == 1)
    {
        //Weight is given by partial sums
        for (int i = 0; i < m_vertices.size(); ++i)
        {
            box_distances[i].resize(pos1.size());
            bool zero_exists = false;
            int zero_val;
            double prod_val = 1.0;
            vector<double> vals;
            vals.resize(pos1.size());
            for (int j = 0; j < p.position().size(); ++j)
            {
                vals[j] = p.project_point_to_box(m_vertices[i].position, pos1[j], size1[j]);
                prod_val *= vals[j];
                if (vals[j] == 0)
                {
                    zero_exists = true;
                    zero_val = j;
                }
            }

            if (zero_exists)
            {
                box_list[i] = zero_val;
                for (int j = 0; j < box_distances[i].size(); ++j)
                {
                    if (j == zero_val)
                        box_distances[i][j] = 1;
                    else box_distances[i][j] = 0;
                }
            }
            else
            {
                double sum_val = 0.0;
                box_list[i] = 0;
                for (int j = 0; j < box_distances[i].size(); ++j)
                {
                    box_distances[i][j] = pow((prod_val / vals[j]), 2);
                    sum_val += box_distances[i][j];
                    if (vals[box_list[i]] < vals[j])
                        box_list[i] = j;
                }

                for (int j = 0; j < box_distances[i].size(); ++j)
                    box_distances[i][j] /= sum_val;
            }
        }
    }
    else if (mode == 2)
    {
        //Weight is given by exp(-d^2)/d
        for (int i = 0; i < m_vertices.size(); ++i)
        {
            box_distances[i].resize(pos1.size());
            bool zero_exists = false;
            int zero_val;
            vector<double> vals;
            vals.resize(pos1.size());
            for (int j = 0; j < p.position().size(); ++j)
            {
                vals[j] = p.project_point_to_box(m_vertices[i].position, pos1[j], size1[j]);
                if (vals[j] == 0)
                {
                    zero_exists = true;
                    zero_val = j;
                }
            }

            if (zero_exists)
            {
                box_list[i] = zero_val;
                for (int j = 0; j < box_distances[i].size(); ++j)
                {
                    if (j == zero_val)
                        box_distances[i][j] = 1;
                    else box_distances[i][j] = 0;
                }
            }
            else
            {
                double sum_val = 0.0;
                box_list[i] = 0;
                for (int j = 0; j < box_distances[i].size(); ++j)
                {
                    box_distances[i][j] = exp(-vals[j] * vals[j]) / (pow(vals[j], 0.5));
                    sum_val += box_distances[i][j];
                    if (vals[box_list[i]] < vals[j])
                        box_list[i] = j;
                }

                for (int j = 0; j < box_distances[i].size(); ++j)
                    box_distances[i][j] /= sum_val;
            }
        }
    }
    else if (mode == 3)
    {
        //Square Exponential Weights
        for (int i = 0; i < m_vertices.size(); ++i)
        {
            box_distances[i].resize(pos1.size());
            double val = p.project_point_to_box(m_vertices[i].position, pos1[0], size1[0]);
            box_list[i] = 0;
            box_distances[i][0] = exp(-val*val);
            double sum_val = exp(-val*val);
            for (int j = 1; j < p.position().size(); ++j)
            {
                double k = p.project_point_to_box(m_vertices[i].position, pos1[j], size1[j]);
                if (k < val)
                {
                    val = k;
                    box_list[i] = j;
                }
                box_distances[i][j] = exp(-k*k);
                sum_val += exp(-k*k);
            }
            for (int j = 0; j < box_distances[i].size(); ++j)
                box_distances[i][j] /= sum_val;
        }

    }
    else if (mode == 4)
    {
        // Plain Affine Transformation
        vector<Vertex> vertices = m_vertices;
        for (int i = 0; i < m_vertices.size(); ++i)
        {
            int corr_box = box_list[i];
            vec3 new_point = (vertices[i].position - pos1[corr_box]);
            new_point = vec3(new_point.x*size2[corr_box].x / size1[corr_box].x, new_point.y*size2[corr_box].y / size1[corr_box].y, new_point.z*size2[corr_box].z / size1[corr_box].z);
            vertices[i].position = new_point + pos2[corr_box];
            if (corr_box == 0)
                vertices[i].color = vec4(1.0, 0.0, 0.0, 1.0);
            if (corr_box >= 1 && corr_box <= 2)
                vertices[i].color = vec4(0.0, 1.0, 0.0, 1.0);
            if (corr_box >= 3 && corr_box <= 6)
                vertices[i].color = vec4(0.0, 0.0, 1.0, 1.0);
            if (corr_box >= 7)
                vertices[i].color = vec4(1.0, 1.0, 1.0, 1.0);
        }
    }

    vector<Vertex> vertices = m_vertices;
    vector<int> colorMode = p.colorMode();
    //Exponential distance weights
    for (int i = 0; i < m_vertices.size(); ++i)
    {
        vertices[i].position -= m_vertices[i].position;
        for (int j = 0; j < pos1.size(); ++j)
        {
            vec3 new_point = (m_vertices[i].position - pos1[j]);
            new_point = vec3(new_point.x*size2[j].x / size1[j].x, new_point.y*size2[j].y / size1[j].y, new_point.z*size2[j].z / size1[j].z);
            vertices[i].position += box_distances[i][j]*(new_point + pos2[j]);
        }
    }

    if (mode == 3)
    {
        vector<vec3> qbox = q.bounding_box();
        vector<vec3> vbox = q.bounding_box(vertices);
        for (int i = 0; i < vertices.size(); ++i)
        {
            vertices[i].position = vec3((vertices[i].position.x - vbox[0].x)*(qbox[1].x - qbox[0].x) / (vbox[1].x - vbox[0].x) + qbox[0].x, (vertices[i].position.y - vbox[0].y)*(qbox[1].y - qbox[0].y) / (vbox[1].y - vbox[0].y) + qbox[0].y, (vertices[i].position.z - vbox[0].z)*(qbox[1].z - qbox[0].z) / (vbox[1].z - vbox[0].z) + qbox[0].z);
            double val = (p.project_point_to_box(m_vertices[i].position, pos1[0], size1[0]));
            box_list[i] = 0;
            for (int j = 1; j < p.position().size(); ++j)
            {
                double k = p.project_point_to_box(m_vertices[i].position, pos1[j], size1[j]);
                if (k < val)
                {
                    val = k;
                    box_list[i] = j;
                }
            }

            int corr_box = box_list[i];
            
            switch (colorMode[corr_box])
            {
            case 0:
                vertices[i].color = vec4(1.0f, 0.6f, 0.3f, 1.0f);
                break;
            case 3:
                //vertices[i].color = vec4(1.0f, 1.0f, 0.3f, 1.0f);
                vertices[i].color = vec4(0.3f, 0.3f, 1.0f, 1.0f);
                break;
            case 1:
                vertices[i].color = vec4(1.0f, 1.0f, 0.3f, 1.0f);
                break;
            case 2:
                vertices[i].color = vec4(1.0f, 0.3f, 0.3f, 1.0f);
                break;
            case 4:
                vertices[i].color = vec4(0.3f, 1.0f, 0.0f, 1.0f);
                break;
            case 5:
                vertices[i].color = vec4(0.0f, 0.5f, 1.0f, 1.0f);
                break;
            case 6:
                vertices[i].color = vec4(0.5f, 0.5f, 1.0f, 1.0f);
                break;
            default:
                break;
            }
        }
    }

    m_allVertices_skin = split_vertices(vertices);

    return vertices;
}

vector<vector<Vertex>> Object::split_vertices(vector<Vertex> vertices)
{
    vector<vector<Vertex>> allVertices_new;
    allVertices_new.resize(m_allVertices.size());
    int k = 0;
    for (int i = 0; i < m_allIndices.size(); ++i)
    {
        vector<uint> indices = m_allIndices[i];
        vector<Vertex> vertices_old = m_allVertices[i];
        allVertices_new[i].resize(vertices_old.size());

        for (int j = 0; j < indices.size(); ++j)
        {
            Vertex v = vertices[k];
            k++;
            allVertices_new[i][indices[j]] = v;
        }
    }

    return allVertices_new;
}

void Object::saveObj(QString filename)
{
    saveObj(filename, m_vertices);
}

void Object::saveObj(QString fileName, vector<Vertex> &vertices)
{
    
    vector<Face> faces;
    for (int i = 0; i<vertices.size() - 2; i += 3)
    {
        Face f;
        f.a = i + 0 + 1;
        f.b = i + 1 + 1;
        f.c = i + 2 + 1;

        faces.push_back(f);
    }

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out(&file);

    for (int i = 0; i<vertices.size(); ++i)
    {
        vec3 p = vertices[i].position;
        out << "v " << p.x << " " << p.y << " " << p.z << endl;
    }

    for (int i = 0; i<vertices.size(); ++i)
    {
        vec4 c = vertices[i].color;
        out << "vc " << c.x << " " << c.y << " " << c.z << endl;
    }

    for (int i = 0; i<faces.size(); ++i)
    {
        Face f = faces[i];

        out << "f " << f.a << "/" << f.a << "/" << f.a << " " << f.b << "/" << f.b << "/" << f.b
            << " " << f.c << "/" << f.c << "/" << f.c << endl;
    }

    file.close();
}

void Object::saveObjPC(QString fileName, vector<Vertex> &vertices)
{


    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out(&file);

    for (int i = 0; i<vertices.size(); ++i)
    {
        vec3 p = vertices[i].position;
        out << "v " << p.x << " " << p.y << " " << p.z << endl;
    }
    file.close();
}

void Object::loadOFF(QString fileName, vector<double> feature)
{
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;

    m_vertices.clear();
    
    QTextStream in(&file);
    while (!file.atEnd())
    {
        QString line = file.readLine();
        QStringList k = line.split(' ');
        for (int i = 0; i < k.length(); ++i)
        {
            if ((k.at(0)).toStdString().compare("v"))
            {
                Vertex v;
            }
        }
    }
}