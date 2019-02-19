#ifndef OBJECT_H
#define OBJECT_H

#include "Headers.h"
#include "Geometry.h"
#include "ObjLoader.h"
#include "Mesh.h"
#include "Proc_Model.h"

class Shader;
class VertexBufferObject;

class Object
{
public:
    struct Face
    {
        uint a, b, c;
    };

    Object(const QString &fileName, bool normalize = true, bool buildLineVBO = false, bool buildNormalVBO = false, const vec3 &pos = vec3(), const vec3 &scale = vec3(1, 1, 1), const vec4 &rot = vec4(), const vec4 &color = vec4(1, 1, 1, 1));
    ~Object();

    void render(const Transform &trans);
	void renderScan(const Transform &trans);
    void renderDepth(const Transform &trans);

    bool m_isSelected;
    vec3 m_position;
    vec4 m_rotation;
    vec3 m_scale;
    vec4 m_color;
    bool m_normalize;
    bool m_lines;
    bool m_normals;

    vector<Vertex> m_vertices;
    vector<vector<Vertex>> m_allVertices, m_allVertices_skin;
    vector<vector<uint>> m_allIndices;

    BoundingBox m_bb;
    void move(int x, int y, const vec3 &dir, const vec3 &up, const vec3 &right, const vec3 &pos);

private:
    //void prepareData(const QString &fileName);
    void buildVBOMesh(vector<Vertex> &vertices, vector<uint> &indices);
    void buildVBOLines(vector<Vertex> &vertices, vector<uint> &indices);
    void buildVBONormals(vector<Vertex> &vertices, vector<uint> &indices);
    void buildVBOMesh_skin(vector<Vertex> &vertices, vector<uint> &indices);
    void buildVBOLines_skin(vector<Vertex> &vertices, vector<uint> &indices);
    void buildVBONormals_skin(vector<Vertex> &vertices, vector<uint> &indices);
    void normalizeGeometry(vector<vector<Vertex>> &vertices, const vec3 &translate, const vec3 &scale, const vec4 &rotate);        

    

public:
    void prepareData(const QString &fileName);
    QString m_fileName;

    vector<VertexBufferObject *> m_vbosTriangles;
    vector<VertexBufferObject *> m_vbosLines;
    vector<VertexBufferObject *> m_vbosNormals;

    vector<VertexBufferObject *> m_vbosTriangles_skin;
    vector<VertexBufferObject *> m_vbosLines_skin;
    vector<VertexBufferObject *> m_vbosNormals_skin;

    VertexBufferObject *m_vboComplete;

    int m_nrTriangles;
    int m_nrVertices;    

    vector<Material> m_materials;

    float selected(Picking &pick, const Transform &trans, int sw, int sh, int mx, int my);

    vector<Vertex> transform(Proc_Model p, Proc_Model q);
    vector<vector<Vertex>> split_vertices(vector<Vertex> vertices);
    void buildVBOMeshComplete(vector<Vertex> vertices);

    void saveObj(QString fileName);
    void saveObj(QString fileName, vector<Vertex> &vertices);
    void saveObjPC(QString fileName, vector<Vertex> &vertices);
    void loadOFF(QString fileName, vector<double> feature);
};

#endif