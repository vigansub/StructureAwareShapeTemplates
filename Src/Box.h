#ifndef BOX_H
#define BOX_H

#include "Headers.h"
#include "VertexBufferObject.h"

class Shader;
class VertexBufferObject;

class Box
{
public:
    enum Type
    {
        FRONT = 0,
        BACK,
        LEFT,
        RIGHT,
        TOP,
        BOTTOM
    };

    struct Face
    {
        uint a, b, c, d;
        vec4 color;
        vec3 normal;
        Type type;
    };

    Box();
    ~Box();

    void init();
    void buildVBO();
    void buildVBOLines();
    void render(const Transform &trans, const mat4 &model = mat4::identitiy(), int colorMode = 0, int shaderSelector = 0, bool isConnector = false);
    void renderDepth(const Transform &trans, const mat4 &model = mat4::identitiy(), int colorMode = 0, int shaderSelector = 0);
    vector<Vertex> vertices();

private:
    float m_size;

    vector<Face> m_faces;
    vector<vec3> m_vertices;

    VertexBufferObject *m_vbo;
    VertexBufferObject *m_vboLines;

};

#endif