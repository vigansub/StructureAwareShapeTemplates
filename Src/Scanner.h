#ifndef SCANNER_H
#define SCANNER_H

#include "Headers.h"
#include "VertexBufferObject.h"

#include <qhash.h>

class Shader;
class VertexBufferObjectAttribs;
class Scene;
class FrameBufferObject;

class Scanner
{
public:
    Scanner(Scene *scene);
    ~Scanner();

    void renderScan(const Transform &trans, bool useView);
    void storeView();
    void resize(int w, int h);

    void buildVBO();
    void updateBuffer();
    void renderPointCloud(const Transform &trans);
    void renderPointCloudDepth(const Transform &trans);

    int setToNextPos();
    void clearBuffer();

    void scan();
    void autoScan();

    FrameBufferObject *m_fbo;
    bool m_scanning;

    // TODO or TOREMOVE
    vector<vec3> scan_data();

private:
    Scene *m_scene;
    int m_bufferWidth;
    int m_bufferHeight;    

    VertexBufferObject *m_vbo;
    QHash<QString, Vertex> m_vertices;

    int m_nrBufferVertices;
    bool m_moveAutomatic;

    float m_fcp;
    float m_ncp;
    float m_aspect;
    float m_fov;
    float m_rotHeight;
    float m_zoom;

    vec2 m_activePos;    
};

#endif
