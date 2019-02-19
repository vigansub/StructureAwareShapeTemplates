#ifndef DISTANCEFIELD
#define DISTANCEFIELD

#include "Headers.h"

class VertexBufferObjectAttribs;
class Shader;
class Texture;

class DistanceField
{
public:
    DistanceField(int w, int h, int d);
    ~DistanceField();   

	struct Face
	{
		uint a, b, c, d;
	};


    void computeDistance(const vector<vec3> &points);
    float dist(int x, int y, int z);

    void buildVBO();
    void buildVBOLines();
    void buildVBOPlanes1();
    void buildVBOPlanes2();
    void render(const Transform &trans);
    void writeObj(vector<vec3> &positions, vector<vec3> &normals);

    void buildTextures();

    float width();
    float height();
    float depth();

    void save(const QString &path);
    void load(const QString &path);

    vector<float> distances();

    float computeMax();
    float max();

    int nr();
    void print();

private:
    float *m_data;
    int m_w, m_h, m_d;
    float m_maxDist;

    VertexBufferObjectAttribs *m_vbo;
    VertexBufferObjectAttribs *m_vboLines;
    VertexBufferObjectAttribs *m_vboPlanes1;
    VertexBufferObjectAttribs *m_vboPlanes2;

    Texture *m_tex1;
    Texture *m_tex2;

    float m_gridOffset;

};

#endif

