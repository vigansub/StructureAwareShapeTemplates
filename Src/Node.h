#ifndef NODe_H
#define NODE_H

#include "Headers.h"
#include "Connector.h"

class Node
{
public:
    Node(Box *box, const QString &name);
    ~Node();

    void addParent(Connector *c);
    bool addChild(Connector *c);

    void setSize(const vec3 &size);
    void setSizeRnd(const vec3 &mi, const vec3 &ma);
    void setTranslation(const vec3 &trans);

    vec3 size();
    vec3 translation();
    mat4 tranform();

    void render(const Transform &trans, int shaderSelector = 0);
    void renderDepth(const Transform &trans, int shaderSelector = 0);

    QString name();
    QString id();

    void setColorMode(int mode);
    void setConstrainSize(int constrain);
    int colorMode() { return m_colorMode; }

private:
    vector<Connector *> m_parents;
    vector<Connector *> m_children;

    vec3 m_size;    
    vec3 m_translation;

    mat4 m_transformation;

    Box *m_box;

    QString m_name;
    int m_constrainSize;
    bool m_avoidOverlapping;

    int m_colorMode;

};

#endif