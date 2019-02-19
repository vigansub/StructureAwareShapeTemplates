#ifndef CONNECTOR_H
#define CONNECTOR_H

#include "Headers.h"

class Box;
class Node;

class Connector
{
   
public:
    enum Location 
    {
        FRONT = 0, 
        BACK, 
        LEFT, 
        RIGHT, 
        TOP, 
        BOTTOM
    };

    Connector(Box *box, Node *p, Node *c);
    ~Connector();
	
   void render(const Transform &trans, int shaderSelector = 0);
   void setSize(float x, float y);
   void setTranslation();
   void setOffset(float s, float t);
   void setOffsetRnd(const vec3 &mi, const vec3 &ma);
   void setSizeRnd(const vec3 &mi, const vec3 &ma);
   void setLocation(Location location);
   bool intersects(Connector *c);
   void adjustOffset(float x, float y);

   vec3 translation();
   vec3 offset();

   vec3 size();
   
   Location location();

   Node* from();
   Node* to();

private:
    vec3 m_size;    
    vec3 m_translation;
    vec3 m_offset;
    mat4 m_transformation;
    vec2 m_offsetValues;

    Node *m_from;
    Node *m_to;

    Box *m_box;

    float m_thickness;

    Location m_location;
};

#endif