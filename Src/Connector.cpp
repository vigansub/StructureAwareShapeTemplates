#include "Connector.h"
#include "Box.h"
#include "Node.h"

Connector::Connector(Box *box, Node *p, Node *c)
: m_box(box), 
  m_location(TOP), 
  m_to(c),
  m_from(p),
  m_size(1, 0.005, 1),
  m_translation(0, 0, 0),
  m_offset(0, 0, 0),
  m_transformation(mat4::identitiy()), 
  m_thickness(0.001), 
  m_offsetValues(0.0f, 0.0f)
{
   setTranslation();
}

Connector::~Connector()
{
}

void Connector::render(const Transform &trans, int shaderSelector)
{ 
    mat4 m = mat4::translate(vec3(0, 0.0, 0));
    m_transformation = m*mat4::translate(m_translation) * mat4::translate(m_offset) * mat4::scale(m_size);
    m_box->render(trans, m_transformation, 1, shaderSelector, true);
}

vec3 Connector::translation()
{
    return m_translation;
}

vec3 Connector::offset()
{
    return m_offset;
}

vec3 Connector::size()
{
    return m_size;
}

void Connector::setSize(float x, float y)
{
    vec3 sizeParent  = m_from->size();

    switch(m_location)
    {
        case TOP:
        case BOTTOM:
            m_size = vec3(min(x, sizeParent.x), m_thickness, min(y, sizeParent.z));
        break;
        case LEFT:
        case RIGHT:
            m_size = vec3(m_thickness, min(y, sizeParent.y), min(x, sizeParent.z));
        break;
        case FRONT:
        case BACK:
            m_size = vec3(min(x, sizeParent.x), min(y, sizeParent.y), m_thickness);
        break;
    } 
}

void Connector::setTranslation()
{
    vec3 sizeParent  = m_from->size();
    vec3 transParent = m_from->translation();

    switch(m_location)
    {
        case TOP:
            m_translation = transParent + toTopOffset(sizeParent);            
        break;
        case BOTTOM:
            m_translation = transParent + toBottomOffset(sizeParent);
        break;
        case LEFT:
            m_translation = transParent + toLeftOffset(sizeParent);
        break;
        case RIGHT:
            m_translation = transParent + toRightOffset(sizeParent);
        break;
        case FRONT:
            m_translation = transParent + toFrontOffset(sizeParent);
        break;
        case BACK:
            m_translation = transParent + toBackOffset(sizeParent);
        break;
    } 
}

void Connector::setOffset(float s, float t)
{
    t = max(min(1.0f, t), 0.0f);
    s = max(min(1.0f, s), 0.0f);

    m_offsetValues = vec2(t, s);

    vec3 sp = m_from->size();
    vec3 trans = vec3();//m_translation;

    if(m_location == TOP)
    {        
        vec3 mi = vec3(trans.x - sp.x*0.5, 0.0, trans.z - sp.z * 0.5) + vec3(m_size.x*0.5, 0.0, m_size.z*0.5);
        vec3 ma = vec3(trans.x + sp.x*0.5, 0.0, trans.z + sp.z * 0.5) - vec3(m_size.x*0.5, 0.0, m_size.z*0.5);
        m_offset = vec3(lerp(mi.x, ma.x, s), 0.0, lerp(mi.z, ma.z, 1-t));
    }

    if(m_location == BOTTOM)
    {        
        vec3 mi = vec3(trans.x - sp.x*0.5, 0.0, trans.z - sp.z * 0.5) + vec3(m_size.x*0.5, 0.0, m_size.z*0.5);
        vec3 ma = vec3(trans.x + sp.x*0.5, 0.0, trans.z + sp.z * 0.5) - vec3(m_size.x*0.5, 0.0, m_size.z*0.5);
        m_offset = vec3(lerp(mi.x, ma.x, s), 0.0, lerp(mi.z, ma.z, 1-t));
    }

    if(m_location == LEFT)
    {
        vec3 mi = vec3(0.0, trans.y - sp.y*0.5, trans.z - sp.z * 0.5) + vec3(0.0, m_size.y*0.5, m_size.z*0.5);
        vec3 ma = vec3(0.0, trans.y + sp.y*0.5, trans.z + sp.z * 0.5) - vec3(0.0, m_size.y*0.5, m_size.z*0.5);
        m_offset = vec3(0.0, lerp(mi.y, ma.y, t), lerp(mi.z, ma.z, s));
    }

    if(m_location == RIGHT)
    {
        vec3 mi = vec3(0.0, trans.y - sp.y*0.5, trans.z - sp.z * 0.5) + vec3(0.0, m_size.y*0.5, m_size.z*0.5);
        vec3 ma = vec3(0.0, trans.y + sp.y*0.5, trans.z + sp.z * 0.5) - vec3(0.0, m_size.y*0.5, m_size.z*0.5);
        
        m_offset = vec3(0.0, lerp(mi.y, ma.y, t), lerp(mi.z, ma.z, 1-s));
    }

    if(m_location == FRONT)
    {
        vec3 mi = vec3(trans.x - sp.x*0.5, trans.y - sp.y * 0.5, 0.0) + vec3(m_size.x*0.5, m_size.y*0.5, 0.0);
        vec3 ma = vec3(trans.x + sp.x*0.5, trans.y + sp.y * 0.5, 0.0) - vec3(m_size.x*0.5, m_size.y*0.5, 0.0);
        m_offset = vec3(lerp(mi.x, ma.x, s), lerp(mi.y, ma.y, t), 0.0);
    }

    if(m_location == BACK)
    {
        vec3 mi = vec3(trans.x - sp.x*0.5, trans.y - sp.y * 0.5, 0.0) + vec3(m_size.x*0.5, m_size.y*0.5, 0.0);
        vec3 ma = vec3(trans.x + sp.x*0.5, trans.y + sp.y * 0.5, 0.0) - vec3(m_size.x*0.5, m_size.y*0.5, 0.0);
        m_offset = vec3(lerp(mi.x, ma.x, 1-s), lerp(mi.y, ma.y, t), 0.0);
    }

}

void Connector::adjustOffset(float x, float y)
{
    m_offsetValues.x = max(min(1.0f, x), 0.0f);
    m_offsetValues.y = max(min(1.0f, y), 0.0f);

    setOffset(m_offsetValues.x, m_offsetValues.y);
}

Connector::Location Connector::location()
{
    return m_location;
}

Node *Connector::from()
{
    return m_from;
}

Node* Connector::to()
{
    return m_to;
}

void Connector::setLocation(Location location)
{
    m_location = location;
    setTranslation();
}

void Connector::setOffsetRnd(const vec3 &mi, const vec3 &ma)
{
    this->setOffset(rand(mi.x, ma.x), rand(mi.y, ma.y));
}

void Connector::setSizeRnd(const vec3 &mi, const vec3 &ma)
{
    this->setSize(rand(mi.x, ma.x), rand(mi.y, ma.y));
}

bool Connector::intersects(Connector *c)
{
    if(this->location() != c->location())
        return false;

    vec3 ao = this->offset();
    vec3 as = this->size();
                
    vec3 co = c->offset();
    vec3 cs = c->size();

    float tw, th, rw, rh, tx, ty, rx, ry;

    if(this->location() == TOP || this->location() == BOTTOM)
    {
        tw = as.x;
        th = as.z;
        rw = cs.x;
        rh = cs.z;

        tx = ao.x;
        ty = ao.z;
        rx = co.x;
        ry = co.z;
    }

    if(this->location() == FRONT || this->location() == BACK)
    {
        tw = as.x;
        th = as.y;
        rw = cs.x;
        rh = cs.y;

        tx = ao.x;
        ty = ao.y;
        rx = co.x;
        ry = co.y;
    }

    if(this->location() == LEFT || this->location() == RIGHT)
    {
        tw = as.z;
        th = as.y;
        rw = cs.z;
        rh = cs.y;

        tx = ao.z;
        ty = ao.y;
        rx = co.z;
        ry = co.y;
    }

    if (rw <= 0 || rh <= 0 || tw <= 0 || th <= 0) 
    {
        return false;
    }

    rw += rx;
    rh += ry;
    tw += tx;
    th += ty;

    //      overflow || intersect
    return ((rw < rx || rw > tx) &&
            (rh < ry || rh > ty) &&
            (tw < tx || tw > rx) &&
            (th < ty || th > ry));

}