#include "Node.h"
#include "Box.h"
#include "Connector.h"

Node::Node(Box *box, const QString &name) 
: m_box(box), 
  m_size(vec3(1, 1, 1)), 
  m_translation(vec3(0, -0.5, 0)),
  m_transformation(mat4::identitiy()), 
  m_name(name), 
  m_constrainSize(0),
  m_avoidOverlapping(true), 
  m_colorMode(0)
{    
}

Node::~Node()
{
}

void Node::render(const Transform &trans, int shaderSelector)
{
    mat4 m = mat4::translate(vec3(0, 0, 0));
    m_transformation = m * mat4::translate(m_translation) * mat4::scale(m_size);
    m_box->render(trans, m_transformation, m_colorMode, shaderSelector, false);
}

void Node::renderDepth(const Transform &trans, int shaderSelector)
{    
    m_transformation = mat4::translate(m_translation) * mat4::scale(m_size);
    m_box->renderDepth(trans, m_transformation, m_colorMode, shaderSelector);
}

void Node::setSize(const vec3 &size)
{
    if(m_parents.size() == 0)
    {
        m_size = size;
    }
    else
    {
        Connector *c = m_parents[0];
        vec3 s = c->size();    

        switch(c->location())
        {
            case Connector::TOP:            
            case Connector::BOTTOM:
                m_size = vec3(s.x, size.y, s.z);
            break;
            case Connector::LEFT:
            case Connector::RIGHT:
                m_size = vec3(size.x, s.y, s.z);
            break;
            case Connector::FRONT:
            case Connector::BACK:
                m_size = vec3(s.x, s.y, size.z);
            break;
        }

    }
}

void Node::setSizeRnd(const vec3 &mi, const vec3 &ma)
{
    this->setSize(vec3(rand(mi.x, ma.x), rand(mi.y, ma.y), rand(mi.z, ma.z)));
}

void Node::setTranslation(const vec3 &trans)
{
    m_translation = trans;    
}

void Node::addParent(Connector *c)
{
    if(m_parents.size() == 0)
    {
        m_parents.push_back(c);

        Connector::Location location = c->location();
        
        vec3 s = c->size();
        vec3 t = c->translation();
        vec3 o = c->offset();

        switch(location)
        {
            case Connector::TOP:     
                if(m_constrainSize == 0)
                    m_size = vec3(s.x, m_size.y, s.z);
                else if(m_constrainSize == 1)
                    m_size = vec3(s.x, m_size.y, m_size.z);
                else if(m_constrainSize == 2)
                    m_size = vec3(m_size.x, m_size.y, s.z);
                else if(m_constrainSize == 3)
                    m_size = vec3(m_size.x, m_size.y, m_size.z);

                m_translation = t + o + toTopOffset(m_size);
            break;
            case Connector::BOTTOM:
                if(m_constrainSize == 0)
                    m_size = vec3(s.x, m_size.y, s.z);
                else if(m_constrainSize == 1)
                    m_size = vec3(s.x, m_size.y, m_size.z);
                else if(m_constrainSize == 2)
                    m_size = vec3(m_size.x, m_size.y, s.z);
                else if(m_constrainSize == 3)
                    m_size = vec3(m_size.x, m_size.y, m_size.z);

                m_translation = t + o + toBottomOffset(m_size);
            break;
            case Connector::LEFT:
                if(m_constrainSize == 0)
                    m_size = vec3(m_size.x, s.y, s.z);
                else if(m_constrainSize == 1)
                    m_size = vec3(m_size.x, m_size.y, s.z);
                else if(m_constrainSize == 2)
                    m_size = vec3(m_size.x, s.y, m_size.z);
                else if(m_constrainSize == 3)
                    m_size = vec3(m_size.x, m_size.y, m_size.z);

                m_translation = t + o + toLeftOffset(m_size);
            break;
            case Connector::RIGHT:
                if(m_constrainSize == 0)
                    m_size = vec3(m_size.x, s.y, s.z);
                else if(m_constrainSize == 1)
                    m_size = vec3(m_size.x, m_size.y, s.z);
                else if(m_constrainSize == 2)
                    m_size = vec3(m_size.x, s.y, m_size.z);
                else if(m_constrainSize == 3)
                    m_size = vec3(m_size.x, m_size.y, m_size.z);

                m_translation = t + o + toRightOffset(m_size);
            break;
            case Connector::FRONT:
                if(m_constrainSize == 0)
                    m_size = vec3(s.x, s.y, m_size.z);
                else if(m_constrainSize == 1)
                    m_size = vec3(s.x, m_size.y, m_size.z);
                else if(m_constrainSize == 2)
                    m_size = vec3(m_size.x, s.y, m_size.z);
                else if(m_constrainSize == 3)
                    m_size = vec3(m_size.x, m_size.y, m_size.z);

                m_translation = t + o + toFrontOffset(m_size);
            break;
            case Connector::BACK:
                if(m_constrainSize == 0)
                    m_size = vec3(s.x, s.y, m_size.z);
                else if(m_constrainSize == 1)
                    m_size = vec3(s.x, m_size.y, m_size.z);
                else if(m_constrainSize == 2)
                    m_size = vec3(m_size.x, s.y, m_size.z);
                else if(m_constrainSize == 3)
                    m_size = vec3(m_size.x, m_size.y, m_size.z);

                m_translation =  t + o + toBackOffset(m_size);
            break;
        }

        c->setTranslation();
    }
    else if(m_parents.size() == 1)
    {
        m_parents.push_back(c);

        Connector *p1 = m_parents[0];
        Connector *p2 = m_parents[1];

        Connector::Location t1 = p1->location();
        Connector::Location t2 = p2->location();

        vec3 size1 = p1->size();
        vec3 translation1 = p1->translation();
        vec3 offset1 = p1->offset();

        vec3 size2 = p2->size();
        vec3 translation2 = p2->translation();
        vec3 offset2 = p2->offset();

        vec3 r1 = translation1 + offset1;
        vec3 r2 = translation2 + offset2;

        vec3 d = normalize(r1-r2);

        int dot1 = abs(dot(d, vec3(1, 0, 0)));
        int dot2 = abs(dot(d, vec3(0, 1, 0)));
        int dot3 = abs(dot(d, vec3(0, 0, 1)));

        if(dot1 == 1 || dot2 == 1 || dot3 == 1)
        {
            if(t1 == Connector::TOP && t2 == Connector::BOTTOM)
            {
                m_size.y = length(r1-r2);
                m_translation = translation1 + offset1 + toTopOffset(m_size);               
            }
            if(t1 == Connector::BOTTOM && t2 == Connector::TOP)
            {
                m_size.y = length(r1-r2);
                m_translation = translation2 + offset2 + toTopOffset(m_size);
            }

            if(t1 == Connector::LEFT && t2 == Connector::RIGHT)
            {
                m_size.x = length(r1-r2);
                m_translation = translation1 + offset1 + toLeftOffset(m_size);
            }

            if(t1 == Connector::RIGHT && t2 == Connector::LEFT)
            {
                m_size.x = length(r1-r2);
                m_translation = translation1 + offset1 + toRightOffset(m_size);
            }

            if(t1 == Connector::FRONT && t2 == Connector::BACK)
            {
                m_size.z = length(r1-r2);
                m_translation = translation1 + offset1 + toFrontOffset(m_size);
            }

            if(t1 == Connector::BACK && t2 == Connector::FRONT)
            {
                m_size.z = length(r1-r2);
                m_translation = translation1 + offset1 + toBackOffset(m_size);
            }
        }
    }
}

bool Node::addChild(Connector *c)
{
    //if(m_avoidOverlapping)
    //{
    //    bool add = true;
    //    for(int i=0; i<m_children.size(); ++i)
    //    {
    //        Connector *a = m_children[i];

    //        if(a->intersects(c))
    //        {
    //            add = false;
    //            return false;
    //        }
    //    }

    //    if(add == true)
    //    {
    //        m_children.push_back(c);
    //        return true;
    //    }        
    //}
    //else
    //{
        m_children.push_back(c);
        return true;
    //}

}

vec3 Node::size()
{
    return m_size;
}

vec3 Node::translation()
{
    return m_translation;
}

QString Node::name()
{
    return m_name;
}

QString Node::id()
{
    Node *n = this;
//    int address = (int)this;
//    return address;
//    return QString((char*)this);
    return QString();
}

void Node::setColorMode(int mode)
{
    m_colorMode = mode;
}

void Node::setConstrainSize(int constrain)
{
    m_constrainSize = constrain;
}

mat4 Node::tranform()
{
    return mat4::translate(m_translation) * mat4::scale(m_size);
}
