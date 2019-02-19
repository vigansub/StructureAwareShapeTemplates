#include "postoffice.h"

GLWidget *widget_ = nullptr;
Renderer *renderer_ = nullptr;


GLWidget *Postoffice::glwidget()
{
    return widget_;
}

void Postoffice::setGLWidget(GLWidget *_w)
{
    widget_ = _w;
}

Renderer *Postoffice::renderer()
{
    return renderer_;
}

void Postoffice::setRenderer(Renderer *_w)
{
    renderer_ = _w;
}
