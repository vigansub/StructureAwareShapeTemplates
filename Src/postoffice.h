#ifndef POSTOFFICE_H
#define POSTOFFICE_H
class GLWidget;
class Renderer;

namespace Postoffice
{
GLWidget *glwidget();
void setGLWidget(GLWidget *_w);

Renderer *renderer();
void setRenderer(Renderer *_w);
}

#endif // POSTOFFICE_H
