#include <QApplication>
#include "Src/Headers.h"
#include "Src/GLWidget.h"

//#include "matlabwrapper.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    glutInit(&argc, argv);

    QGLFormat format;
    format.setSampleBuffers(true);
    format.setDoubleBuffer(true);
    format.setRgba(true);
    format.setDirectRendering(true);
    format.setSamples(8);

    QGLContext *context = new QGLContext(format);

    GLWidget widget(context, 2000, 1200);
    widget.setWindowTitle("");
    widget.setGeometry(100, 100, 2000, 1200);
    widget.show();

//    matlab::example_2();
//    matlab::example_3();

    return app.exec();
}
