DEFINES += _USE_MATH_DEFINES
CONFIG += console
QT += widgets core opengl xml
TARGET = Shape_Grammar

!win32 {
QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXXFLAGS_WARN_ON += -Wall
}

win32 {
QMAKE_CXXFLAGS_RELEASE += /Zi
QMAKE_LFLAGS_RELEASE += /DEBUG
QMAKE_CXXFLAGS_WARN_ON += /W0
QMAKE_CXXFLAGS_WARN_ON += /W0

QMAKE_CXXFLAGS += -MP8
QMAKE_LFLAGS += -LARGEADDRESSAWARE

QMAKE_CFLAGS_RELEASE += /MD
QMAKE_CXXFLAGS_RELEASE += /MD
QMAKE_CFLAGS_DEBUG += /MDd
QMAKE_CXXFLAGS_DEBUG += /MDd

DEFINES += WIN32

INCLUDEPATH += \
D:\Source\ShapeGrammars\libs\glew-1.12.0\include \
D:\Source\ShapeGrammars\libs\freeglut-3.0.0\include \
D:\Source\ShapeGrammars\libs\flann-1.7.1\include \

CONFIG(release, debug | release){
LIBS += \
D:\Source\ShapeGrammars\libs\glew-1.12.0\lib\Release\x64\glew32.lib \
D:\Source\ShapeGrammars\libs\freeglut-3.0.0\lib\Release\freeglut.lib \
D:\Source\ShapeGrammars\libs\flann-1.7.1\lib\flann.lib \
}
else {
LIBS += \
D:\Source\ShapeGrammars\libs\glew-1.12.0\lib\Release\x64\glew32.lib \
D:\Source\ShapeGrammars\libs\freeglut-3.0.0\lib\Debug\freeglutd.lib \
D:\Source\ShapeGrammars\libs\flann-1.7.1\lib\flann.lib \
}

}

linux {
HEADERS +=  \
    Src/large_integer.h \
    HEADERS +=  \
        Src/large_integer.h
    DEFINES += LINUX
    LIBS += -L"/usr/local/MATLAB/R2016a/bin/glnxa64/" -leng -lmx
    LIBS += \
    /usr/local/lib/libglut.so \
    /usr/local/lib/libGLEW.so
    INCLUDEPATH += /usr/local/MATLAB/R2016a/extern/include \
    ./ \
    ./Src \
    ./Shader \
    ./../grammar/
}

RESOURCES = main.qrc \
    main.qrc

HEADERS +=  \
    Src/Camera.h \
    Src/CameraManager.h \
    Src/Color.h \
    Src/FrameBufferObject.h \
    Src/Geometry.h \
    Src/Global.h \
    Src/GLWidget.h \
    Src/GUI.h \
    Src/Headers.h \
    Src/Light.h \
    Src/Math.h \
    Src/Mesh.h \
    Src/NiceGrid.h \
    Src/Object.h \
    Src/ObjLoader.h \
    Src/PerlinNoise.h \
    Src/Renderer.h \
    Src/Scene.h \
    Src/Shader.h \
    Src/Statistics.h \
    Src/Texture.h \
    Src/TransformFeedback.h \
    Src/VertexBufferObject.h \
    Src/Grammar.h \
    Src/Box.h \
    Src/Node.h \
    Src/Connector.h \
    Src/desolver.h \
    Src/dssolver.h \
    Src/Optimizer.h \
    Src/Proc_Model.h \
    Src/Scanner.h \
    Src/CMAESWrapper.h \
    Src/cma-es/cmaes.h \
    Src/cma-es/parameters.h \
    Src/cma-es/random.h \
    Src/cma-es/timings.h \
    Src/cma-es/utils.h \
    Src/cma-es/cmaes.h \
    Src/cma-es/parameters.h \
    Src/cma-es/random.h \
    Src/cma-es/timings.h \
    Src/cma-es/utils.h \
    Src/SimulatedAnnealing.h \
    Src/Genetic_Algorithm.h \
    Src/Swarm.h \
    Src/pointcloudio.h \
    Src/postoffice.h \
    Src/Clustering.h \
    Src/flann-1.7.1/include/flann/defines.h \
    Src/flann-1.7.1/include/flann/config.h \
    Src/flann-1.7.1/include/flann/flann.h \
    Src/flann-1.7.1/include/flann/general.h \
    Src/flann-1.7.1/include/flann/flann.hpp \    
    Src/Occupancy.h 

SOURCES +=  \
    Main.cpp \
    Src/Camera.cpp \
    Src/CameraManager.cpp \
    Src/Color.cpp \
    Src/FrameBufferObject.cpp \
    Src/Geometry.cpp \
    Src/GLWidget.cpp \
    Src/GUI.cpp \
    Src/Headers.cpp \
    Src/Light.cpp \
    Src/Math.cpp \
    Src/Mesh.cpp \
    Src/NiceGrid.cpp \
    Src/Object.cpp \
    Src/ObjLoader.cpp \
    Src/PerlinNoise.cpp \
    Src/Renderer.cpp \
    Src/Scene.cpp \
    Src/Shader.cpp \
    Src/Statistics.cpp \
    Src/Texture.cpp \
    Src/TransformFeedback.cpp \
    Src/VertexBufferObject.cpp \
    Src/Box.cpp \
    Src/Node.cpp \
    Src/Connector.cpp \
    Src/Grammar.cpp \
    Src/desolver.cpp \
    Src/dssolver.cpp \
    Src/Optimizer.cpp \
    Src/Proc_Model.cpp \
    Src/Scanner.cpp \
    Src/CMAESWrapper.cpp \
    Src/SimulatedAnnealing.cpp \
    Src/Swarm.cpp \
    Src/Genetic_Algorithm.cpp \
    Src/pointcloudio.cpp \
    Src/postoffice.cpp \
    Src/Clustering.cpp \
    Src/Occupancy.cpp 