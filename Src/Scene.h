#ifndef SCENE_H
#define SCENE_H

#include "Headers.h"
#include "Geometry.h"
#include "CMAESWrapper.h"
#include "Optimizer.h"
#include "Proc_Model.h"
#include "SimulatedAnnealing.h"
#include "Swarm.h"
#include "Genetic_Algorithm.h"
#include "Clustering.h"
#include "flann-1.7.1\include\flann\flann.h"
#include "Occupancy.h"

class NiceGrid;
class Light;
class VertexBufferObject;
class Shader;
class CameraManager;
class Object;
class TransformFeedback;
class Grammar;
class GLWidget;

class Scene : public QObject
{
    Q_OBJECT

public:
    Scene(GLWidget *parent, CameraManager *camManager);
    ~Scene();

	void update(float delta);
    void init();

    void renderWorld(const Transform &trans);  
    void renderObjects(const Transform &trans);  
	void renderObjectsScan(const Transform &trans);
    void renderObjectsDepth(const Transform &trans);

    void select(const Transform &trans, int sw, int sh, int mx, int my);
    void move(const Transform &trans, int x, int y);
    void resetSelection();

    void buildVBO();
	void optimize();
	void buildVBOSolution(const vector<vec3> &positions, const vector<vec3> &sizes);
	void buildVBOTarget(const vector<vec3> &positions, const vector<vec3> &sizes);
    void buildVBOVoxels();
	
	void initializeProceduralModel();

    void buildBox(vector<vec3> &vertices, vector<vec3> &normals, vector<vec3> &verticesLines, vec3 p, float s);

	void testIndependentVariables();
    void test_Genetic();
    void test_Swarm();
    void test_Simulated_Anneal();
	void test_DESolverProcModel();
	void test_DESolverGrammar();
    void test_CMAESGrammar();
    void test_MultipleShapes();
    vector<double> test_CMAESGrammarParallel(int i, vector<vec3> samplevals);
    void test_ParallelRun(vector<vec3> points);

    void generateErrorMatrix();
    void generateErrorMatrix_Parallel();
    void generateNewData();
    void saveAllImages();

public:
    vector<Light *> m_lights;
    Grammar *m_grammar;
    vector<Grammar *> m_grammarvec;
	QHash<QString, Vertex> m_verticesScan;
    const int m_num_parallel = 8;
    vector<vector<int>> m_neighbors;

    vector<vector<double>> m_optimal;
    vector<double> m_mean, m_std;
    vector<vector<double>> m_errormatrix;

    void newGenerate();
    void skinGenerate();
    void renderOptimal();
    void generateNewDatabase();
    void renderOptimalTrue();
    void newShapeResize();
    void writeBoxes_grammars();
    void writeBoxes();
    void obtainBestFitShape();
    void mergeShapes(QString shape1, QString shape2);
    void obtain_gram(); 
    void saveScene(QString path);
    void buildVBOModelDifference(vector<Vertex> &vertices);


    void erasePoints(float dx, float dy, float dh);
    float m_ballEraserRadius;
    vec3 m_ballEraserPos;
    
    vector<float> computeMeshDistance();
    void generateContinuousSkinning();
    void generateContinuousHeatMap();
    void compareMinhyuk();

    void saveObjboxes(QString filename);
    void savePartialResults();
    void sceneCompletion();
    vec3 minimum(vec3 x, vec3 y);
    vec3 maximum(vec3 x, vec3 y);
    double Scene::norm(vector<double> x, vector<double> y);
    void findClosestShape();

private:
    void load_object();
    void load_object(QString path);
    void load_object(QString path, bool normalize);
    
    vector<vector<double>> load_matrix(string filename);

    NiceGrid *m_niceGrid;
    CameraManager *m_cameraManager;
    
    VertexBufferObject *m_vbo;
	VertexBufferObject *m_vboSolution;
	VertexBufferObject *m_vboTarget;
    VertexBufferObject *m_vboVoxels;
    VertexBufferObject *m_vboVoxelsLines;
    vector<Object*> m_objects;

    int m_activeIdx;

	Optimizer *m_optimizer;
    CMAESWrapper *m_cmaes;
    vector<CMAESWrapper *> m_cmaesvec;
    Sim_Anneal *m_sim_anneal;
    Swarm *m_swarm;
    Genetic* m_genetic;

    Clusters *m_clusters;
    
	Proc_Model *m_pmodel;

    Object *m_object = nullptr;
    GLWidget *m_parent;

    double m_error;
    float m_shapevolume;
    vector<vector<vector<int>>> m_grid;
    vector<vec3> m_gridbbox;

    QString m_folder;
    QString m_objectname;

    vector<vector<int>> m_voxels;

    VertexBufferObject *m_vboModelDifference;
    VertexBufferObject *m_vboBallEraser;
    
    

    int m_activeLayer;
public slots:
    void updateGrammarGeometry();
    void updateGrammarGeometry_CMAES();
    void updateGrammarGeometry_SimAnneal();
    void updateGrammarGeometry_Swarm();
};

#endif


