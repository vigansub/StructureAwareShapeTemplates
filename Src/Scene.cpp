#include "Scene.h"
#include "NiceGrid.h"
#include "Light.h"
#include "Shader.h"
#include "VertexBufferObject.h"
#include "Mesh.h"
#include "CameraManager.h"
#include "Object.h"
#include "Geometry.h"
#include "TransformFeedback.h"
#include "Grammar.h"
#include "GLWidget.h"
#include "pointcloudio.h"
#include "postoffice.h"
#include "Renderer.h"
#include "Scanner.h"
#include "QDir"
#include "QDirIterator"
#include <ppl.h>
#include <fstream>
#include "Node.h"
#include <vector>
#include <thread>
#include <omp.h>

Scene::Scene(GLWidget *parent, CameraManager *camManager)
    : m_parent(parent),
    m_cameraManager(camManager),
    m_activeIdx(-1),
    m_vboSolution(nullptr),
    m_vboTarget(nullptr),
    m_optimizer(nullptr),
    m_cmaes(nullptr),
    m_sim_anneal(nullptr),
    m_swarm(nullptr),
    m_genetic(nullptr),
    m_clusters(nullptr),
    m_vboVoxels(nullptr),
    m_vboVoxelsLines(nullptr),
    m_activeLayer(32),
    m_ballEraserPos(vec3(0, 2, 0)),
    m_vboBallEraser(nullptr), 
    m_ballEraserRadius(0.1f), 
    m_vboModelDifference(nullptr)
{
    init();   

    buildVBOVoxels();
}

Scene::~Scene()
{
}

void Scene::init()
{
    m_lights.push_back(new Light(this, Light::SPOT_LIGHT, vec3(0.9f), vec3(-4.3, 7.7, 7.2),  vec3(), vec3(1.2f), vec3(), vec3(0.7f, 0.001f, 0.0001f)));  

	m_niceGrid = new NiceGrid(100.0f, 40.0f);  

    m_grammar = new Grammar();

    m_grammar->softCleanUp();
    m_grammar->derive();
    int p = m_grammar->determineNrIndependentVariables();

    m_parent->updateGL();

    m_error = 1000.0;

    m_grammarvec.resize(m_num_parallel);
    for (int i = 0; i < m_num_parallel; ++i)
        m_grammarvec[i] = new Grammar();

    load_object();

    m_vboBallEraser = Mesh::sphere(m_ballEraserRadius, 4, vec4(1, 0, 0, 1));
}   

void Scene::test_MultipleShapes()
{ // locate the folder
#ifdef WIN32
    m_folder = "Data/Objs/SpecChairs/simple/test/";
#else
    m_folder = "./../_main/Data/Objs/";
#endif
    m_optimal.clear();
    m_mean.clear();
    m_std.clear();
    srand(time(0));

    QDirIterator it(m_folder, QStringList() << "*.obj");
    while (it.hasNext())
    {
        // scan all the obj file
        for (int i = 0; i < m_num_parallel; ++i)
            m_grammarvec[i] = new Grammar();


        qDebug() << it.next();
        m_objectname = it.fileName().left(it.fileName().lastIndexOf("."));
        qDebug() << m_objectname;
        m_verticesScan.clear();
        load_object(m_folder + m_objectname + ".obj");
        QString cloud_name = m_folder + m_objectname + ".cld";
        vector<vec3> points = PointCloudIO::load_point_cloud(cloud_name);
        // check if the corresponding point cloud files exist
        // if the point clouds do not exist, scan and save them
        if (points.size() <= 10)
        {
            points.clear();
            for (int i = 0; i < 10; ++i)
            {
                Postoffice::renderer()->m_scanner->m_scanning = true;
                Postoffice::renderer()->m_scanner->autoScan();
            }
            // save
            QList<Vertex> pList = m_verticesScan.values();
            points.reserve(pList.size());
            std::transform(pList.begin(), pList.end(), back_inserter(points),
                [](Vertex V){return V.position; });
            qDebug() << "Saving Point Cloud";
            PointCloudIO::save_point_cloud(points, cloud_name);
        }
        //qDebug() << it.fileName();

        if (points.size() >= 10)
        {
            qDebug() << "No. of points: " << points.size();
            //Occupancy k(points);
            float volume;
            fstream fin;
            fin.open((m_folder + m_objectname + "__Volume.txt").toStdString());
            fin >> volume;
            fin.close();
            m_shapevolume = volume;
            test_ParallelRun(points);
        }
    }

    m_mean.resize(m_optimal[0].size());
    m_std.resize(m_mean.size());
    for (int i = 0; i < m_optimal.size(); ++i)
    {
        vector<double> optimal_i = m_optimal[i];
        for (int j = 0; j < optimal_i.size(); ++j)
        {
            m_mean[j] += ((optimal_i[j]) / m_optimal.size());
            std::cout << optimal_i[j] << " ";
        }
        std::cout << endl;
    }
    std::cout << "Mean: " << " ";
    for (int i = 0; i < m_mean.size(); ++i)
        std::cout << m_mean[i] << " ";

    for (int i = 0; i < m_optimal.size(); ++i)
    {
        vector<double> optimal_i = m_optimal[i];
        for (int j = 0; j < m_mean.size(); ++j)
        {
            m_std[j] += (((optimal_i[j] - m_mean[j])*(optimal_i[j] - m_mean[j])) / m_optimal.size());
        }
    }

    std::cout << "Std. Deviation: " << " ";
    for (int i = 0; i < m_std.size(); ++i)
    {
        std::cout << m_std[i] << " ";
        //m_std[i] = sqrt(m_std[i]);
        m_std[i] = 0.01;
    }
}

void Scene::generateErrorMatrix()
{
#ifdef WIN32
    m_folder = "C:/Users/Vignesh/Documents/shape-grammars/Data/Grammars/dresser/";
#else
    m_folder = "./../_main/Data/Objs/";
#endif
    m_optimal.clear();
    m_mean.clear();
    m_std.clear();
    srand(time(0));

    QDirIterator it(m_folder, QStringList() << "*.obj");

    ofstream fout, gout;
    fout.open((m_folder + "ErrorMatrix.txt").toStdString());
    gout.open((m_folder + "OptimalVectors.txt").toStdString());

    while (it.hasNext())
    {
        qDebug() << it.next();
        m_objectname = it.fileName().left(it.fileName().lastIndexOf("."));
        qDebug() << m_objectname; 
        m_verticesScan.clear();
        load_object(m_folder + m_objectname + ".obj");
        QString cloud_name = m_folder + m_objectname + ".cld";
        vector<vec3> points = PointCloudIO::load_point_cloud(cloud_name);
        // check if the corresponding point cloud files exist
        // if the point clouds do not exist, scan and save them
        if (points.size() <= 10)
        {
            points.clear();
            for (int i = 0; i < 10; ++i)
            {
                Postoffice::renderer()->m_scanner->m_scanning = true;
                Postoffice::renderer()->m_scanner->autoScan();
            }
            // save
            QList<Vertex> pList = m_verticesScan.values();
            points.reserve(pList.size());
            std::transform(pList.begin(), pList.end(), back_inserter(points),
                [](Vertex V){return V.position; });
            qDebug() << "Saving Point Cloud";
            PointCloudIO::save_point_cloud(points, cloud_name);
        }

        vector<double> error;
        QDirIterator gt(m_folder, QStringList() << "*.gram");
        int g_index = 0;
        // scan all the obj file
        while (gt.hasNext())
        {
            g_index++;
            qDebug() << gt.next();
            qDebug() << gt.fileName();
            m_grammar = new Grammar(m_folder + gt.fileName());
            for (int i = 0; i < m_num_parallel; ++i)
                m_grammarvec[i] = new Grammar(m_folder + gt.fileName());

            if (points.size() >= 10)
            {
                qDebug() << "Number of points: " << points.size();
                qDebug() << m_folder + "Output/Shape_" + it.fileName() + "_Grammar" + g_index;

                float volume;
                fstream fin;
                fin.open((m_folder + m_objectname + "__Volume.txt").toStdString());
                fin >> volume;
                fin.close();
                m_shapevolume = volume;

                test_ParallelRun(points);
                error.push_back(m_error);
                fout << m_error << endl;

                vector<double> optimal = m_optimal[m_optimal.size() - 1];
                if (optimal.size() >= 1)
                    gout << optimal[0];
                for (int j = 1; j < optimal.size(); ++j)
                    gout << " " << optimal[j];

                gout << endl;
                                
                qDebug() << "Best Error is: " << m_error;
                saveFrameBuffer(m_parent, "Shape_" + it.fileName() + "_Grammar" + gt.fileName());
            }
        }
        fout << endl;
        m_errormatrix.push_back(error);
    }

    for (int i = 0; i < m_errormatrix.size(); ++i)
    {
        vector<double> error = m_errormatrix[i];
        for (int j = 0; j < error.size(); ++j)
            qDebug() << error[j];
    }

    fout.close();
    gout.close();
}

void Scene::generateErrorMatrix_Parallel()
{
#ifdef WIN32
    m_folder = "C:/Users/Vignesh/Documents/shape-grammars/Data/Grammars/sample_folder/";
#else
    m_folder = "./../_main/Data/Objs/";
#endif
    m_optimal.clear();
    m_mean.clear();
    m_std.clear();
    srand(time(0));

    QDirIterator it(m_folder, QStringList() << "*.obj", QDir::Files, QDirIterator::Subdirectories);
    vector<QString> obj_files, gram_files;
    vector<vector<vec3>> all_points;

    while (it.hasNext())
    {
        qDebug() << it.next();
        m_objectname = it.filePath().left(it.filePath().lastIndexOf("."));
        obj_files.push_back(m_objectname);
        qDebug() << m_objectname;
        m_verticesScan.clear();
        QString cloud_name = m_objectname + ".cld";
        vector<vec3> points = PointCloudIO::load_point_cloud(cloud_name);
        // check if the corresponding point cloud files exist
        // if the point clouds do not exist, scan and save them
        if (points.size() <= 10)
        {
            points.clear();
            for (int i = 0; i < 10; ++i)
            {
                Postoffice::renderer()->m_scanner->m_scanning = true;
                Postoffice::renderer()->m_scanner->autoScan();
            }
            // save
            QList<Vertex> pList = m_verticesScan.values();
            points.reserve(pList.size());
            std::transform(pList.begin(), pList.end(), back_inserter(points),
                [](Vertex V){return V.position; });
            qDebug() << "Saving Point Cloud";
            PointCloudIO::save_point_cloud(points, cloud_name);
            //Remove later
            //m_parent->updateGL();
            //m_parent->updateGL();
            //m_parent->updateGL();
            //m_parent->updateGL();
        }

        all_points.push_back(points);
    }

    int num_parallel = obj_files.size();

    /*load_object(m_folder + obj_files[0] + ".obj");
    m_parent->updateGL();
    m_parent->updateGL();
    m_parent->updateGL();
    m_parent->updateGL();*/

    qDebug() << (m_folder + obj_files[0] + ".obj");

    QDirIterator gt(m_folder, QStringList() << "*.gram");
    int g_index = 0;
    // scan all the obj file
    while (gt.hasNext())
    {
        g_index++;
        qDebug() << gt.next();
        qDebug() << gt.fileName();
        gram_files.push_back(gt.fileName());
    }

    m_grammarvec.resize(obj_files.size());

    m_optimal.resize(obj_files.size()*gram_files.size());
    m_errormatrix.resize(obj_files.size());
    for (int i = 0; i < m_errormatrix.size(); ++i)
        m_errormatrix[i].resize(gram_files.size());

    ofstream fout, gout;
    fout.open((m_folder + "ErrorMatrix.txt").toStdString());
    gout.open((m_folder + "OptimalVectors.txt").toStdString());

    m_cmaesvec.resize(num_parallel);
    //memset(m_cmaesvec, nullptr, m_num_parallel);
    for (int i = 0; i < num_parallel; ++i)
        m_cmaesvec[i] = nullptr;


/*#pragma omp parallel for schedule(dynamic) reduction
    {
        for (int i = 0; i < num_parallel; ++i)
        {
            QString cloud_name = m_folder + obj_files[i] + ".cld";
            qDebug() << cloud_name;
            vector<vec3> points = PointCloudIO::load_point_cloud(cloud_name);
            // i++; //Rewrite

            for (int gram = 0; gram < gram_files.size(); ++gram)
            {
                qDebug() << i << " " << gram;
                m_grammarvec[i] = new Grammar(m_folder + gram_files[gram]);
                qDebug() << (m_folder + gram_files[gram]);
                vector<double> optimal_vec = test_CMAESGrammarParallel(i, points);
                m_optimal[i*gram_files.size() + gram] = optimal_vec;
                m_errormatrix[i][gram] = m_cmaesvec[i]->m_besterror;

                fout << (m_folder + obj_files[i]).toStdString() << " " << gram_files[gram].toStdString() << endl;
                fout << m_errormatrix[i][gram] << endl;

                gout << (m_folder + obj_files[i]).toStdString() << " " << gram_files[gram].toStdString() << endl;
                for (int k = 0; k < optimal_vec.size(); ++k)
                    gout << optimal_vec[k] << " ";

                gout << endl;

                // Remove following lines when you don't want to see output.
                //m_grammar = m_grammarvec[i];
                //qDebug() << m_grammar->m_nrIndependentVariables;
                //const double *solution = optimal_vec.data();
                //m_grammar->softCleanUp();
                //m_grammar->adjustParameters(solution, m_grammar->m_nrIndependentVariables);
                //m_grammar->derive();
                //m_parent->updateGL();
                //m_parent->updateGL();
                //m_parent->updateGL();
                //m_parent->updateGL();
                //m_parent->updateGL();
                //m_parent->updateGL();
                //m_parent->updateGL();

            }
        }
    }*/

    const int num_threads = std::thread::hardware_concurrency();
    qDebug() << num_threads;

    std::vector<std::thread> threads(num_threads);
    for (int t = 0; t<num_threads; t++)
    {
        threads[t] = std::thread(std::bind(
            [&](const int bi, const int ei, const int t)
        {
            // loop over all items
            for (int i = bi; i<ei; i++)
            {
                // inner loop
                {
                    QString cloud_name = obj_files[i] + ".cld";
                    qDebug() << cloud_name;
                    vector<vec3> points = PointCloudIO::load_point_cloud(cloud_name);
                    // i++; //Rewrite

                    for (int gram = 0; gram < gram_files.size(); ++gram)
                    {
                        QRegExp separator("[/.()lp]");
                        QStringList a = obj_files[i].split(separator);
                        QStringList b = gram_files[gram].split(separator);
                        int a_intval, b_intval;
                        for (int a_iter = 0; a_iter < a.size(); ++a_iter)
                            if (a[a_iter].toInt() != 0)
                                a_intval = a[a_iter].toInt();

                        for (int b_iter = 0; b_iter < b.size(); ++b_iter)
                            if (b[b_iter].toInt() != 0)
                                b_intval = b[b_iter].toInt();

                        qDebug() << i << " " << gram;
                        m_grammarvec[i] = new Grammar(m_folder + gram_files[gram]);
                        qDebug() << (m_folder + gram_files[gram]);
                        vector<double> optimal_vec = test_CMAESGrammarParallel(i, points);
                        m_optimal[i*gram_files.size() + gram] = optimal_vec;
                        m_errormatrix[i][gram] = m_cmaesvec[i]->m_besterror;

                        fout << (m_folder + obj_files[i]).toStdString() << " " << gram_files[gram].toStdString() << endl;
                        fout << m_errormatrix[i][gram] << endl;

                        gout << (m_folder + obj_files[i]).toStdString() << " " << gram_files[gram].toStdString() << endl;
                        for (int k = 0; k < optimal_vec.size(); ++k)
                            gout << optimal_vec[k] << " ";

                        gout << endl;

                        ofstream kout;
                        
                        kout.open((m_folder + "ErrorFiles/" + QString::number(a_intval) + "_" + QString::number(b_intval) + ".txt").toStdString());
                        for (int j = 0; j < m_cmaesvec[i]->local_error.size(); ++j)
                        {
                            kout << m_cmaesvec[i]->local_error[j][0];
                            for (int k = 1; k < m_cmaesvec[i]->local_error[j].size(); ++k)
                            {
                                kout << " " << m_cmaesvec[i]->local_error[j][k];
                            }
                            kout << "\n";
                        }
                        kout.close();
                    }

                }
            }
        }, t*num_parallel / num_threads, (t + 1) == num_threads ? num_parallel : (t + 1)*num_parallel / num_threads, t));
    }
    std::for_each(threads.begin(), threads.end(), [](std::thread& x){x.join(); });


   /* Concurrency::parallel_for(
        0, num_parallel, [&](int i)
    {
        QString cloud_name = m_folder + obj_files[i] + ".cld";
        qDebug() << cloud_name;
        vector<vec3> points = PointCloudIO::load_point_cloud(cloud_name);
        // i++; //Rewrite

        for (int gram = 0; gram < gram_files.size(); ++gram)
        {
            qDebug() << i << " " << gram;
            m_grammarvec[i] = new Grammar(m_folder + gram_files[gram]);
            qDebug() << (m_folder + gram_files[gram]);
            vector<double> optimal_vec = test_CMAESGrammarParallel(i, points);
            m_optimal[i*gram_files.size() + gram] = optimal_vec;
            m_errormatrix[i][gram] = m_cmaesvec[i]->m_besterror;

            fout << (m_folder + obj_files[i]).toStdString() << " " << gram_files[gram].toStdString() << endl;
            fout << m_errormatrix[i][gram] << endl;

            gout << (m_folder + obj_files[i]).toStdString() << " " << gram_files[gram].toStdString() << endl;
            for (int k = 0; k < optimal_vec.size(); ++k)
                gout << optimal_vec[k] << " ";

            gout << endl;

            // Remove following lines when you don't want to see output.
            //m_grammar = m_grammarvec[i];
            //qDebug() << m_grammar->m_nrIndependentVariables;
            //const double *solution = optimal_vec.data();
            //m_grammar->softCleanUp();
            //m_grammar->adjustParameters(solution, m_grammar->m_nrIndependentVariables);
            //m_grammar->derive();
            //m_parent->updateGL();
            //m_parent->updateGL();
            //m_parent->updateGL();
            //m_parent->updateGL();
            //m_parent->updateGL();
            //m_parent->updateGL();
            //m_parent->updateGL();

        }
    }
    ); */

    fout.close();
    gout.close();
}

void Scene::generateNewData()
{
    //QString m_folder = "C:/Users/Vignesh/Downloads/benchmark_results/benchmark_results/assembly_chairs/output/35/";
    //vector<vector<vec3>> p = PointCloudIO::load_point_cloud_ply(m_folder + "35_input.ply");
    //vector<vector<double>> camera_pos = load_matrix((m_folder + "occlusion_pose.txt").toStdString());
    //vector<vec4> cp;
    //cp.resize(4);
    //for (int i = 0; i < 4; ++i)
    //{
    //    cp[i].x = camera_pos[i*4][0];
    //    cp[i].y = camera_pos[i*4+1][0];
    //    cp[i].z = camera_pos[i*4+2][0];
    //    cp[i].w = camera_pos[i*4+3][0];
    //}
    //m_verticesScan.clear();
    //for (int i = 0; i < p[0].size(); ++i)
    //{
    //    Vertex v;
    //    v.position = p[0][i];
    //    v.normal = p[1][i];
    //    m_verticesScan.insert(QString::number(i), v);
    //}
    
    //for (int i = 0; i < 4; ++i)
    //    qDebug() << cp[i].x << cp[i].y << cp[i].z << cp[i].w << endl;

    QString m_folder = "C:/Users/Vignesh/Downloads/for_vig/for_vig/table/";
    vector<vector<double>> p = load_matrix((m_folder+"scan.off").toStdString());
    m_verticesScan.clear();
    for (int i = 0; i < p.size(); ++i)
    {
        Vertex v;
        v.position = vec4(p[i][0], p[i][1]+0.04 , p[i][2], 1.0);
        m_verticesScan.insert(QString::number(i), v);
    }
    
}

void Scene::saveAllImages()
{
    QString m_prefix = "C:/Users/Vignesh/Downloads/objects/TestSet/chair/chair_scan_0002";
    QString m_suffix = ".obj";
    for (int i = 100; i < 478; ++i)
    {
        load_object(m_prefix + QString::number(i) + m_suffix);
        m_parent->updateGL();
        m_parent->updateGL();
        m_parent->updateGL();
        saveScene(QString::number(i+21574));
    }
}

void Scene::load_object()
{
#ifdef WIN32
   load_object("C:/Users/Vignesh/Dropbox/reqd/chair0515/model.obj");
#else
    load_object("./../_main/Data/Objs/chair_04.obj");
#endif
}

void Scene::load_object(QString path)
{
    delete m_object;
    // Translation, Scaling and Rotation (lateral, vertical = y axis) information provided
    m_object = new Object(path, true, true, true, vec3(0.0f, 0.0f, 0.0f), vec3(1.0, 1.0, 1.0), vec4(90.0f, 0.0f, 1.0f, 0.0f), vec4(1.0f, 1.0f, 1.0f, 1.0f));
}

void Scene::load_object(QString path, bool normalize)
{
    delete m_object;
    m_object = new Object(path, normalize, true, true, vec3(0.0f, 0.0f, 0.0f), vec3(1.0, 1.0, 1.0), vec4(310.0f, 0.0f, 1.0f, 0.0f), vec4(1.0f, 1.0f, 1.0f, 1.0f));
    //m_object = new Object(path, true, true, true, vec3(2.05f, 0.69f, 2.09f), vec3(1.0, 1.0, 1.0), vec4(90.0f, 0.0f, 1.0f, 0.0f), vec4(1.0f, 1.0f, 1.0f, 1.0f));
    //m_object = new Object(path, true, true, true, vec3(0.0f, 0.50f, 0.0f), vec3(1.0, 1.0, 1.0), vec4(90.0f, 1.0f, 0.0f, 0.0f), vec4(1.0f, 1.0f, 1.0f, 1.0f));
}

void Scene::obtain_gram()
{
    // Given a set of template parameters, this is used to load it to that particular template.
    vector<vector<double>> m = load_matrix("C:/Users/Vignesh/Downloads/objects/TestSet/desk/desk_scan_0005_4.txt");
 
    m_grammar->softCleanUp();
    m_grammar->adjustParameters(m[0].data(), m_grammar->determineNrIndependentVariables());
    m_grammar->derive();
}

vector<vector<double>> Scene::load_matrix(string filename)
{
    ifstream fin;
    string line;
    fin.open(filename, std::ifstream::in);
    vector<vector<double>> matrix;
    getline(fin, line, '\n');
    int i = 1;
    while (!fin.eof())
    {
        vector<double> matrix_row;

        std::stringstream ss(line);
        double myDouble = 0;
        while (ss >> myDouble)
            matrix_row.push_back(myDouble);

        matrix.push_back(matrix_row);
        getline(fin, line, '\n');

    }

    return matrix;
}

void Scene::initializeProceduralModel()
{
	//Test Proc_Model
	m_pmodel = new Proc_Model();
	//p.new_proc(num_boxes);
	vector<vec3> targetPositions, targetSizes;
	targetPositions.push_back(vec3(1, 5, 3));
	targetPositions.push_back(vec3(7, 5, 3));
	targetPositions.push_back(vec3(8, 2.6, 3));
	targetSizes.push_back(vec3(6, 1, 0.2));
	targetSizes.push_back(vec3(4, 2.2, 0.1));
	targetSizes.push_back(vec3(3, 2.4, 2.1));
	m_pmodel->new_proc(targetPositions, targetSizes);

	buildVBOTarget(targetPositions, targetSizes);
}

void Scene::optimize()
{
    if (params::inst()->optTechnique == 0)
        test_CMAESGrammar();
    else if (params::inst()->optTechnique == 1)
        test_DESolverGrammar();
    else if (params::inst()->optTechnique == 2)
        test_Simulated_Anneal();
    else if (params::inst()->optTechnique == 3)
        test_Swarm();
    else test_Genetic();
}

void Scene::test_Genetic()
{
    delete m_genetic;
    int population = params::inst()->optPopulation;
    int nr = m_grammar->determineNrIndependentVariables();
    int iterations = params::inst()->optIterations;
    vector<double> mi, ma;
    mi.resize(nr);
    ma.resize(nr);
    for (int i = 0; i < nr; i++)
    {
        mi[i] = 0.0;
        ma[i] = 2.0;
    }
    m_genetic = new Genetic(this, iterations, params::inst()->optBits, population, mi, ma, params::inst()->optCrossover, params::inst()->optMutation);

    QList<Vertex> list = m_verticesScan.values();
    vector<vec3> samplevals;
    for (auto iter = list.begin(); iter != list.end(); ++iter)
        samplevals.push_back(iter->position);

    m_genetic->setSamples(samplevals);
    m_genetic->setGrammar(m_grammar);

    m_genetic->optimize();

    const double *solution = m_genetic->global_opt().data();
}

void Scene::test_Swarm()
{
    delete m_swarm;
    int num_particles = params::inst()->optParticles;
    int nr = m_grammar->determineNrIndependentVariables();
    int iterations = params::inst()->optIterations;
    vector<double> mi, ma;
    double w = params::inst()->optW , phi_p = params::inst()->optphip, phi_g = params::inst()->optphig;
    mi.resize(nr);
    ma.resize(nr);
    for (int i = 0; i < nr; i++)
    {
        mi[i] = 0.0;
        ma[i] = 2.0;
    }
    m_swarm = new Swarm(this, m_grammar, num_particles, iterations, mi, ma, w, phi_p, phi_g);

    QList<Vertex> list = m_verticesScan.values();
    vector<vec3> samplevals;
    for (auto iter = list.begin(); iter != list.end(); ++iter)
        samplevals.push_back(iter->position);
    
    m_swarm->setSamples(samplevals);
    m_swarm->setGrammar(m_grammar);

    m_swarm->optimize();

    const double *solution = m_swarm->global_opt().data();
}

void Scene::test_Simulated_Anneal()
{
    delete m_sim_anneal;

    int nr = m_grammar->determineNrIndependentVariables();
    m_sim_anneal = new Sim_Anneal(nr, this);

    QList<Vertex> list = m_verticesScan.values();
    vector<vec3> samplevals;
    for (auto iter = list.begin(); iter != list.end(); ++iter)
        samplevals.push_back(iter->position);

    m_sim_anneal->setSampleVals(samplevals);
    m_sim_anneal->setGrammar(m_grammar);

    const double *solution = m_sim_anneal->optimize().data();
}

void Scene::test_DESolverProcModel()
{
	delete m_optimizer;

	const int num_boxes = 3, nr = 12;// nr = 6*num_boxes;
    m_optimizer = new Optimizer(this, params::inst()->optScale, params::inst()->optProbability, params::inst()->optIterations, nr);
    //connect(m_optimizer, SIGNAL(bestCurrentSolution(double)), this, SLOT(updateGrammarGeometry(double)));

	vector<vec3> positions = m_pmodel->position();
	vector<vec3> sizes = m_pmodel->size();
	vector<vec3> samplevals;
	std::default_random_engine unif_generator;
	std::uniform_real_distribution<double> unif_distribution(0.0, 1.0);

	for (int i = 0; i < 200; i++)
	{
		double d1 = unif_distribution(unif_generator);
		double d2 = unif_distribution(unif_generator);
		vec3 samplevals_new;
		for (int j = 0; j < num_boxes; j++)
		{
			vec3 pos = positions[j], size = sizes[j];

			//Bottom face
			samplevals_new.z = pos.z;
			samplevals_new.x = pos.x + d1*size.x;
			samplevals_new.y = pos.y + d2*size.y;
			samplevals.push_back(samplevals_new);

			//Top face
			samplevals_new.z += size.z;
			samplevals.push_back(samplevals_new);

			//Side face x-
			samplevals_new.x = pos.x;
			samplevals_new.y = pos.y + d1*size.y;
			samplevals_new.z = pos.z + d2*size.z;
			samplevals.push_back(samplevals_new);

			//Side face x+
			samplevals_new.x += size.x;
			samplevals.push_back(samplevals_new);

			//Side face y-
			samplevals_new.y = pos.y;
			samplevals_new.x = pos.x + d2*size.x;
			samplevals_new.z = pos.z + d1*size.z;
			samplevals.push_back(samplevals_new);

			//Side face y+
			samplevals_new.y += size.y;
			samplevals.push_back(samplevals_new);
		}
	}

	m_pmodel->create_default_distance_field();

	vector<double> targetDf = m_pmodel->df();//create_point_cloud_distance_field(samplevals);

	vector<double> mi, ma;

	for (int i = 0; i<nr; ++i)
	{
		mi.push_back(0.0);
		ma.push_back(16.0);
	}

	m_optimizer->setRange(mi, ma);
	m_optimizer->setDistanceField(targetDf);
	m_optimizer->setSampleVals(samplevals);
	//m_optimizer->setGrammar(grammar);    

	//Run Optimization
	const double *solution = m_optimizer->solve();

	//Output Target
	for (int i = 0; i < num_boxes; i++)
        cout << positions[i].x << " " << positions[i].y << " " << positions[i].z << " " << sizes[i].x << " " << sizes[i].y << " " << sizes[i].z << endl;

	//Output Solution
	printf("\n\n");
	for (int i = 0; i < nr; ++i)
	{
		printf("%f ", solution[i]);
	}

	vector<vec3> solution_positions, solution_sizes;
	//for (int i = 0; i < nr / 6; i++)
	for (int i = 0; i < num_boxes; i++)
	{
		//solution_positions.push_back(vec3(solution[6 * i], solution[6 * i + 1], solution[6 * i + 2]));
		//solution_sizes.push_back(vec3(solution[6 * i + 3], solution[6 * i + 4], solution[6 * i + 5]));

		if (i == 0)
		{
			solution_positions.push_back(vec3(solution[0], solution[1], solution[2]));
			solution_sizes.push_back(vec3(solution[3], solution[4], solution[5]));
		}
		else if (i == 1)
		{
			solution_positions.push_back(vec3(solution_positions[0].x + solution_sizes[0].x, solution_positions[0].y, solution_positions[0].z));
			solution_sizes.push_back(vec3(solution[6], solution[7], solution[8]));
		}
		else if (i == 2)
		{
			solution_positions.push_back(vec3(solution_positions[1].x + solution_sizes[1].x - solution[9], solution_positions[1].y - solution[10], solution_positions[1].z));
			solution_sizes.push_back(vec3(solution[9], solution[10], solution[11]));
		}
	}

	buildVBOSolution(solution_positions, solution_sizes);
}

void Scene::test_DESolverGrammar()
{
	delete m_optimizer;	

	vector<double> mi, ma;

	int nr = m_grammar->determineNrIndependentVariables();
    m_optimizer = new Optimizer(this, params::inst()->optScale, params::inst()->optProbability, params::inst()->optIterations, nr);

	for (int i = 0; i<nr; ++i)
	{
        mi.push_back(-5.0);
        ma.push_back(5.0);
	}

	m_optimizer->setRange(mi, ma);
	
    QList<Vertex> list= m_verticesScan.values();
	vector<vec3> samplevals;
	for (auto iter = list.begin(); iter != list.end(); ++iter)
		samplevals.push_back(iter->position);	

	m_optimizer->setSampleVals(samplevals);
	m_optimizer->setGrammar(m_grammar);

	const double *solution = m_optimizer->solve();
    
	//buildVBOSolution(solution_positions, solution_sizes);

    
}

void Scene::test_CMAESGrammar()
{
    delete m_cmaes;

    vector<double> mi, ma;

    int nr = m_grammar->determineNrIndependentVariables();
    m_cmaes = new CMAESWrapper(nr, this);

    QList<Vertex> list = m_verticesScan.values();
    vector<vec3> samplevals;
    for (auto iter = list.begin(); iter != list.end(); ++iter)
        samplevals.push_back(iter->position);

    m_cmaes->setSampleVals(samplevals);
    m_cmaes->setGrammar(m_grammar);
    m_cmaes->setVolume(m_shapevolume);

    //const double *solution = m_cmaes->optimize();
    vector<double> solution = m_cmaes->optimize();

    //buildVBOSolution(solution_positions, solution_sizes);
}

vector<double> Scene::test_CMAESGrammarParallel(int i, vector<vec3> samplevals)
{
    
    delete m_cmaesvec[i];

    vector<double> mi, ma;

    int nr = m_grammarvec[i]->determineNrIndependentVariables();
    m_cmaesvec[i] = new CMAESWrapper(nr, this);

    //QList<Vertex> list = m_verticesScan.values();
    //vector<vec3> samplevals;
    //for (auto iter = list.begin(); iter != list.end(); ++iter)
    //    samplevals.push_back(iter->position);

    m_cmaesvec[i]->setSampleVals(samplevals);
    m_cmaesvec[i]->setGrammar(m_grammarvec[i]);
    m_cmaesvec[i]->setVolume(m_shapevolume);

    //m_grammarvec[i]->update();

    //const double *val = m_cmaesvec[i]->optimize();
    vector<double> solution = m_cmaesvec[i]->optimize();

    qDebug() << "CMA-ES Iteration Number: " << i;
    return solution;
    //buildVBOSolution(solution_positions, solution_sizes);
    
}

void Scene::test_ParallelRun(vector<vec3> points)
{
    
    vector<vector<double>> xfinal;
    xfinal.resize(m_num_parallel);
    double dim = m_grammar->determineNrIndependentVariables();
    for (int i = 0; i < m_num_parallel; ++i)
        xfinal[i].resize(dim);

    m_cmaesvec.resize(m_num_parallel);
    //memset(m_cmaesvec, nullptr, m_num_parallel);
    for (int i = 0; i < m_num_parallel; ++i)
        m_cmaesvec[i] = nullptr;

    Concurrency::parallel_for(
        0, m_num_parallel, [&](int i)
    {
        vector<double> sol_i = test_CMAESGrammarParallel(i, points);
        for (int j = 0; j < dim; ++j)
            xfinal[i][j] = sol_i[j];
    }
        );

    int bestindex = 0;
    m_error = 1000.0;
    for (int i = 0; i < m_num_parallel; ++i)
    {
        if (m_cmaesvec[i]->m_besterror < m_error)
        {
            bestindex = i;
            m_error = m_cmaesvec[i]->m_besterror;
        }
    }

    m_grammar = m_grammarvec[bestindex];
    m_optimal.push_back(xfinal[bestindex]);
    const double *solution = xfinal[bestindex].data();
    qDebug() << m_error;
    qDebug() << bestindex;
    m_grammar->softCleanUp();
    m_grammar->adjustParameters(solution, m_grammar->m_nrIndependentVariables);
    m_grammar->derive();

    m_parent->updateGL();
    m_parent->updateGL();
    m_parent->updateGL();
    m_parent->updateGL();
    m_parent->updateGL();
    m_parent->updateGL();
    m_parent->updateGL();
    m_parent->updateGL();
    m_parent->updateGL();
    m_parent->updateGL();
    m_parent->updateGL();
    m_parent->updateGL();

}

void Scene::newGenerate()
{
    if (m_mean.size() == 0)
    {
        m_mean = m_grammar->setParameterVals();
        m_std.resize(m_mean.size());
        for (int i = 0; i < m_mean.size(); ++i)
            m_std[i] = 0.00*sqrt(i);
    }
    for (int i = 0; i < m_std.size(); ++i)
        qDebug() << m_mean[i] << m_std[i];
    
    Proc_Model p;
    QMultiMap<QString, Node*> nodes = m_grammar->nodes();

    vector<vec3> positions, sizes, boxEnds;
    vector<int> colors;
    for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
    {
        Node *n = iter.value();
        positions.push_back(n->translation());
        sizes.push_back(n->size());
        colors.push_back(n->colorMode());
    }

    for (int i = 0; i < positions.size(); i++)
    {
        boxEnds.push_back(positions[i] - sizes[i] / 2);
        boxEnds.push_back(positions[i] + sizes[i] / 2);
        positions[i] -= sizes[i] / 2;
    }

    p.new_proc(positions, sizes);
    p.setColor(colors);

    qDebug() << "Done 1";

    //m_clusters->OrganizeClusters();
    srand(time(0));
    vector<double> new_params;
    new_params.resize(m_mean.size());
    for (int i = 0; i < m_mean.size(); ++i)
        new_params[i] = (m_mean[i] + m_std[i] * rand(0.0, 1.0));//*rand(0.0,1.0);

    vector<vector<double>> new_params_1 = load_matrix("test/test_params.txt");
    new_params = new_params_1[0];

    m_grammar->softCleanUp();
    //m_grammar->adjustParameters(m_mean.data(), m_grammar->determineNrIndependentVariables());
    m_grammar->adjustParameters(new_params.data(), m_grammar->determineNrIndependentVariables());
    m_grammar->derive();

    Proc_Model q;
    nodes = m_grammar->nodes();

    positions.clear();
    sizes.clear();
    boxEnds.clear();
    for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
    {
        Node *n = iter.value();
        positions.push_back(n->translation());
        sizes.push_back(n->size());
    }

    for (int i = 0; i < positions.size(); i++)
    {
        boxEnds.push_back(positions[i] - sizes[i] / 2);
        boxEnds.push_back(positions[i] + sizes[i] / 2);
        positions[i] -= sizes[i] / 2;
    }

    q.new_proc(positions, sizes);
    q.setColor(colors);

    qDebug() << "Done 2";

    vector<Vertex> vertices = m_object->transform(p, q);
    //QString fname = "test/lamp87_skinned.obj";
    //m_object->saveObj(fname, vertices);
    m_object->buildVBOMeshComplete(vertices);

    qDebug() << "Done 3";
}

void Scene::skinGenerate()
{
    QString m_folder = "C:/Users/Vignesh/Documents/shape-grammars/Data/Grammars/chairs/chairs/chairs/";
    QString m_grammar_folder = "D:/Source/ShapeGrammars/code/_main/Data/Grammars/chair/";
    vector<vector<double>> optimals = load_matrix((m_folder + "chairs_all.txt").toStdString());
    vector<vector<double>> endpts = load_matrix((m_folder + "endpts_chairs.txt").toStdString());

    int shape_number;
    std::cout << "Enter a shape number: " << std::endl;
    std::cin >> shape_number;
    shape_number--;

    int req_gram;
    std::cout << "Enter the template you're interested in: " << std::endl;
    std::cin >> req_gram;
    req_gram--;

    m_grammar = new Grammar(m_grammar_folder + "chair_gram" + QString::number(req_gram+1) + ".gram");

    int start_val = endpts[req_gram][0], end_val = endpts[req_gram][1];

    if (shape_number > optimals.size()-1)
        std::cout << "Sorry, this is beyond the available limit.";
    else
    {
        QString file1 = m_folder + "model (" + QString::number(shape_number + 1) + ")/model.obj";
        
        load_object(file1);
        //
        //qDebug() << "Number of vertices are: " << m_object->m_vertices.size();

        //QString file1;
        //if (shape_number + 1 < 10)
        //    file1 = m_folder + "model000" + QString::number(shape_number + 1) + ".obj";
        //else if (shape_number + 1 < 100)
        //    file1 = m_folder + "model00" + QString::number(shape_number + 1) + ".obj";
        //else if (shape_number + 1 < 1000)
        //    file1 = m_folder + "model0" + QString::number(shape_number + 1) + ".obj";
        //else file1 = m_folder + "model" + QString::number(shape_number + 1) + ".obj";
        //load_object(file1);

        vector<double> original_params;

        for (int i = start_val; i <= end_val; ++i)
        {
            original_params.push_back(optimals[shape_number][i - 1]);
        }

        //m_grammar = new Grammar(m_grammar_folder + "desk_gram4.gram");
        //vector<vector<double>> orig = load_matrix("Paper_Data/Matthias_Comparison/desk351_4.txt");
        //original_params = orig[0];

        m_grammar->softCleanUp();
        m_grammar->adjustParameters(original_params.data(), m_grammar->determineNrIndependentVariables());
        m_grammar->derive();

        Proc_Model p;
        QMultiMap<QString, Node*> nodes = m_grammar->nodes();

        vector<vec3> positions, sizes, boxEnds;
        vector<int> colors;
        for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
        {
            Node *n = iter.value();
            positions.push_back(n->translation());
            sizes.push_back(n->size());
            colors.push_back(n->colorMode());
        }

        for (int i = 0; i < positions.size(); i++)
        {
            boxEnds.push_back(positions[i] - sizes[i] / 2);
            boxEnds.push_back(positions[i] + sizes[i] / 2);
            positions[i] -= sizes[i] / 2;
        }

        p.new_proc(positions, sizes);
        p.setColor(colors);

        //string filename;
        //std::cout << "Input file name: ";
        //std::cin >> filename;
        //vector<vector<double>> new_params_1 = load_matrix(filename);
        //vector<vector<double>> new_params_1 = load_matrix("test/test_params.txt");
        vector<vector<double>> new_params_1 = load_matrix("C:/Users/Vignesh/Downloads/objects/TestSet/chair/chair_scan_0001_6.txt");
        vector<double> new_params = new_params_1[0];

        m_grammar->softCleanUp();
        //m_grammar->adjustParameters(m_mean.data(), m_grammar->determineNrIndependentVariables());
        m_grammar->adjustParameters(new_params.data(), m_grammar->determineNrIndependentVariables());
        m_grammar->derive();

        Proc_Model q;
        nodes = m_grammar->nodes();

        positions.clear();
        sizes.clear();
        boxEnds.clear();
        for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
        {
            Node *n = iter.value();
            positions.push_back(n->translation());
            sizes.push_back(n->size());
        }

        for (int i = 0; i < positions.size(); i++)
        {
            boxEnds.push_back(positions[i] - sizes[i] / 2);
            boxEnds.push_back(positions[i] + sizes[i] / 2);
            positions[i] -= sizes[i] / 2;
        }

        q.new_proc(positions, sizes);
        q.setColor(colors);

        vector<Vertex> vertices = m_object->transform(p, q);
        
        QString fname = "C:/Users/Vignesh/Downloads/objects/TestSet/chair/Results/chair"+QString::number(shape_number+1)+"_skinned.obj";
        m_object->saveObj(fname, vertices);
        m_object->buildVBOMeshComplete(vertices);

        m_object->m_rotation = vec4(315.0f, 0.0f, 1.0f, 0.0f);
        m_object->m_vbosLines.clear();
        m_object->m_vbosTriangles.clear();
        m_object->m_vbosNormals.clear();
    }

}

void Scene::renderOptimal()
{
    ifstream fin_gt, fin_predictions, fin_endpts, fin_grams;
    vector<vector<double>> all_optimals;
    QString m_folder = "C:/Users/Vignesh/Documents/shape-grammars/Data/Grammars/desk/";
    QString m_grammar_folder = "D:/Source/ShapeGrammars/code/_main/Data/Grammars/desk/";
    vector<vector<double>> endpts = load_matrix((m_folder + "endpts_desks.txt").toStdString());
    vector<vector<double>> optimals = load_matrix((m_folder + "desks_all.txt").toStdString());

    int shape_number, grammar_number;
    std::cout << "Enter a shape number: " << std::endl;
    std::cin >> shape_number;
    shape_number--;

    std::cout << "Enter a grammar number: " << std::endl;
    std::cin >> grammar_number;
    grammar_number--;

    if (shape_number > optimals.size())
        std::cout << "Sorry, this is beyond the limit available." << std::endl;
    else
    {
        vector<double> gt_dims;
        int start_val = endpts[grammar_number][0], end_val = endpts[grammar_number][1];
        for (int i = start_val; i <= end_val; ++i)
            gt_dims.push_back(optimals[shape_number][i-1]);

        m_grammar = new Grammar(m_grammar_folder + "desk_gram" + QString::number(grammar_number + 1) + ".gram");
        m_grammar->softCleanUp();
        m_grammar->adjustParameters(gt_dims.data(), m_grammar->determineNrIndependentVariables());
        m_grammar->derive();
        QString file1;
        if (shape_number+1 < 10)
            file1 = m_folder + "model000" + QString::number(shape_number+1) + ".obj";
        else if (shape_number+1 < 100)
            file1 = m_folder + "model00" + QString::number(shape_number + 1) + ".obj";
        else if (shape_number+1 < 1000)
            file1 = m_folder + "model0" + QString::number(shape_number + 1) + ".obj";
        else file1 = m_folder + "model" + QString::number(shape_number + 1) + ".obj";
        load_object(file1);
                
    }
}

void Scene::renderOptimalTrue()
{
    ifstream fin_gt, fin_predictions, fin_endpts, fin_grams;
    vector<vector<double>> all_optimals;
    QString m_folder = "C:/Users/Vignesh/Documents/shape-grammars/Data/Grammars/chairs/chairs/chairs/";
    QString m_grammar_folder = "D:/Source/ShapeGrammars/code/_main/Data/Grammars/chair/";
    vector<vector<double>> predictions = load_matrix((m_folder + "chair_predictions_2.txt").toStdString());
    vector<vector<double>> endpts = load_matrix((m_folder + "endpts.txt").toStdString());
    vector<vector<double>> optimals = load_matrix((m_folder + "chair_optimals_2.txt").toStdString());
    vector<vector<double>> grams = load_matrix((m_folder + "test_grams_2.txt").toStdString());

    int shape_number;
    std::cout << "Enter a shape number: " << std::endl;
    std::cin >> shape_number;
    shape_number--;

    if (shape_number > grams.size())
        std::cout << "Sorry, this is beyond the limit available." << std::endl;
    else
    {
        int actual_shape = grams[shape_number][0];
        int grammar_shape = grams[shape_number][1];

        vector<double> prediction_dims, gt_dims;
        int start_val = endpts[grammar_shape - 1][0], end_val = endpts[grammar_shape - 1][1];
        for (int i = start_val; i <= end_val; ++i)
        {
            prediction_dims.push_back(predictions[shape_number][i - 1]);
            gt_dims.push_back(optimals[shape_number][i - 1]);
        }

        m_grammar = new Grammar(m_grammar_folder + "chair_gram" + QString::number(grammar_shape) + ".gram");
        std::cout << "Enter 1 for ground truth and 2 for prediction: ";
        int dec;
        std::cin >> dec;
        m_grammar->softCleanUp();
        if (dec == 1)
            m_grammar->adjustParameters(gt_dims.data(), m_grammar->determineNrIndependentVariables());
        else m_grammar->adjustParameters(prediction_dims.data(), m_grammar->determineNrIndependentVariables());
        m_grammar->derive();
        QString file1 = m_folder + "model (" + QString::number(actual_shape) + ")/model.obj";
        load_object(file1);
    }
}

void Scene::generateNewDatabase()
{
    ifstream fin_gt, fin_predictions, fin_endpts, fin_grams;
    vector<vector<double>> all_optimals;
    QString m_folder = "C:/Users/Vignesh/Documents/shape-grammars/Data/Grammars/sofas_all/";
    QString m_grammar_folder = "D:/Source/ShapeGrammars/code/_main/Data/Grammars/sofa/";
    vector<vector<double>> endpts = load_matrix((m_folder + "endpts_sofas.txt").toStdString());
    vector<vector<double>> optimals = load_matrix((m_folder + "sofas_all.txt").toStdString());
    vector<vector<double>> best_gram = load_matrix("C:/Users/Vignesh/Downloads/sofa_grams.txt");

    for (int i = 0; i < 3173; ++i)
    {
        if ((i != 3619) && (i != 5050))
        {
            int shape = i + 1;
            QString file1 = m_folder + "model (" + QString::number(shape) + ").obj";
            int grammar_number = best_gram[i][1] - 1;
            load_object(file1);

            vector<double> gt_dims;
            int start_val = endpts[grammar_number][0], end_val = endpts[grammar_number][1];
            for (int j = start_val; j <= end_val; ++j)
                gt_dims.push_back(optimals[i][j - 1]);

            m_grammar = new Grammar(m_grammar_folder + "sofa_gram" + QString::number(grammar_number + 1) + ".gram");
            m_grammar->softCleanUp();
            m_grammar->adjustParameters(gt_dims.data(), m_grammar->determineNrIndependentVariables());
            m_grammar->derive();

            int K = 0;
            while (K < 10)
            {
                if (m_mean.size() == 0)
                {
                    m_mean = m_grammar->setParameterVals();
                    m_std.resize(m_mean.size());
                    for (int j = 0; j < m_mean.size(); ++j)
                        m_std[j] = 0.2;
                }
               
                Proc_Model p;
                QMultiMap<QString, Node*> nodes = m_grammar->nodes();

                vector<vec3> positions, sizes, boxEnds;
                vector<int> colors;
                for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
                {
                    Node *n = iter.value();
                    positions.push_back(n->translation());
                    sizes.push_back(n->size());
                    colors.push_back(n->colorMode());
                }

                for (int j = 0; j < positions.size(); j++)
                {
                    boxEnds.push_back(positions[j] - sizes[j] / 2);
                    boxEnds.push_back(positions[j] + sizes[j] / 2);
                    positions[j] -= sizes[j] / 2;
                }

                p.new_proc(positions, sizes);
                p.setColor(colors);

                qDebug() << "Done 1";

                //m_clusters->OrganizeClusters();
                srand(time(0));
                vector<double> new_params;
                new_params.resize(m_mean.size());
                for (int j = 0; j < m_mean.size(); ++j)
                    new_params[j] = (m_mean[j] + m_std[j] * rand(0.0, 1.0));//*rand(0.0,1.0);

                //vector<vector<double>> new_params_1 = load_matrix("test/test_params.txt");
                //new_params = new_params_1[0];

                m_grammar->softCleanUp();
                //m_grammar->adjustParameters(m_mean.data(), m_grammar->determineNrIndependentVariables());
                m_grammar->adjustParameters(new_params.data(), m_grammar->determineNrIndependentVariables());
                m_grammar->derive();

                Proc_Model q;
                nodes = m_grammar->nodes();

                positions.clear();
                sizes.clear();
                boxEnds.clear();
                for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
                {
                    Node *n = iter.value();
                    positions.push_back(n->translation());
                    sizes.push_back(n->size());
                }

                for (int j = 0; j < positions.size(); j++)
                {
                    boxEnds.push_back(positions[j] - sizes[j] / 2);
                    boxEnds.push_back(positions[j] + sizes[j] / 2);
                    positions[j] -= sizes[j] / 2;
                }

                q.new_proc(positions, sizes);
                q.setColor(colors);

                qDebug() << "Done 2";

                vector<Vertex> vertices = m_object->transform(p, q);
                QString fname = "Skinned_Shapes/sofa_" + QString::number(i) + "_" + QString::number(K) + "_skinned.obj";
                m_object->saveObj(fname, vertices);
                K++;
                qDebug() << i << K;
            }
        }
    }
}

void Scene::newShapeResize()
{
    //Compute Optimal Boxes for both shapes
    QString m_folder = "C:/Users/Vignesh/Documents/shape-grammars/Data/Grammars/chairs/chairs/chairs/";
    QString m_opts_file = m_folder + "opts.txt"; //"Data/Objs/SpecChairs/simple/test/opts.txt";
    
    vector<vector<double>> grams = load_matrix((m_folder + "chair_grams.txt").toStdString());
    vector<vector<double>> endpts = load_matrix((m_folder + "chair_endpts.txt").toStdString());
    vector<vector<double>> optimals = load_matrix((m_folder + "chair_opts.txt").toStdString());

    //For mug ONLY
    //vector<vector<double>> grams;
    //vector<vector<double>> endpts;
    //vector<double> a;
    //a.push_back(1);
    //a.push_back(9);
    //endpts.push_back(a);
    //for (int i = 1; i < 215; ++i)
    //{
    //    vector<double> b;
    //    b.push_back(i);
    //    b.push_back(1);
    //    grams.push_back(b);
    //}
    //QString grammar = "mug_gram2.gram";

    int ind1, ind2;

    std::cout << "Enter shape 1: ";
    std::cin >> ind1;
    std::cout << "Enter shape 2: ";
    std::cin >> ind2;

    //int ind1 = shape1.right(shape1.lastIndexOf("l")).toInt()-1, ind2 = shape2.right(shape2.lastIndexOf("l")).toInt()-1;
    //int ind1 = (shape1.right(4)).left(3).toInt()-1, ind2 = (shape2.right(4)).left(3).toInt()-1;

    //QString file1 = m_folder + "model (" + QString::number(ind1) + ").obj", file2 = m_folder + "model (" + QString::number(ind2) + ").obj";
    //for chair ONLY
    QString file1 = m_folder + "model (" + QString::number(ind1) + ")/model.obj", file2 = m_folder + "model (" + QString::number(ind2) + ")/model.obj";
    //for remaining guys with 000
    //QString file1, file2;
    //QString num1 = QString::number(ind1), num2 = QString::number(ind2);
    //if (ind1 < 10)
    //    file1 = m_folder + "model000" + num1 + ".obj";
    //else if (ind1 < 100)
    //    file1 = m_folder + "model00" + num1 + ".obj";
    //else if (ind1 < 1000)
    //    file1 = m_folder + "model0" + num1 + ".obj";
    //else file1 = m_folder + "model" + num1 + ".obj";
    //if (ind2 < 10)
    //    file2 = m_folder + "model000" + num2 + ".obj";
    //else if (ind1 < 100)
    //    file2 = m_folder + "model00" + num2 + ".obj";
    //else if (ind1 < 1000)
    //    file2 = m_folder + "model0" + num2 + ".obj";
    //else file2 = m_folder + "model" + num2 + ".obj";

    load_object(file1);
    qDebug() << file1;
    
    vector<double> shape1_dims, shape2_dims;
    int gram1 = grams[ind1-1][1], gram2 = grams[ind2-1][1];
    
    QString grammar = "chair_gram" + QString::number(gram1) + ".gram";
    m_grammar = new Grammar(m_folder + grammar);

    int start_val = endpts[gram1 - 1][0], end_val = endpts[gram1 - 1][1];

    for (int i = 0; i < endpts.size(); ++i)
    {
        qDebug() << endpts[i][0] << endpts[i][1];
    }

    for (int i = start_val; i <= end_val; ++i)
    {
        shape1_dims.push_back(optimals[ind1-1][i - 1]);
        shape2_dims.push_back(optimals[ind2-1][i - 1]);
    }

    qDebug() << "The dimensions of shape 1 are: " << endl;

    for (int i = 0; i < shape1_dims.size(); ++i)
        qDebug() << shape1_dims[i];

    qDebug() << "The dimensions of shape 2 are: " << endl;

    for (int i = 0; i < shape2_dims.size(); ++i)
        qDebug() << shape2_dims[i];
    
    Proc_Model p, q;
    m_grammar->softCleanUp();
    m_grammar->adjustParameters(shape1_dims.data(), m_grammar->determineNrIndependentVariables());
    m_grammar->derive();
    
    vector<vec3> positions, sizes, boxEnds;
    vector<int> colors;
    QMultiMap<QString, Node*> nodes = m_grammar->nodes();
    for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
    {
        Node *n = iter.value();
        positions.push_back(n->translation());
        sizes.push_back(n->size());
        colors.push_back(n->colorMode());
    }
    
    for (int i = 0; i < positions.size(); i++)
    {
        boxEnds.push_back(positions[i] - sizes[i] / 2);
        boxEnds.push_back(positions[i] + sizes[i] / 2);
        positions[i] -= sizes[i] / 2;
    }
    
    p.new_proc(positions, sizes);
    p.setColor(colors);

    positions.clear();
    sizes.clear();
    boxEnds.clear();
    
    m_grammar->softCleanUp();
    m_grammar->adjustParameters(shape2_dims.data(), m_grammar->determineNrIndependentVariables());
    m_grammar->derive();
    nodes = m_grammar->nodes();
    
    for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
    {
        Node *n = iter.value();
        positions.push_back(n->translation());
        sizes.push_back(n->size());
    }
    
    for (int i = 0; i < positions.size(); i++)
    {
        boxEnds.push_back(positions[i] - sizes[i] / 2);
        boxEnds.push_back(positions[i] + sizes[i] / 2);
        positions[i] -= sizes[i] / 2;
    }
    
    q.new_proc(positions, sizes);
    q.setColor(colors);

    m_object->buildVBOMeshComplete(m_object->transform(p, q));
}

void Scene::writeBoxes_grammars()
{
    QString m_folder = "C:/Users/Vignesh/Documents/shape-grammars/Data/Grammars/lamp/";
    QString m_grammar_folder = "D:/Source/ShapeGrammars/code/_main/Data/Grammars/lamp/";
    vector<vector<double>> optimals = load_matrix((m_folder + "lamps_all.txt").toStdString());
    vector<vector<double>> endpts = load_matrix((m_folder + "endpts_lamps.txt").toStdString());

    QMap<int, QString> mp;
    ifstream fin;
    fin.open((m_folder + "lamp_equals.txt").toStdString(), std::ifstream::in);
    string line;
    getline(fin, line, '\n');
    while (!fin.eof())
    {
        int num;
        QString rqd_line = QString::fromStdString(line);
        QStringList l = rqd_line.split(" ");
        num = l[0].toInt();
        mp.insert(num, l[1]);
        qDebug() << num;
        getline(fin, line, '\n');
    }
    
    int gram_type = 2;

    vector<double> vals = endpts[gram_type - 1];
    int start_val = vals[0], end_val = vals[1];
    qDebug() << start_val << end_val;
    
    vector<double> ground_truth;
    ground_truth.resize(end_val - start_val + 1);

    m_grammar = new Grammar(m_grammar_folder + "lamp_gram" + QString::number(gram_type) + ".gram");
    
    for (int i = 0; i < optimals.size(); ++i)
    {
        QString filename = m_folder + QString::number(gram_type) + "/" + QString::number(i) + ".txt";
        ofstream fout;
        string line;
        fout.open(filename.toStdString());

        for (int j = start_val; j <= end_val; ++j)
        {
            ground_truth[j - start_val] = optimals[i][j-1];
        }

        m_grammar->softCleanUp();
        m_grammar->adjustParameters(ground_truth.data(), m_grammar->determineNrIndependentVariables());
        m_grammar->derive();

        QString file1 = m_folder + mp[i+1] + "/model.obj";
        //load_object(file1);

        QMultiMap<QString, Node*> nodes = m_grammar->nodes();
        vector<vec3> positions, sizes, boxEnds;
        for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
        {
            Node *n = iter.value();
            positions.push_back(n->translation());
            sizes.push_back(n->size());
        }

        for (int i = 0; i < positions.size(); i++)
        {
            boxEnds.push_back(positions[i] - sizes[i] / 2);
            boxEnds.push_back(positions[i] + sizes[i] / 2);
            positions[i] -= sizes[i] / 2;
        }

        for (int j = 0; j < positions.size(); ++j)
        {
            fout << positions[j].x << " " << positions[j].y << " " << positions[j].z << endl;
            fout << sizes[j].x << " " << sizes[j].y << " " << sizes[j].z << endl;
        }
        fout.close();

        //m_parent->updateGL();
    }
}

void Scene::writeBoxes()
{
    QMultiMap<QString, Node*> nodes = m_grammar->nodes();

    int write_mode;
    std::cout << "Enter write mode: " << std::endl;
    std::cin >> write_mode;

    QString name;
    (write_mode == 1) ? (name = "test_params") : (name = "test_boxes");
    QString filename = "test/" + name + ".txt";
    ofstream fout;
    fout.open(filename.toStdString());

    if (write_mode == 1)
    {
        vector<double> xstart = m_grammar->setParameterVals();
        for (int i = 0; i < xstart.size(); ++i)
            fout << xstart[i] << " ";

        fout << endl;
    }
    else
    {
        vector<vec3> positions, sizes, boxEnds;
        for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
        {
            Node *n = iter.value();
            positions.push_back(n->translation());
            sizes.push_back(n->size());
        }

        for (int i = 0; i < positions.size(); i++)
        {
            boxEnds.push_back(positions[i] - sizes[i] / 2);
            boxEnds.push_back(positions[i] + sizes[i] / 2);
            positions[i] -= sizes[i] / 2;
        }

        string line;
        for (int i = 0; i < positions.size(); ++i)
        {
            fout << positions[i].x << " " << positions[i].y << " " << positions[i].z << endl;
            fout << sizes[i].x << " " << sizes[i].y << " " << sizes[i].z << endl;
        }
    }


    fout.close();
}

void Scene::obtainBestFitShape()
{
    QMultiMap<QString, Node*> nodes = m_grammar->nodes();

    vector<vec3> positions, sizes, boxEnds;
    for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
    {
        Node *n = iter.value();
        positions.push_back(n->translation());
        sizes.push_back(n->size());
    }

    for (int i = 0; i < positions.size(); i++)
    {
        boxEnds.push_back(positions[i] - sizes[i] / 2);
        boxEnds.push_back(positions[i] + sizes[i] / 2);
        positions[i] -= sizes[i] / 2;
    }

    double err = 100000000.0;
    int gram = 6;
    int best = 0;
    for (int s = 0; s < 6778; ++s)
    {
        double tot = 0.0;
        //qDebug() << ("C:/Users/Vignesh/Documents/shape-grammars/Data/Grammars/chairs/chairs/chairs/" + QString::number(gram) + "/" + QString::number(s) + ".txt");
        vector<vector<double>> k = load_matrix(("C:/Users/Vignesh/Documents/shape-grammars/Data/Grammars/chairs/chairs/chairs/" + QString::number(gram) + "/" + QString::number(s) + ".txt").toStdString());
        
        for (int i = 1; i < k.size(); i += 2)
        {
            tot += (k[i][0] - sizes[(i-1)/2].x)*(k[i][0] - sizes[(i-1)/2].x) + (k[i][1] - sizes[(i-1)/2].y)*(k[i][1] - sizes[(i-1)/2].y) + (k[i][2] - sizes[(i-1)/2].z)*(k[i][2] - sizes[(i-1)/2].z);
        }
        //qDebug() << tot;
        if (tot < err)
        {
            err = tot;
            best = s;
            qDebug() << s;
        }
    }

    qDebug() << "The best fit shape is: " << best;
}

void Scene::mergeShapes(QString shape1, QString shape2)
{

    //Compute Optimal Boxes for both shapes
    QString object1 = shape1, object2 = shape2;
    QString m_folder = "Data/Objs/SpecChairs/simple/test/";

    QString cloud_name1 = m_folder + object1 + ".cld";
    QString cloud_name2 = m_folder + object2 + ".cld";
    vector<vec3> points1 = PointCloudIO::load_point_cloud(cloud_name1);
    vector<vec3> points2 = PointCloudIO::load_point_cloud(cloud_name2);
    // check if the corresponding point cloud files exist
    // if the point clouds do not exist, scan and save them
    m_verticesScan.clear();
    load_object(m_folder + object1 + ".obj");
    if (points1.size() <= 10)
    {
        points1.clear();
        for (int i = 0; i < 10; ++i)
        {
            Postoffice::renderer()->m_scanner->m_scanning = true;
            Postoffice::renderer()->m_scanner->autoScan();
        }
        // save
        QList<Vertex> pList = m_verticesScan.values();
        points1.reserve(pList.size());
        std::transform(pList.begin(), pList.end(), back_inserter(points1),
            [](Vertex V){return V.position; });
        qDebug() << "Saving Point Cloud";
        PointCloudIO::save_point_cloud(points1, cloud_name1);
    }
    else if (points1.size() > 10)
    {
        qDebug() << "No. of points: " << points1.size();
        //Occupancy k(points);
        test_ParallelRun(points1);
    }

    m_verticesScan.clear();
    load_object(m_folder + object2 + ".obj");
    if (points2.size() <= 10)
    {
        points2.clear();
        for (int i = 0; i < 10; ++i)
        {
            Postoffice::renderer()->m_scanner->m_scanning = true;
            Postoffice::renderer()->m_scanner->autoScan();
        }
        // save
        QList<Vertex> pList = m_verticesScan.values();
        points2.reserve(pList.size());
        std::transform(pList.begin(), pList.end(), back_inserter(points2),
            [](Vertex V){return V.position; });
        qDebug() << "Saving Point Cloud";
        PointCloudIO::save_point_cloud(points2, cloud_name2);
    }
    else if (points2.size() > 10)
    {
        qDebug() << "No. of points: " << points2.size();
        //Occupancy k(points);
        test_ParallelRun(points2);
    }

    qDebug() << "Size of m_optimal: " << m_optimal.size();

    //Compute mean and standard deviation of solution space
    m_mean.resize(m_optimal[0].size());
    m_std.resize(m_mean.size());
    for (int i = 0; i < m_optimal.size(); ++i)
    {
        vector<double> optimal_i = m_optimal[i];
        for (int j = 0; j < optimal_i.size(); ++j)
        {
            m_mean[j] += ((optimal_i[j]) / m_optimal.size());
            std::cout << optimal_i[j] << " ";
        }
        std::cout << endl;
    }
    std::cout << "Mean: " << " ";
    for (int i = 0; i < m_mean.size(); ++i)
        std::cout << m_mean[i] << " ";

    for (int i = 0; i < m_optimal.size(); ++i)
    {
        vector<double> optimal_i = m_optimal[i];
        for (int j = 0; j < m_mean.size(); ++j)
        {
            m_std[j] += (((optimal_i[j] - m_mean[j])*(optimal_i[j] - m_mean[j])) / m_optimal.size());
        }
    }

    std::cout << "Std. Deviation: " << " ";
    for (int i = 0; i < m_std.size(); ++i)
    {
        std::cout << m_std[i] << " ";
        //m_std[i] = sqrt(m_std[i]);
        m_std[i] = 0.01;
    }

    //Compute boxes corresponding to both shapes
    Proc_Model p, q;
    m_grammar->softCleanUp();
    m_grammar->adjustParameters(m_optimal[0].data(), m_grammar->m_nrIndependentVariables);
    m_grammar->derive();
    QMultiMap<QString, Node*> nodes = m_grammar->nodes();
    vector<vec3> positions, sizes, boxEnds;
    for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
    {
        Node *n = iter.value();
        positions.push_back(n->translation());
        sizes.push_back(n->size());
    }

    for (int i = 0; i < positions.size(); i++)
    {
        boxEnds.push_back(positions[i] - sizes[i] / 2);
        boxEnds.push_back(positions[i] + sizes[i] / 2);
        positions[i] -= sizes[i] / 2;
    }

    p.new_proc(positions, sizes);

    m_grammar->softCleanUp();
    m_grammar->adjustParameters(m_optimal[1].data(), m_grammar->m_nrIndependentVariables);
    m_grammar->derive();

    nodes = m_grammar->nodes();
    positions.clear(); sizes.clear(); boxEnds.clear();
    for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
    {
        Node *n = iter.value();
        positions.push_back(n->translation());
        sizes.push_back(n->size());
    }

    for (int i = 0; i < positions.size(); i++)
    {
        boxEnds.push_back(positions[i] - sizes[i] / 2);
        boxEnds.push_back(positions[i] + sizes[i] / 2);
        positions[i] -= sizes[i] / 2;
    }
    q.new_proc(positions, sizes);

    vector<int> p_boxes, q_boxes;
    p_boxes.resize(points1.size());
    q_boxes.resize(points2.size());
    vector<vec3> p_pos = p.position(), p_size = p.size();
    for (int i = 0; i < points1.size(); ++i)
    {
        double k = p.project_point_to_box(points1[i], p_pos[0], p_size[0]);
        int box_req = 0;
        for (int box_val = 1; box_val < p.position().size(); ++box_val)
        {
            double new_val = p.project_point_to_box(points1[i], p_pos[box_val], p_size[box_val]);
            if (new_val < k)
            {
                k = new_val;
                box_req = box_val;
            }
        }
        p_boxes[i] = box_req;
    }

    vector<vec3> q_pos = q.position(), q_size = q.size();
    for (int i = 0; i < points2.size(); ++i)
    {
        double k = q.project_point_to_box(points2[i], q_pos[0], q_size[0]);
        int box_req = 0;
        for (int box_val = 1; box_val < q.position().size(); ++box_val)
        {
            double new_val = q.project_point_to_box(points2[i], q_pos[box_val], q_size[box_val]);
            if (new_val < k)
            {
                k = new_val;
                box_req = box_val;
            }
        }
        q_boxes[i] = box_req;
    }

    //Compute new shape parametrized by current shape dimensions
    vector<double> new_params;
    new_params.resize(m_mean.size());
    for (int i = 0; i < m_mean.size(); ++i)
        new_params[i] = (m_mean[i] + m_std[i] * rand(0.0, 1.0));//*rand(0.0,1.0);

    m_grammar->softCleanUp();
    m_grammar->adjustParameters(new_params.data(), m_grammar->determineNrIndependentVariables());
    m_grammar->derive();

    Proc_Model r;
    nodes = m_grammar->nodes();

    positions.clear();
    sizes.clear();
    boxEnds.clear();
    for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
    {
        Node *n = iter.value();
        positions.push_back(n->translation());
        sizes.push_back(n->size());
    }

    for (int i = 0; i < positions.size(); i++)
    {
        boxEnds.push_back(positions[i] - sizes[i] / 2);
        boxEnds.push_back(positions[i] + sizes[i] / 2);
        positions[i] -= sizes[i] / 2;
    }

    r.new_proc(positions, sizes);
    vector<vec3> r_pos = r.position(), r_size = r.size();

    //Generate split for boxes to help merging and neighbors of boxes in grammar.
    vector<int> split;
    m_neighbors.clear();
    m_neighbors.resize(r_pos.size());
    for (int i = 0; i < m_neighbors.size(); ++i)
    {
        m_neighbors[i].resize(r_pos.size());
        for (int j = 0; j < r_pos.size(); ++j)
            m_neighbors[i][j] = 0;
    }
    m_neighbors[0][5] = 1; m_neighbors[5][0] = 1;
    m_neighbors[1][5] = 1; m_neighbors[5][1] = 1;
    m_neighbors[2][5] = 1; m_neighbors[5][2] = 1;
    m_neighbors[3][5] = 1; m_neighbors[5][3] = 1;
    m_neighbors[4][5] = 1; m_neighbors[5][4] = 1;
    split.resize(6);
    split[0] = 0, split[1] = 1, split[2] = 1, split[3] = 1, split[4] = 1, split[5] = 0;

    //Use split to build new point cloud.
    vector<vec3> op_pointcloud;
    vector<int> pc_source, box_source;
    for (int i = 0; i < p_boxes.size(); ++i)
    {
        if (split[p_boxes[i]] == 0)
        {
            vec3 trans_pt;
            trans_pt = (points1[i] - p_pos[p_boxes[i]])*r_size[p_boxes[i]] / p_size[p_boxes[i]] + r_pos[p_boxes[i]];
            op_pointcloud.push_back(trans_pt);
            pc_source.push_back(0);
            box_source.push_back(p_boxes[i]);
        }
    }
    for (int i = 0; i < q_boxes.size(); ++i)
    {
        if (split[q_boxes[i]] == 1)
        {
            vec3 trans_pt;
            trans_pt = (points2[i] - q_pos[q_boxes[i]])*r_size[q_boxes[i]] / q_size[q_boxes[i]] + r_pos[q_boxes[i]];
            op_pointcloud.push_back(trans_pt);
            pc_source.push_back(1);
            box_source.push_back(q_boxes[i]);
        }
    }

    vector<vector<int>> points_ref;
    points_ref.resize(p_pos.size());
    for (int i = 0; i < box_source.size(); ++i)
    {
        points_ref[box_source[i]].push_back(i);
    }
    for (int i = 0; i < 6; ++i)
        qDebug() << "Number of points in box " << i << ":" << points_ref[i].size();

    //Extrapolation from existing information.
    for (int i = 0; i < m_neighbors.size(); ++i)
    {      
        for (int j = i + 1; j < m_neighbors.size(); ++j)
        {
            if ((m_neighbors[i][j] == 1) && (split[i] != split[j]))
            {
                vector<vec3> p1, p2;
                p1.resize(points_ref[i].size());
                p2.resize(points_ref[j].size());
                for (int k = 0; k < p1.size(); ++k)
                    p1[k] = op_pointcloud[points_ref[i][k]];
                for (int k = 0; k < p2.size(); ++k)
                    p2[k] = op_pointcloud[points_ref[j][k]];

                int closest_k = 20;
                double min_1 = dot(p1[0]-p1[1], p1[0]-p1[1]), min_2 = dot(p2[0]-p2[1], p2[0]-p2[1]), min_across = dot(p1[0]-p2[0],p1[0]-p2[0]);
                for (int point1 = 0; point1 < p1.size(); ++point1)
                    for (int point2 = point1 + 1; point2 < p1.size(); ++point2)
                    {
                        double k = dot(p1[point1] - p1[point2], p1[point1] - p1[point2]);
                        if (k < min_1)
                            min_1 = k;
                    }
                for (int point1 = 0; point1 < p2.size(); ++point1)
                    for (int point2 = point1 + 1; point2 < p2.size(); ++point2)
                    {
                        double k = dot(p2[point1] - p2[point2], p2[point1] - p2[point2]);
                        if (k < min_2)
                            min_2 = k;
                    }

                for (int point1 = 0; point1 < p1.size(); ++point1)
                    for (int point2 = 0; point2 < p2.size(); ++point2)
                    {
                        double k = dot(p1[point1] - p2[point2], p1[point1] - p2[point2]);
                        if (k < min_across)
                            min_across = k;
                    }
                qDebug() << "Min_1:" << min_1 << ", Min2: " << min_2 << ", Min_across: " << min_across;
                if ((min_1 < min_across) || (min_2 < min_across))
                {
                    vector<vec3> best_k1, best_k2;
                    /*for (int point1 = 0; point1 < p1.size(); ++point1)
                    {
                        vec3 full_point1 = p1[point1];
                        for (int point2 = 0; point2 < p2.size(); ++point2)
                        {
                            vec3 full_point2 = p2[point2];
                            if (best_k1.size() < closest_k)
                            {
                                best_k1.push_back(full_point1), best_k2.push_back(full_point2);
                                int k = best_k1.size() - 1;
                                while ((k > 0) && (dot(best_k1[k] - best_k2[k], best_k1[k] - best_k2[k]) < dot(best_k1[k - 1] - best_k2[k - 1], best_k1[k - 1] - best_k2[k - 1])))
                                {
                                    vec3 t1, t2;
                                    t1 = best_k1[k], t2 = best_k2[k];
                                    best_k1[k] = best_k1[k - 1], best_k2[k] = best_k2[k - 1];
                                    best_k1[k - 1] = t1, best_k2[k - 1] = t2;
                                    k--;
                                }
                            }
                            else
                            {
                                double target = dot(full_point1 - full_point2, full_point1 - full_point2);
                                int k = closest_k - 1;
                                if (target < dot(best_k1[k] - best_k2[k], best_k1[k] - best_k2[k]))
                                {
                                    while (k > 0)
                                    {
                                        best_k1[k] = best_k1[k - 1];
                                        best_k2[k] = best_k2[k - 1];
                                        k--;
                                    }
                                    best_k1[k] = full_point1, best_k2[k] = full_point2;
                                }
                            }
                        }
                    }*/

                    best_k1.resize(closest_k), best_k2.resize(closest_k);
                    int count = 0;
                    vector<vec3> p1_dup = p1, p2_dup = p2;
                    while (count < closest_k)
                    {
                        double best_val = dot(p1_dup[0] - p2_dup[0], p1_dup[0] - p2_dup[0]);
                        int min_ind_1 = 0, min_ind_2 = 0;
                        for (int point1 = 0; point1 < p1_dup.size(); ++point1)
                        {
                            for (int point2 = 0; point2 < p2_dup.size(); ++point2)
                            {
                                if (dot(p1_dup[point1] - p2_dup[point2], p1_dup[point1] - p2_dup[point2]) < best_val)
                                {
                                    min_ind_1 = point1;
                                    min_ind_2 = point2;
                                    best_val = dot(p1_dup[point1] - p2_dup[point2], p1_dup[point1] - p2_dup[point2]);
                                }
                            }
                        }
                        best_k1[count] = p1_dup[min_ind_1];
                        best_k2[count] = p2_dup[min_ind_2];
                        p1_dup.erase(p1_dup.begin() + min_ind_1);
                        p2_dup.erase(p2_dup.begin() + min_ind_2);
                        ++count;
                    }

                    qDebug() << "Best is: " << dot(best_k1[0] - best_k2[0], best_k1[0] - best_k2[0]);
                    
                    for (int k = 0; k < best_k1.size(); ++k)
                    {
                        vec3 fp1 = best_k1[k];
                        //for (int l = 0; l < best_k2.size(); ++l)
                        //{
                        //    vec3 fp2 = best_k2[l];
                            vec3 fp2 = best_k2[k];
                            int tot = 3;
                            for (int m = 1; m < tot; ++m)
                            {
                                vec3 pt = (1.0/tot)*(fp1*m + fp2*(tot - m));
                                op_pointcloud.push_back(pt);
                                pc_source.push_back(2);
                            }
                        //}
                    }
                }
            }
        }
    }

    vector<Vertex> vertices;
    vertices.resize(op_pointcloud.size());
    for (int i = 0; i < vertices.size(); ++i)
    {
        vertices[i].position = op_pointcloud[i];
        if (pc_source[i]==0)
            vertices[i].color = vec4(1.0, 0.0, 0.0, 1.0);
        else if(pc_source[i] ==1)
            vertices[i].color = vec4(0.0, 1.0, 0.0, 1.0);
        else vertices[i].color = vec4(0.0, 0.0, 1.0, 1.0);
    }

    m_object->buildVBOMeshComplete(vertices);
}

void Scene::saveScene(QString path)
{
    m_parent->updateGL();
    m_parent->updateGL();
    m_parent->updateGL();
    saveFrameBuffer(m_parent, path);
}

void Scene::testIndependentVariables()
{
	qDebug() << "-----------";
	int nrIntVars = m_grammar->determineNrIndependentVariables();

	double *trial = new double[nrIntVars];

	for (int i = 0; i < nrIntVars; ++i)
	{
		trial[i] = 0.0f;
	}

	m_grammar->adjustParameters(trial, nrIntVars);
	m_grammar->printParameters();
}

void Scene::updateGrammarGeometry()
{
    m_grammar->softCleanUp();
    m_grammar->derive();

    m_parent->updateGL();
}

void Scene::updateGrammarGeometry_CMAES()
{
    m_grammar->softCleanUp();
    m_grammar->derive();

    m_parent->updateGL();
}

void Scene::updateGrammarGeometry_SimAnneal()
{
    m_grammar->softCleanUp();
    m_grammar->derive();

    m_parent->updateGL();
}

void Scene::updateGrammarGeometry_Swarm()
{
    m_grammar->softCleanUp();
    m_grammar->derive();

    m_parent->updateGL();
}

void Scene::renderWorld(const Transform &trans)
{
    m_niceGrid->render(trans);

    if(params::inst()->renderMisc)
    {	         
        for(int i=0; i<m_lights.size(); ++i)
        {   
            m_lights[i]->render(trans);
        }

        m_cameraManager->renderCameras(trans);
    }
}

void Scene::renderObjects(const Transform &trans)
{
    if (params::inst()->renderCurrent)
	    m_grammar->render(trans);

    mat4 model = mat4::translate(0.0f, 0.0f, 0.0f) *mat4::scale(vec3(0.5));

    /*Target Rendering
	//Shader *shader = shaders::inst()->defaultLight;
	//shader->bind();

	//	shader->set3f("lightPos", params::inst()->lights[0]->position());
	//	shader->setMatrices(trans, model, true, true, true);

	//	if (m_vboSolution && params::inst()->renderCurrent)
	//		m_vboSolution->render();

	//	if (m_vboTarget && params::inst()->renderTarget)
	//		m_vboTarget->render();

	//shader->release();
    */
    
    // Voxel Scanner Rendering 
    /*glPointSize(3.0f);
    Shader *shader = shaders::inst()->defaultLight;
    shader->bind();

    	shader->set3f("lightPos", params::inst()->lights[0]->position());
    	shader->setMatrices(trans, model, true, true, true);

        m_vboVoxels->render();

    shader->release();



    shader = shaders::inst()->default;
    shader->bind();

        shader->set3f("lightPos", params::inst()->lights[0]->position());
        shader->setMatrices(trans, model, true, true, true);

        m_vboVoxelsLines->render();

    shader->release();


    m_activeLayer = 64;*/
	
    //m_activeLayer++;
    //
    //if (m_activeLayer >= 64)
    //    m_activeLayer = 0;
    //
    //buildVBOVoxels();
    
   
    m_object->render(trans);

    if (m_vboModelDifference)
    {
        mat4 model = mat4::translate(m_object->m_position) * mat4::rotateY(m_object->m_rotation.y) * mat4::scale(m_object->m_scale);
        Shader *shader = shaders::inst()->defaultLight;
        shader->bind();
            shader->setMatrices(trans, model, true, true, true);
            shader->set3f("lightPos", params::inst()->lights[0]->position());
            m_vboModelDifference->render();
        shader->release();
    }

    mat4 modelEraser = mat4::translate(m_ballEraserPos);
    Shader *shader = shaders::inst()->default;
    shader->bind();
        shader->setMatrices(trans, modelEraser, true, true, true);
        m_vboBallEraser->render();        
    shader->release();

    //m_clusters->render(trans);
}

void Scene::renderObjectsScan(const Transform &trans)
{    
	m_object->renderScan(trans);
}

void Scene::renderObjectsDepth(const Transform &trans)
{
    m_object->renderDepth(trans);

    if (params::inst()->renderCurrent)
        m_grammar->renderDepth(trans);

   // if (params::inst()->renderGenerations)
   //     m_object->renderDepth(trans);


}
 
void Scene::update(float delta)
{
    for(int i=0; i<m_lights.size(); ++i)
    {
        m_lights[i]->update(delta);
    }
}

void Scene::select(const Transform &trans, int sw, int sh, int mx, int my)
{
    float minDist = math_maxfloat;
    int idx = -1;
    
    Picking pick;
    for(int i=0; i<m_objects.size(); ++i)
    {
        float t = m_objects[i]->selected(pick, trans, sw, sh, mx, my);
        if( t > 0.0f && t<minDist)
        {
            minDist = t;
            idx = i;
        }
    }

    if(idx >= 0)
    {
        m_activeIdx = idx;
        m_objects[m_activeIdx]->m_isSelected = true;
    }
}

void Scene::move(const Transform &trans, int x, int y)
{
	vec3 dir, right, up, pos;
    getCameraFrame(trans, dir, up, right, pos);

    if(m_activeIdx >= 0)
    {
        m_objects[m_activeIdx]->move(x, y, dir, up, right, pos);
    }
}

void Scene::resetSelection()
{
    for(int i=0; i<m_objects.size(); ++i)
    {
        m_objects[i]->m_isSelected = false;
    }

    m_activeIdx = -1;
}

void Scene::buildVBOSolution(const vector<vec3> &positions, const vector<vec3> &sizes)
{
	delete m_vboSolution;

	vector<vec3> vertices, normals;
	for (int i = 0; i < positions.size(); ++i)
	{
		vec3 p = positions[i];
		vec3 s = sizes[i];

		vec3 p1 = p;
		vec3 p2 = p + vec3(s.x, 0.0f, 0.0f);
		vec3 p3 = p + vec3(s.x, s.y, 0.0f);
		vec3 p4 = p + vec3(0.0f, s.y, 0.0f);
		vec3 p5 = p + vec3(0.0f, 0.0f, s.z);
		vec3 p6 = p + vec3(s.x, 0.0f, s.z);
		vec3 p7 = p + s;
		vec3 p8 = p + vec3(0.0f, s.y, s.z);

		vertices.push_back(p1);
		vertices.push_back(p2);
		vertices.push_back(p3);
		vertices.push_back(p4);

		normals.push_back(vec3(0.0f, 0.0f, -1.0f));
		normals.push_back(vec3(0.0f, 0.0f, -1.0f));
		normals.push_back(vec3(0.0f, 0.0f, -1.0f));
		normals.push_back(vec3(0.0f, 0.0f, -1.0f));

		vertices.push_back(p1);
		vertices.push_back(p2);
		vertices.push_back(p6);
		vertices.push_back(p5);

		normals.push_back(vec3(0.0f, -1.0f, 0.0f));
		normals.push_back(vec3(0.0f, -1.0f, 0.0f));
		normals.push_back(vec3(0.0f, -1.0f, 0.0f));
		normals.push_back(vec3(0.0f, -1.0f, 0.0f));

		vertices.push_back(p2);
		vertices.push_back(p3);
		vertices.push_back(p7);
		vertices.push_back(p6);

		normals.push_back(vec3(1.0f, 0.0f, 0.0f));
		normals.push_back(vec3(1.0f, 0.0f, 0.0f));
		normals.push_back(vec3(1.0f, 0.0f, 0.0f));
		normals.push_back(vec3(1.0f, 0.0f, 0.0f));

		vertices.push_back(p1);
		vertices.push_back(p4);
		vertices.push_back(p8);
		vertices.push_back(p5);

		normals.push_back(vec3(-1.0f, 0.0f, 0.0f));
		normals.push_back(vec3(-1.0f, 0.0f, 0.0f));
		normals.push_back(vec3(-1.0f, 0.0f, 0.0f));
		normals.push_back(vec3(-1.0f, 0.0f, 0.0f));

		vertices.push_back(p5);
		vertices.push_back(p6);
		vertices.push_back(p7);
		vertices.push_back(p8);

		normals.push_back(vec3(0.0f, 0.0f, 1.0f));
		normals.push_back(vec3(0.0f, 0.0f, 1.0f));
		normals.push_back(vec3(0.0f, 0.0f, 1.0f));
		normals.push_back(vec3(0.0f, 0.0f, 1.0f));

		vertices.push_back(p4);
		vertices.push_back(p3);
		vertices.push_back(p7);
		vertices.push_back(p8);

		normals.push_back(vec3(0.0f, 1.0f, 0.0f));
		normals.push_back(vec3(0.0f, 1.0f, 0.0f));
		normals.push_back(vec3(0.0f, 1.0f, 0.0f));
		normals.push_back(vec3(0.0f, 1.0f, 0.0f));
	}

	uint nrVertices = vertices.size();
	VertexBufferObject::DATA *attrData = new VertexBufferObject::DATA[nrVertices];

	for (uint i = 0; i<nrVertices; ++i)
	{
		vec3 v = vertices[i];
		vec3 n = normals[i];

		attrData[i].vx = v.x;
		attrData[i].vy = v.y;
		attrData[i].vz = v.z;
		attrData[i].vw = 1.0f;

		attrData[i].nx = n.x;
		attrData[i].ny = n.y;
		attrData[i].nz = n.z;
		attrData[i].nw = 0.0f;

		attrData[i].cx = 1.0f;
		attrData[i].cy = 0.0f;
		attrData[i].cz = 0.0f;
		attrData[i].cw = 1.0f;

		attrData[i].tx = 0.0f;
		attrData[i].ty = 0.0f;
		attrData[i].tz = 0.0f;
		attrData[i].tw = 0.0f;
	}

	m_vboSolution = new VertexBufferObject();
	m_vboSolution->setData(attrData, GL_STATIC_DRAW, nrVertices, GL_QUADS);
	m_vboSolution->bindDefaultAttribs();

	delete[] attrData;
}

void Scene::buildVBOTarget(const vector<vec3> &positions, const vector<vec3> &sizes)
{
	delete m_vboTarget;

	vector<vec3> vertices, normals;
	for (int i = 0; i < positions.size(); ++i)
	{
		vec3 p = positions[i];
		vec3 s = sizes[i];

		vec3 p1 = p;
		vec3 p2 = p + vec3(s.x, 0.0f, 0.0f);
		vec3 p3 = p + vec3(s.x, s.y, 0.0f);
		vec3 p4 = p + vec3(0.0f, s.y, 0.0f);
		vec3 p5 = p + vec3(0.0f, 0.0f, s.z);
		vec3 p6 = p + vec3(s.x, 0.0f, s.z);
		vec3 p7 = p + s;
		vec3 p8 = p + vec3(0.0f, s.y, s.z);

		vertices.push_back(p1);
		vertices.push_back(p2);
		vertices.push_back(p3);
		vertices.push_back(p4);

		normals.push_back(vec3(0.0f, 0.0f, -1.0f));
		normals.push_back(vec3(0.0f, 0.0f, -1.0f));
		normals.push_back(vec3(0.0f, 0.0f, -1.0f));
		normals.push_back(vec3(0.0f, 0.0f, -1.0f));

		vertices.push_back(p1);
		vertices.push_back(p2);
		vertices.push_back(p6);
		vertices.push_back(p5);

		normals.push_back(vec3(0.0f, -1.0f, 0.0f));
		normals.push_back(vec3(0.0f, -1.0f, 0.0f));
		normals.push_back(vec3(0.0f, -1.0f, 0.0f));
		normals.push_back(vec3(0.0f, -1.0f, 0.0f));

		vertices.push_back(p2);
		vertices.push_back(p3);
		vertices.push_back(p7);
		vertices.push_back(p6);

		normals.push_back(vec3(1.0f, 0.0f, 0.0f));
		normals.push_back(vec3(1.0f, 0.0f, 0.0f));
		normals.push_back(vec3(1.0f, 0.0f, 0.0f));
		normals.push_back(vec3(1.0f, 0.0f, 0.0f));

		vertices.push_back(p1);
		vertices.push_back(p4);
		vertices.push_back(p8);
		vertices.push_back(p5);

		normals.push_back(vec3(-1.0f, 0.0f, 0.0f));
		normals.push_back(vec3(-1.0f, 0.0f, 0.0f));
		normals.push_back(vec3(-1.0f, 0.0f, 0.0f));
		normals.push_back(vec3(-1.0f, 0.0f, 0.0f));

		vertices.push_back(p5);
		vertices.push_back(p6);
		vertices.push_back(p7);
		vertices.push_back(p8);

		normals.push_back(vec3(0.0f, 0.0f, 1.0f));
		normals.push_back(vec3(0.0f, 0.0f, 1.0f));
		normals.push_back(vec3(0.0f, 0.0f, 1.0f));
		normals.push_back(vec3(0.0f, 0.0f, 1.0f));

		vertices.push_back(p4);
		vertices.push_back(p3);
		vertices.push_back(p7);
		vertices.push_back(p8);

		normals.push_back(vec3(0.0f, 1.0f, 0.0f));
		normals.push_back(vec3(0.0f, 1.0f, 0.0f));
		normals.push_back(vec3(0.0f, 1.0f, 0.0f));
		normals.push_back(vec3(0.0f, 1.0f, 0.0f));
	}

	uint nrVertices = vertices.size();
	VertexBufferObject::DATA *attrData = new VertexBufferObject::DATA[nrVertices];

	for (uint i = 0; i<nrVertices; ++i)
	{
		vec3 v = vertices[i];
		vec3 n = normals[i];

		attrData[i].vx = v.x;
		attrData[i].vy = v.y;
		attrData[i].vz = v.z;
		attrData[i].vw = 1.0f;

		attrData[i].nx = n.x;
		attrData[i].ny = n.y;
		attrData[i].nz = n.z;
		attrData[i].nw = 0.0f;

		attrData[i].cx = 0.0f;
		attrData[i].cy = 0.0f;
		attrData[i].cz = 1.0f;
		attrData[i].cw = 1.0f;

		attrData[i].tx = 0.0f;
		attrData[i].ty = 0.0f;
		attrData[i].tz = 0.0f;
		attrData[i].tw = 0.0f;
	}

	m_vboTarget = new VertexBufferObject();
	m_vboTarget->setData(attrData, GL_STATIC_DRAW, nrVertices, GL_QUADS);
	m_vboTarget->bindDefaultAttribs();

	delete[] attrData;
}

void Scene::buildBox(vector<vec3> &vertices, vector<vec3> &normals, vector<vec3> &verticesLines, vec3 p, float s)
{
    float t = s / 2.0f;
    vec3 a = p + vec3(-t, -t, -t);
    vec3 b = p + vec3( t, -t, -t);
    vec3 c = p + vec3( t, -t,  t);
    vec3 d = p + vec3(-t, -t,  t);
    vec3 e = p + vec3(-t,  t, -t);
    vec3 f = p + vec3( t,  t, -t);
    vec3 g = p + vec3( t,  t,  t);
    vec3 h = p + vec3(-t,  t,  t);

    vertices.push_back(a);
    vertices.push_back(b);
    vertices.push_back(c);
    vertices.push_back(d);

    vertices.push_back(b);
    vertices.push_back(f);
    vertices.push_back(g);
    vertices.push_back(c);

    vertices.push_back(a);
    vertices.push_back(e);
    vertices.push_back(h);
    vertices.push_back(d);

    vertices.push_back(e);
    vertices.push_back(f);
    vertices.push_back(g);
    vertices.push_back(h);

    vertices.push_back(d);
    vertices.push_back(c);
    vertices.push_back(g);
    vertices.push_back(h);
    
    vertices.push_back(a);
    vertices.push_back(b);
    vertices.push_back(f);
    vertices.push_back(e);


    normals.push_back(vec3(0, 0, -1));
    normals.push_back(vec3(0, 0, -1));
    normals.push_back(vec3(0, 0, -1));
    normals.push_back(vec3(0, 0, -1));

    normals.push_back(vec3(1, 0, 0));
    normals.push_back(vec3(1, 0, 0));
    normals.push_back(vec3(1, 0, 0));
    normals.push_back(vec3(1, 0, 0));

    normals.push_back(vec3(-1, 0, 0));
    normals.push_back(vec3(-1, 0, 0));
    normals.push_back(vec3(-1, 0, 0));
    normals.push_back(vec3(-1, 0, 0));

    normals.push_back(vec3(0, 0, 1));
    normals.push_back(vec3(0, 0, 1));
    normals.push_back(vec3(0, 0, 1));
    normals.push_back(vec3(0, 0, 1));

    normals.push_back(vec3(0, 1, 0));
    normals.push_back(vec3(0, 1, 0));
    normals.push_back(vec3(0, 1, 0));
    normals.push_back(vec3(0, 1, 0));
 
    normals.push_back(vec3(0, -1, 0));
    normals.push_back(vec3(0, -1, 0));
    normals.push_back(vec3(0, -1, 0));
    normals.push_back(vec3(0, -1, 0));


    verticesLines.push_back(a);
    verticesLines.push_back(b);
    verticesLines.push_back(b);
    verticesLines.push_back(c);
    verticesLines.push_back(c);
    verticesLines.push_back(d);
    verticesLines.push_back(d);
    verticesLines.push_back(a);

    verticesLines.push_back(e);
    verticesLines.push_back(f);
    verticesLines.push_back(f);
    verticesLines.push_back(g);
    verticesLines.push_back(g);
    verticesLines.push_back(h);
    verticesLines.push_back(h);
    verticesLines.push_back(e);

    verticesLines.push_back(a);
    verticesLines.push_back(d);
    verticesLines.push_back(e);
    verticesLines.push_back(h);
    verticesLines.push_back(b);
    verticesLines.push_back(c);
    verticesLines.push_back(f);
    verticesLines.push_back(g);


}

void Scene::buildVBOVoxels()
{
    delete m_vboVoxels;
    delete m_vboVoxelsLines;

    ifstream fin;
    fin.open("C:/Users/Vignesh/Downloads/test_txts/test_txts_3/d9432798bd2aa338ad5067eac75a07f7__1__.txt");
    //fin.open("C:/Users/Vignesh/Downloads/test_txts/test_txts/hello/vignesh.txt");

    int s = 32;    

    vector<vector<vector<int>>> voxels;
    voxels.resize(s);
    int term;

    if (!fin.eof())
    {
        for (int i = 0; i < s; ++i)
        {
            vector<vector<int>> data;
            data.resize(s);
            for (int j = 0; j < s; ++j)
            {
                vector<int> single_dim;
                single_dim.resize(s);
                for (int k = 0; k < s; ++k)
                {
                    fin >> term;
                    single_dim[k] = term;
                }
                data[j] = single_dim;
            }

            voxels[i] = data;
        }
    }
    fin.close();


    float domainSize = 1.0f;
    float boxSize = domainSize / s;

    vector<vec3> vertices, verticesLines, normals;
    vec3 boxCenter = vec3(boxSize / 2, boxSize / 2, boxSize / 2);
    vec3 stepSize = vec3(boxSize);

    int layer = 0;
    for (int i = 0; i < voxels.size(); ++i)
    {      
        for (int j = 0; j < s; ++j)
        {
            for (int k = 0; k < s; ++k)
            {                
                if (layer <= m_activeLayer)
                {
                    if (voxels[i][j][k])
                    {
                        //vec3 p = boxCenter + vec3(-j, i, k) * stepSize;
                        vec3 p = boxCenter + vec3(-i, j, k) * stepSize;
                        buildBox(vertices, normals, verticesLines, p, boxSize);
                    }
                }
            }
        }

        layer++;
    }

    uint nrVertices = vertices.size();
    VertexBufferObject::DATA *attrData = new VertexBufferObject::DATA[nrVertices];

    for (uint i = 0; i<nrVertices; ++i)
    {
        vec3 v = vertices[i];
        vec3 n = normals[i];

        attrData[i].vx = v.x;
        attrData[i].vy = v.y;
        attrData[i].vz = v.z;
        attrData[i].vw = 1.0f;

        attrData[i].nx = n.x;
        attrData[i].ny = n.y;
        attrData[i].nz = n.z;
        attrData[i].nw = 1.0f;

        attrData[i].cx = 1.0f;
        attrData[i].cy = 1.0f;
        attrData[i].cz = 1.0f;
        attrData[i].cw = 1.0f;

        attrData[i].tx = 0.0f;
        attrData[i].ty = 0.0f;
        attrData[i].tz = 0.0f;
        attrData[i].tw = 0.0f;
    }

    m_vboVoxels = new VertexBufferObject();
    m_vboVoxels->setData(attrData, GL_STATIC_DRAW, nrVertices, GL_QUADS);
    m_vboVoxels->bindDefaultAttribs();


    uint nrVerticesLines = vertices.size();
    VertexBufferObject::DATA *attrDataLines = new VertexBufferObject::DATA[nrVerticesLines];

    for (uint i = 0; i<nrVerticesLines; ++i)
    {
        vec3 v = verticesLines[i];

        attrDataLines[i].vx = v.x;
        attrDataLines[i].vy = v.y;
        attrDataLines[i].vz = v.z;
        attrDataLines[i].vw = 1.0f;
        
        attrDataLines[i].nx = 0.0f;
        attrDataLines[i].ny = 0.0f;
        attrDataLines[i].nz = 0.0f;
        attrDataLines[i].nw = 0.0f;
        
        attrDataLines[i].cx = 0.0f;
        attrDataLines[i].cy = 0.0f;
        attrDataLines[i].cz = 0.0f;
        attrDataLines[i].cw = 1.0f;
        
        attrDataLines[i].tx = 0.0f;
        attrDataLines[i].ty = 0.0f;
        attrDataLines[i].tz = 0.0f;
        attrDataLines[i].tw = 0.0f;
    }

    m_vboVoxelsLines = new VertexBufferObject();
    m_vboVoxelsLines->setData(attrDataLines, GL_STATIC_DRAW, nrVerticesLines, GL_LINES);
    m_vboVoxelsLines->bindDefaultAttribs();

    delete[] attrDataLines;
}

void Scene::erasePoints(float dx, float dy, float dh)
{
    m_ballEraserPos.x += dx;
    m_ballEraserPos.z += dy;
    m_ballEraserPos.y += dh;    

    for (QHash<QString, Vertex>::iterator i = m_verticesScan.begin(); i != m_verticesScan.end(); ++i)
    {
        if (length(i.value().position - m_ballEraserPos) < (m_ballEraserRadius + rand<float>(-0.05f, 0.05f)))
        {
            m_verticesScan.erase(i);
        }        
    }

    m_parent->m_renderer->m_scanner->updateBuffer();
    
    qDebug() << m_verticesScan.size();
}

vector<float> Scene::computeMeshDistance()
{

    load_object("test/lamp87_skinned.obj");
    vector<Vertex> shape2_vertices = m_object->m_vertices;
    
    load_object("C:/Users/Vignesh/Documents/shape-grammars/Data/Grammars/lamp/19fc4044912f5cc4869594a32151bfdf/model.obj");
    vector<Vertex> shape1_vertices = m_object->m_vertices;
    
    vector<float> diffs;
    float max_diff = 0.0f;
    
    if (shape1_vertices.size() == shape2_vertices.size())
    {
        diffs.resize(shape1_vertices.size());
        for (int i = 0; i < diffs.size(); ++i)
        {
            vec3 pos1 = shape1_vertices[i].position, pos2 = shape2_vertices[i].position;
            diffs[i] = length(pos1 - pos2);
            if (max_diff < diffs[i])
                max_diff = diffs[i];
        }

        for (int i = 0; i < diffs.size(); ++i)
        {
            diffs[i] /= max_diff;
            //m_object->m_vertices[i].color = (1 - diffs[i])*vec4(1.0, 0.0, 0.0, 1.0) + diffs[i] * vec4(0.0, 0.0, 1.0, 1.0);
            //shape1_vertices[i].color = (1 - diffs[i])*vec4(1.0, 0.0, 0.0, 1.0) + diffs[i] * vec4(0.0, 0.0, 1.0, 1.0);

            float c[3];
            colorMap(diffs[i], c, colorJet);
            shape1_vertices[i].color = vec4(c[0], c[1], c[2], 1.0f);

            shape1_vertices[i].normal = m_object->m_vertices[i].normal;
        }


        
    }
    else
        qDebug() << "Number of vertices don't match.";

    buildVBOModelDifference(shape1_vertices);

    return diffs;
}

void Scene::buildVBOModelDifference(vector<Vertex> &vertices)
{
    VertexBufferObject::DATA *data = new VertexBufferObject::DATA[vertices.size()];

    vector<vec3> normals;
    for (uint i = 0; i < vertices.size(); i += 3)
    {
        vec3 a = vertices[i].position;
        vec3 b = vertices[i + 1].position;
        vec3 c = vertices[i + 2].position;

        vec3 s = a - b;
        vec3 t = a - c;

        vec3 n = normalize(cross(s, t));

        normals.push_back(n);
        normals.push_back(n);
        normals.push_back(n);
    }

    for (uint i = 0; i<vertices.size(); ++i)
    {
        vec3 v = vertices[i].position;
        vec3 n = normals[i];
        vec4 c = vertices[i].color;

        data[i].vx = v.x;
        data[i].vy = v.y;
        data[i].vz = v.z;
        data[i].vw = 0.0f;

        data[i].cx = c.x;
        data[i].cy = c.y;
        data[i].cz = c.z;
        data[i].cw = c.w;

        data[i].nx = n.x;
        data[i].ny = n.y;
        data[i].nz = n.z;
        data[i].nw = 0.0f;

        data[i].tx = 0.0f;
        data[i].ty = 0.0f;
        data[i].tz = 0.0f;
        data[i].tw = 0.0f;
    }

    m_vboModelDifference = new VertexBufferObject();
    m_vboModelDifference->setData(data, GL_STATIC_DRAW, vertices.size(), GL_TRIANGLES);
    m_vboModelDifference->bindDefaultAttribs();

    delete[] data;
}

void Scene::generateContinuousSkinning()
{
    m_mean.clear();
    if (m_mean.size() == 0)
    {
        m_mean = m_grammar->setParameterVals();
        m_std.resize(m_mean.size());
        for (int i = 0; i < m_mean.size(); ++i)
            m_std[i] = 0.00*sqrt(i);
    }
    for (int i = 0; i < m_std.size(); ++i)
        qDebug() << m_mean[i] << m_std[i];

    Proc_Model p;
    QMultiMap<QString, Node*> nodes = m_grammar->nodes();

    vector<vec3> positions, sizes, boxEnds;
    vector<int> colors;
    for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
    {
        Node *n = iter.value();
        positions.push_back(n->translation());
        sizes.push_back(n->size());
        colors.push_back(n->colorMode());
    }

    for (int i = 0; i < positions.size(); i++)
    {
        boxEnds.push_back(positions[i] - sizes[i] / 2);
        boxEnds.push_back(positions[i] + sizes[i] / 2);
        positions[i] -= sizes[i] / 2;
    }

    p.new_proc(positions, sizes);
    p.setColor(colors);

    qDebug() << "Done 1";

    //m_clusters->OrganizeClusters();
    srand(time(0));
    vector<double> new_params;
    new_params.resize(m_mean.size());
    for (int i = 0; i < m_mean.size(); ++i)
        new_params[i] = (m_mean[i] + m_std[i] * rand(0.0, 1.0));//*rand(0.0,1.0);

    vector<vector<double>> new_params_1 = load_matrix("test/test_params.txt");
    new_params = new_params_1[0];

    vector<double> solns;
    solns.resize(new_params.size());

    for (int nums = 1; nums < 101; ++nums)
    {
        for (int i = 0; i < new_params.size(); ++i)
            solns[i] = m_mean[i] + (new_params[i] - m_mean[i])*nums / 100;

        m_grammar->softCleanUp();
        //m_grammar->adjustParameters(m_mean.data(), m_grammar->determineNrIndependentVariables());
        m_grammar->adjustParameters(solns.data(), m_grammar->determineNrIndependentVariables());
        m_grammar->derive();

        Proc_Model q;
        nodes = m_grammar->nodes();

        positions.clear();
        sizes.clear();
        boxEnds.clear();
        for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
        {
            Node *n = iter.value();
            positions.push_back(n->translation());
            sizes.push_back(n->size());
        }

        for (int i = 0; i < positions.size(); i++)
        {
            boxEnds.push_back(positions[i] - sizes[i] / 2);
            boxEnds.push_back(positions[i] + sizes[i] / 2);
            positions[i] -= sizes[i] / 2;
        }

        q.new_proc(positions, sizes);
        q.setColor(colors);

        qDebug() << "Done 2";

        vector<Vertex> vertices = m_object->transform(p, q);
        //QString fname = "test/teaser_1_skinned.obj";
        //m_object->saveObj(fname, vertices);
        m_object->m_vbosTriangles_skin.clear();
        m_object->buildVBOMeshComplete(vertices);
        m_parent->updateGL();
        if (nums < 10)
            saveFrameBuffer(m_parent, "0" + QString::number(nums-1));
        else
            saveFrameBuffer(m_parent, QString::number(nums - 1));

        qDebug() << "Done 3";
    }
}

void Scene::generateContinuousHeatMap()
{
    m_mean.clear();
    if (m_mean.size() == 0)
    {
        m_mean = m_grammar->setParameterVals();
        m_std.resize(m_mean.size());
        for (int i = 0; i < m_mean.size(); ++i)
            m_std[i] = 0.00*sqrt(i);
    }
    for (int i = 0; i < m_std.size(); ++i)
        qDebug() << m_mean[i] << m_std[i];

    Proc_Model p;
    QMultiMap<QString, Node*> nodes = m_grammar->nodes();

    vector<vec3> positions, sizes, boxEnds;
    vector<int> colors;
    for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
    {
        Node *n = iter.value();
        positions.push_back(n->translation());
        sizes.push_back(n->size());
        colors.push_back(n->colorMode());
    }

    for (int i = 0; i < positions.size(); i++)
    {
        boxEnds.push_back(positions[i] - sizes[i] / 2);
        boxEnds.push_back(positions[i] + sizes[i] / 2);
        positions[i] -= sizes[i] / 2;
    }

    p.new_proc(positions, sizes);
    p.setColor(colors);

    qDebug() << "Done 1";

    //m_clusters->OrganizeClusters();
    srand(time(0));
    vector<double> new_params;
    new_params.resize(m_mean.size());
    for (int i = 0; i < m_mean.size(); ++i)
        new_params[i] = (m_mean[i] + m_std[i] * rand(0.0, 1.0));//*rand(0.0,1.0);

    vector<vector<double>> new_params_1 = load_matrix("test/test_params.txt");
    new_params = new_params_1[0];

    vector<double> solns;
    solns.resize(new_params.size());

    for (int nums = 1; nums < 101; ++nums)
    {
        for (int i = 0; i < new_params.size(); ++i)
            solns[i] = m_mean[i] + (new_params[i] - m_mean[i])*nums / 100;

        m_grammar->softCleanUp();
        //m_grammar->adjustParameters(m_mean.data(), m_grammar->determineNrIndependentVariables());
        m_grammar->adjustParameters(solns.data(), m_grammar->determineNrIndependentVariables());
        m_grammar->derive();

        Proc_Model q;
        nodes = m_grammar->nodes();

        positions.clear();
        sizes.clear();
        boxEnds.clear();
        for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
        {
            Node *n = iter.value();
            positions.push_back(n->translation());
            sizes.push_back(n->size());
        }

        for (int i = 0; i < positions.size(); i++)
        {
            boxEnds.push_back(positions[i] - sizes[i] / 2);
            boxEnds.push_back(positions[i] + sizes[i] / 2);
            positions[i] -= sizes[i] / 2;
        }

        q.new_proc(positions, sizes);
        q.setColor(colors);

        qDebug() << "Done 2";

        vector<Vertex> shape2_vertices = m_object->transform(p, q);
        vector<Vertex> shape1_vertices = m_object->m_vertices;
        vector<float> diffs;
        double max_diff = 0.0;

        if (shape1_vertices.size() == shape2_vertices.size())
        {
            diffs.resize(shape1_vertices.size());
            for (int i = 0; i < diffs.size(); ++i)
            {
                vec3 pos1 = shape1_vertices[i].position, pos2 = shape2_vertices[i].position;
                diffs[i] = length(pos1 - pos2);
                if (max_diff < diffs[i])
                    max_diff = diffs[i];
            }

            for (int i = 0; i < diffs.size(); ++i)
            {
                diffs[i] /= max_diff;
                //m_object->m_vertices[i].color = (1 - diffs[i])*vec4(1.0, 0.0, 0.0, 1.0) + diffs[i] * vec4(0.0, 0.0, 1.0, 1.0);
                //shape1_vertices[i].color = (1 - diffs[i])*vec4(1.0, 0.0, 0.0, 1.0) + diffs[i] * vec4(0.0, 0.0, 1.0, 1.0);

                float c[3];
                colorMap(diffs[i], c, colorJet);
                shape1_vertices[i].color = vec4(c[0], c[1], c[2], 1.0f);

                shape1_vertices[i].normal = m_object->m_vertices[i].normal;
            }
            buildVBOModelDifference(shape1_vertices);

            m_parent->updateGL();
            if (nums <= 10)
                saveFrameBuffer(m_parent, "0" + QString::number(nums-1));
            else
                saveFrameBuffer(m_parent, QString::number(nums-1));
        }
        qDebug() << "Done 3";
    }
}

void Scene::compareMinhyuk()
{
    int number = 7008;
    //QString filename = "C:/Users/Vignesh/Downloads/to_test/coseg_chairs/output/" + QString::number(number) + "/" + QString::number(number) + "_input.ply";
    QString filename = "C:/Users/Vignesh/Downloads/to_test/assembly_chairs/output/" + QString::number(number) + "/" + QString::number(number) + "_input.ply";
    vector<vector<vec3>> p = PointCloudIO::load_point_cloud_ply(filename);

    m_verticesScan.clear();
    for (int i = 0; i < p[0].size(); ++i)
    {
        Vertex v;
        v.position = p[0][i];
        v.normal = p[1][i];
        m_verticesScan.insert(QString::number(i), v);
    }
}

void Scene::saveObjboxes(QString filename)
{
    QMultiMap<QString, Node*> nodes = m_grammar->nodes();

    ofstream fout;
    fout.open(filename.toStdString());

    vector<vec3> positions, sizes, boxEnds;
    vector<int> colorModes;
    for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
    {
        Node *n = iter.value();
        positions.push_back(n->translation());
        sizes.push_back(n->size());
        colorModes.push_back(n->colorMode());
    }

    vector<vec3> vertices, colors; 
    vector<vector<int>> faces;
    for (int i = 0; i < positions.size(); ++i)
    {
        boxEnds.resize(2);
        boxEnds[0] = positions[i] - 0.5*sizes[i];
        boxEnds[1] = positions[i] + 0.5*sizes[i];
        
        //populate vertices
        vertices.push_back(boxEnds[0]);
        vertices.push_back(vec3(boxEnds[0].x, boxEnds[0].y, boxEnds[1].z));
        vertices.push_back(vec3(boxEnds[0].x, boxEnds[1].y, boxEnds[1].z));
        vertices.push_back(vec3(boxEnds[0].x, boxEnds[1].y, boxEnds[0].z));
        vertices.push_back(vec3(boxEnds[1].x, boxEnds[0].y, boxEnds[0].z));
        vertices.push_back(vec3(boxEnds[1].x, boxEnds[0].y, boxEnds[1].z));
        vertices.push_back(boxEnds[1]);
        vertices.push_back(vec3(boxEnds[1].x, boxEnds[1].y, boxEnds[0].z));

        //populate faces
        int t = 8 * i;
        vector<int> v;
        v.push_back(t + 1);
        v.push_back(t + 2);
        v.push_back(t + 3);
        v.push_back(t + 4);
        faces.push_back(v);
        v.clear();

        v.push_back(t + 5);
        v.push_back(t + 6);
        v.push_back(t + 7);
        v.push_back(t + 8);
        faces.push_back(v);
        v.clear();

        v.push_back(t + 1);
        v.push_back(t + 5);
        v.push_back(t + 8);
        v.push_back(t + 4);
        faces.push_back(v);
        v.clear();

        v.push_back(t + 2);
        v.push_back(t + 6);
        v.push_back(t + 7);
        v.push_back(t + 3);
        faces.push_back(v);
        v.clear();

        v.push_back(t + 1);
        v.push_back(t + 5);
        v.push_back(t + 6);
        v.push_back(t + 2);
        faces.push_back(v);
        v.clear();

        v.push_back(t + 4);
        v.push_back(t + 8);
        v.push_back(t + 7);
        v.push_back(t + 3);
        faces.push_back(v);
        v.clear();

        vec4 color;
        switch (colorModes[i])
        {
        case 0:
            color = vec4(1.0f, 0.6f, 0.3f, 1.0f);
            break;
        case 1:
            color = vec4(0.3f, 0.3f, 1.0f, 1.0f);
            break;
        case 2:
            color = vec4(1.0f, 1.0f, 0.3f, 1.0f);
            break;
        case 3:
            color = vec4(1.0f, 0.3f, 0.3f, 1.0f);
            break;
        case 4:
            color = vec4(0.3f, 1.0f, 0.0f, 1.0f);
            break;
        case 5:
            color = vec4(0.0f, 0.5f, 1.0f, 1.0f);
            break;
        case 6:
            color = vec4(0.5f, 0.5f, 1.0f, 1.0f);
            break;
        default:
            break;
        }

        for (int j = 0; j < 8; ++j)
            colors.push_back(vec3(color.x, color.y, color.z));
    }

    for (int i = 0; i < vertices.size(); ++i)
    {
        fout << "v " << vertices[i].x << " " << vertices[i].y << " " << vertices[i].z << " " << colors[i].x << " " << colors[i].y << " " << colors[i].z << endl;
    }

    for (int i = 0; i < faces.size(); ++i)
        fout << "f " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << " " << faces[i][3] << endl;

    fout.close();
}

void Scene::savePartialResults()
{
    int p_no;
    cout << "Partial Scan number:" << endl;
    cin >> p_no;

    QString m_folder = "C:/Users/Vignesh/Documents/shape-grammars/Data/Grammars/monitor/";
    QString m_grammar_folder = "D:/Source/ShapeGrammars/code/_main/Data/Grammars/monitor/";
    vector<vector<double>> optimals = load_matrix((m_folder + "monitors_all.txt").toStdString());
    vector<vector<double>> endpts = load_matrix((m_folder + "endpts_monitors.txt").toStdString());

    int shape_number;
    std::cout << "Enter a shape number: " << std::endl;
    std::cin >> shape_number;
    shape_number--;

    int req_gram;
    std::cout << "Enter the template you're interested in: " << std::endl;
    std::cin >> req_gram;
    req_gram--;

    m_grammar = new Grammar(m_grammar_folder + "monitor_gram" + QString::number(req_gram + 1) + ".gram");

    int start_val = endpts[req_gram][0], end_val = endpts[req_gram][1];

    if (shape_number > optimals.size() - 1)
        std::cout << "Sorry, this is beyond the available limit.";
    else
    {
        //QString file1 = m_folder + "model (" + QString::number(shape_number + 1) + ").obj";
        //
        //load_object(file1);
        
        //qDebug() << "Number of vertices are: " << m_object->m_vertices.size();

        QString file1;
        if (shape_number + 1 < 10)
            file1 = m_folder + "model000" + QString::number(shape_number + 1) + ".obj";
        else if (shape_number + 1 < 100)
            file1 = m_folder + "model00" + QString::number(shape_number + 1) + ".obj";
        else if (shape_number + 1 < 1000)
            file1 = m_folder + "model0" + QString::number(shape_number + 1) + ".obj";
        else file1 = m_folder + "model" + QString::number(shape_number + 1) + ".obj";
        load_object(file1);

        vector<double> original_params;

        for (int i = start_val; i <= end_val; ++i)
        {
            original_params.push_back(optimals[shape_number][i - 1]);
        }

        //vector<vector<double>> original_params1 = load_matrix("C:/Users/Vignesh/Downloads/objects/TestSet/desk/desk0961_fit_3.txt");
        //original_params = original_params1[0];

        //m_grammar = new Grammar(m_grammar_folder + "desk_gram4.gram");
        //vector<vector<double>> orig = load_matrix("Paper_Data/Matthias_Comparison/desk351_4.txt");
        //original_params = orig[0];

        m_grammar->softCleanUp();
        m_grammar->adjustParameters(original_params.data(), m_grammar->determineNrIndependentVariables());
        m_grammar->derive();

        Proc_Model p;
        QMultiMap<QString, Node*> nodes = m_grammar->nodes();

        vector<vec3> positions, sizes, boxEnds;
        vector<int> colors;
        for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
        {
            Node *n = iter.value();
            positions.push_back(n->translation());
            sizes.push_back(n->size());
            colors.push_back(n->colorMode());
        }

        for (int i = 0; i < positions.size(); i++)
        {
            boxEnds.push_back(positions[i] - sizes[i] / 2);
            boxEnds.push_back(positions[i] + sizes[i] / 2);
            positions[i] -= sizes[i] / 2;
        }

        p.new_proc(positions, sizes);
        p.setColor(colors);

        //QString box_filename_0 = "C:/Users/Vignesh/Downloads/to_vignesh/975/" + QString::number(p_no) + "_inp_boxes.obj";
        //saveObjboxes(box_filename_0);
        //string filename;
        //std::cout << "Input file name: ";
        //std::cin >> filename;
        //vector<vector<double>> new_params_1 = load_matrix(filename);
        //vector<vector<double>> new_params_1 = load_matrix("test/test_params.txt");
        vector<vector<double>> new_params_1 = load_matrix("C:/Users/Vignesh/Downloads/objects/TestSet/monitor/monitor_scan_0034_1.txt");
        //vector<vector<double>> new_params_1 = load_matrix("C:/Users/Vignesh/Downloads/to_vignesh/975/975_box0_table1.txt");
        vector<double> new_params = new_params_1[0];

        m_grammar->softCleanUp();
        //m_grammar->adjustParameters(m_mean.data(), m_grammar->determineNrIndependentVariables());
        m_grammar->adjustParameters(new_params.data(), m_grammar->determineNrIndependentVariables());
        m_grammar->derive();

        Proc_Model q;
        nodes = m_grammar->nodes();

        positions.clear();
        sizes.clear();
        boxEnds.clear();
        for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
        {
            Node *n = iter.value();
            positions.push_back(n->translation());
            sizes.push_back(n->size());
        }

        for (int i = 0; i < positions.size(); i++)
        {
            boxEnds.push_back(positions[i] - sizes[i] / 2);
            boxEnds.push_back(positions[i] + sizes[i] / 2);
            positions[i] -= sizes[i] / 2;
        }

        q.new_proc(positions, sizes);
        q.setColor(colors);

        vector<Vertex> vertices = m_object->transform(p, q);
        //m_object = new Object("C:/Users/Vignesh/Downloads/objects/TestSet/dresser/dresser_scan_0004.obj", true, true, true, vec3(-0.0f, 0.0f, 0.0f), vec3(1.0, 1.0, 1.0), vec4(330.0f, 0.0f, 1.0f, 0.0f), vec4(1.0f, 1.0f, 1.0f, 1.0f));
        QString fname_orig = "C:/Users/Vignesh/Downloads/objects/TestSet/Results/monitor/" + QString::number(p_no) + "/monitor" + QString::number(shape_number + 1) + ".obj";
        //QString fname_orig = "C:/Users/Vignesh/Downloads/to_vignesh/975/" + QString::number(shape_number + 1) + ".obj";
        m_object->saveObj(fname_orig);
        QString fname = "C:/Users/Vignesh/Downloads/objects/TestSet/Results/monitor/" + QString::number(p_no) + "/monitor" + QString::number(shape_number + 1) + "_scan00" + QString::number(p_no) + "_skinned.obj";
        //QString fname = "C:/Users/Vignesh/Downloads/to_vignesh/975/" + QString::number(shape_number + 1) + "_skinned_box0.obj";
        m_object->m_vertices = vertices;
        m_object->saveObj(fname, vertices);
        //load_object(fname);

        //QString box_filename = "C:/Users/Vignesh/Downloads/to_vignesh/975/" + QString::number(p_no) + "_boxes.obj";
        QString box_filename = "C:/Users/Vignesh/Downloads/objects/TestSet/Results/monitor/" + QString::number(p_no) + "/monitor_scan00" + QString::number(p_no) + "_boxes.obj";
        saveObjboxes(box_filename);
    }
}

void Scene::sceneCompletion()
{
    QString input_folder = "C:/Users/Vignesh/Downloads/to_vignesh/";
    int scene = 975;
    QString bbox_name = input_folder + "box3d_" + QString::number(scene) + ".txt";
    QString pc_name = input_folder + "pc_" + QString::number(scene) + ".txt";
    vector<vector<double>> bbox = load_matrix(bbox_name.toStdString());
    vector<vector<double>> pc = load_matrix(pc_name.toStdString());
    cout << "Size of bbox: " << bbox.size() << endl;
    cout << "Size of pc: " << pc.size() << endl;

    vector<vector<vec3>> bbox_vals;
    bbox_vals.resize(bbox.size());
    for (int i = 0; i < bbox.size(); ++i)
    {
        bbox_vals[i].resize(8);
        for (int j = 0; j < 8; ++j)
        {
            bbox_vals[i][j] = vec3(bbox[i][3 * j], bbox[i][3 * j + 1], bbox[i][3 * j + 2]);
        }
    }
    
    m_verticesScan.clear();
    int bbox_num = 4;
    vector<vec3> centers;
    centers.resize(bbox_vals.size());
    for (int i = 0; i < centers.size(); ++i)
    {
        for (int j = 0; j < 8; ++j)
            centers[i] = centers[i] + bbox_vals[i][j];

        centers[i] = centers[i] / 8.0f;
    }

    vector<vec3> bbox_act;
    bbox_act.push_back(bbox_vals[bbox_num][0]);
    bbox_act.push_back(bbox_vals[bbox_num][0]);
    for (int j = 0; j < 8; ++j)
    {
        bbox_act[0] = minimum(bbox_act[0], bbox_vals[bbox_num][j]);
        bbox_act[1] = maximum(bbox_act[1], bbox_vals[bbox_num][j]);
    }
    cout << bbox_act[0].x << " " << bbox_act[0].y << " " << bbox_act[0].z << endl;
    cout << bbox_act[1].x << " " << bbox_act[1].y << " " << bbox_act[1].z << endl;

    double angle = 85.0;
    double theta = angle*(6.28f / 360.0f);

    vector<Vertex> pc_vec;
    pc_vec.resize(pc.size());
    for (int i = 0; i < pc.size(); ++i)
    {
        Vertex v;
        v.position = vec3(pc[i][0], pc[i][1], pc[i][2]);
        v.normal = vec3(0.0f, 0.0f, 0.0f);
        //m_verticesScan.insert(QString::number(i), v);
        //pc_vec[i] = v;
        
        /*if (!((v.position.x >= bbox_act[0].x) && (v.position.x <= bbox_act[1].x) && (v.position.y >= bbox_act[0].y) && (v.position.y <= bbox_act[1].y) && (v.position.z >= bbox_act[0].z) && (v.position.z <= bbox_act[1].z)))
        {
            m_verticesScan.insert(QString::number(i), v);
            pc_vec[i] = v;
        }*/

        
        /*if ((v.position.x >= bbox_act[0].x) && (v.position.x <= bbox_act[1].x) && (v.position.y >= bbox_act[0].y) && (v.position.y <= bbox_act[1].y) && (v.position.z >= bbox_act[0].z) && (v.position.z <= bbox_act[1].z))
        {
            v.position = vec3(0.0f, centers[bbox_num].y - 0.2f, 0.0f) - (v.position - centers[bbox_num]);
            float sx = v.position.x, sz = v.position.z;
            v.position.x = sx*cos(theta) + sz*sin(theta);
            v.position.z = sz*cos(theta) - sx*sin(theta);
           m_verticesScan.insert(QString::number(i), v);
        }*/
        //cout << i << endl;
    }
    //m_object->saveObjPC(input_folder + QString::number(scene) + ".obj", pc_vec);

    m_parent->m_renderer->m_scanner->updateBuffer();

    int saveObject = 0;
    cout << "Do you want to save the object?" << endl;
    cin >> saveObject;

    if (saveObject)
    {
        load_object(input_folder + QString::number(scene) + "/1064_skinned_box" + QString::number(bbox_num) + ".obj", false);
        vector<Vertex> vertices = m_object->m_vertices;
        for (int i = 0; i < vertices.size(); ++i)
        {
            Vertex v = vertices[i];
            float sx = v.position.x, sz = v.position.z;
            v.position.x = sx*cos(theta) - sz*sin(theta);
            v.position.z = sx*sin(theta) + sz*cos(theta);
            //cout << v.position.y << endl;
            v.position = vec3(0.0f, centers[bbox_num].y - 0.2f, 0.0f) - (v.position - centers[bbox_num]);
            vertices[i] = v;
        }
        QString fname = input_folder + QString::number(scene) + "/1064_skinned_box" + QString::number(bbox_num) + "_rescaled.obj";
        m_object->m_vertices = vertices;
        m_object->saveObj(fname, vertices);
    }

    //Box-removal code
    pc_vec.clear();
    pc_vec.resize(pc.size());
    bbox_act.clear();
    bbox_act.resize(2*pc.size());
    int num_boxes = 4;
    for (int i = 0; i < num_boxes; ++i)
    {
        bbox_act.push_back(bbox_vals[i][0]);
        bbox_act.push_back(bbox_vals[i][0]);
        for (int j = 0; j < 8; ++j)
        {
            bbox_act[2*i] = minimum(bbox_act[2*i], bbox_vals[i][j]);
            bbox_act[2*i + 1] = maximum(bbox_act[2*i + 1], bbox_vals[i][j]);
        }
    }
    for (int i = 0; i < pc.size(); ++i)
    {
        Vertex v;
        v.position = vec3(pc[i][0], pc[i][1], pc[i][2]);
        v.normal = vec3(0.0f, 0.0f, 0.0f);
        
        bool k = true;
        for (int j = 0; j < num_boxes; ++j)
        {
            k = (k&(!((v.position.x >= bbox_act[2 * j].x) && (v.position.x <= bbox_act[2 * j + 1].x) && (v.position.y >= bbox_act[2 * j].y) && (v.position.y <= bbox_act[2 * j + 1].y) && (v.position.z >= bbox_act[2 * j].z) && (v.position.z <= bbox_act[2 * j + 1].z))));
        }
        if (k)
        {
            m_verticesScan.insert(QString::number(i), v);
            pc_vec[i] = v;
        }
    }
    m_object->saveObjPC(input_folder + QString::number(scene) + "_box_removed.obj", pc_vec);

}

vec3 Scene::minimum(vec3 X, vec3 Y)
{
    return vec3(min(X.x,Y.x), min(X.y,Y.y), min(X.z,Y.z));
}

vec3 Scene::maximum(vec3 X, vec3 Y)
{
    return vec3(max(X.x, Y.x), max(X.y, Y.y), max(X.z, Y.z));
}

double Scene::norm(vector<double> x, vector<double> y)
{
    double n = 0.0;
    for (int i = 0; i < x.size(); ++i)
        n += (x[i] - y[i])*(x[i] - y[i]);

    return sqrt(n);
}

void Scene::findClosestShape()
{
    QString m_folder = "C:/Users/Vignesh/Documents/shape-grammars/Data/Grammars/monitor/";
    QString m_grammar_folder = "D:/Source/ShapeGrammars/code/_main/Data/Grammars/monitor/";
    vector<vector<double>> optimals = load_matrix((m_folder + "monitors_all.txt").toStdString());
    vector<vector<double>> endpts = load_matrix((m_folder + "endpts_monitors.txt").toStdString());

    int req_gram;
    std::cout << "Enter the template you're interested in: " << std::endl;
    std::cin >> req_gram;
    req_gram--;

    int start_val = endpts[req_gram][0], end_val = endpts[req_gram][1];
    vector<vector<double>> new_params_1 = load_matrix("C:/Users/Vignesh/Downloads/objects/TestSet/monitor/monitor_scan_0034_1.txt");
    vector<double> new_params = new_params_1[0];
    double norm_diff = 0.0;
    int corr_shape = 0;
    for (int shape_number = 0; shape_number < optimals.size(); ++shape_number)
    {
        vector<double> original_params;
        for (int i = start_val; i <= end_val; ++i)
        {
            original_params.push_back(optimals[shape_number][i - 1]);
        }

        if ((shape_number == 0) || (norm(original_params, new_params) < norm_diff))
        {
            norm_diff = norm(original_params, new_params);
            corr_shape = shape_number;
        }
    }

    cout << corr_shape+1 << endl;
}