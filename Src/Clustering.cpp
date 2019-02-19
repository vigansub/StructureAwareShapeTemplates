#include "Clustering.h"
#include <time.h>
#include "Scene.h"
#include "Node.h"
#include "CMAESWrapper.h"

#include "QDir"
#include "QDirIterator"


Clusters::Clusters(Scene *scene) : m_object(nullptr)
{
    m_scene = scene;
    m_folder = "../test/test/";

    QDirIterator gt(m_folder, QStringList() << "*.gram");
    while (gt.hasNext())
    {
        gt.next();
        m_grammarlist.push_back(m_folder + gt.fileName());
    }
    m_totalgrams = m_grammarlist.size();
}

Clusters::~Clusters()
{

}

vector<int> Clusters::obtain_clusters()
{
    return m_clusters;
}


void Clusters::load_error()
{
    m_errormatrix.clear();
    QString path = m_folder + "ErrorMatrix.txt";
    ifstream fin;

    int total_grams = 0;
    
    fin.open(path.toStdString(), std::ifstream::in);

    double error;
    fin >> error;
    while (!fin.eof())
    {
        vector<double> error_shape;
        error_shape.resize(m_totalgrams);
        error_shape[0] = error;
        for (int i = 1; i < m_totalgrams; ++i)
        {
            fin >> error;
            error_shape[i] = error;
        }
        fin >> error;
        m_errormatrix.push_back(error_shape);
    }
}

void Clusters::load_optimal()
{
    m_optimal.clear();

    QString path = m_folder + "OptimalVectors.txt";
    ifstream fin;
    fin.open(path.toStdString(), std::ifstream::in);
    double parameter;
    fin >> parameter;
    while (!fin.eof())
    {
        QDirIterator gt(m_folder, QStringList() << "*.gram");
        while (gt.hasNext())
        {
            gt.next();
            m_grammar = new Grammar(m_folder + gt.fileName());
            int nr = m_grammar->determineNrIndependentVariables();
            vector<double> parameter_vals;
            parameter_vals.resize(nr);
            parameter_vals[0] = parameter;
            for (int i = 1; i < nr; ++i)
            {
                fin >> parameter;
                parameter_vals[i] = parameter;
            }
            fin >> parameter;
            m_optimal.push_back(parameter_vals);
        }
    }
}

void Clusters::load_boxes()
{
    m_boxes.resize(m_optimal.size());
    QDirIterator gt(m_folder, QStringList() << "*.gram");
    for (int i = 0; i < m_optimal.size(); ++i)
    {
        m_grammar = new Grammar(m_grammarlist[i%m_totalgrams]);
        m_grammar->softCleanUp();
        m_grammar->adjustParameters(m_optimal[i].data(), m_grammar->m_nrIndependentVariables);
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

        m_boxes[i].new_proc(positions, sizes);
    }
}

void Clusters::load_labels()
{
    m_grammar_labels.resize(m_optimal.size()/m_totalgrams);
    QString path = m_folder + "Classification.txt";
    ifstream fin;
    fin.open(path.toStdString(), std::ifstream::in);
    int label;
    int i = 0;
    while (!fin.eof())
    {
        fin >> label;
        m_grammar_labels[i] = label;
        i++;
    }
    fin.close();
}

vector<vector<double>> Clusters::error()
{
    return m_errormatrix;
}

vector<vector<vector<double>>> Clusters::create_training_input()
{
    vector<vector<vector<double>>> train_input;
    train_input.resize(m_totalgrams);
    for (int i = 0; i < m_totalgrams; ++i)
    {
        train_input[i].resize(m_numtraining);
        for (int j = 0; j < m_numtraining; ++j)
        {
            train_input[i][j] = m_optimal[i + j*m_totalgrams];
        }
    }

    return train_input;
}

vector<vector<int>> Clusters::create_training_output()
{
    vector<vector<int>> train_output;
    train_output.resize(m_totalgrams);
    for (int i = 0; i < m_totalgrams; ++i)
    {
        train_output[i].resize(m_numtraining); 
        for (int j = 0; j < m_numtraining; ++j)
            if (m_grammar_labels[j] == i+1)
                train_output[i][j] = 1;
            else train_output[i][j] = 0;
    }

    return train_output;
}

vector<double> Clusters::binary_logistic_classifier_train(vector<vector<double>> train_input, vector<int> train_output)
{
    int train_size = train_output.size();
    double alpha = 0.1;
    int data_dimension = train_input[0].size();
    vector<double> train_parameters;
    
    train_parameters.resize(data_dimension);
    
    for (int i = 0; i < data_dimension; ++i)
        train_parameters[i] = rand(0.0, 0.1);

    double error = 100.0;

    while (error > 0.01)
    {
        vector<double> theta = train_parameters;
        for (int i = 0; i < data_dimension; ++i)
        {
            for (int j = 0; j < train_size; ++j)
            {
                train_parameters[i] += alpha*(train_output[j] - exp(log_sigmoid(train_input[j], theta)))*train_input[j][i];
            }
        }

        error = 0.0;

        for (int i = 0; i < data_dimension; ++i)
            error += (train_parameters[i] - theta[i])*(train_parameters[i] - theta[i]);
    }

    return train_parameters;
}

vector<vector<double>> Clusters::classify_logistic_grammars(vector<vector<vector<double>>> train_input, vector<vector<int>> train_output)
{
    vector<vector<double>> train_parameters;
    train_parameters.resize(m_totalgrams);
    for (int i = 0; i < m_totalgrams; ++i)
        train_parameters[i] = binary_logistic_classifier_train(train_input[i], train_output[i]);

    return train_parameters;
}

vector<vector<double>> Clusters::create_LR_coeffs()
{
    vector<vector<vector<double>>> train_input = create_training_input();
    vector<vector<int>> train_output = create_training_output();
    vector<vector<double>> parameters = classify_logistic_grammars(train_input, train_output);

    vector<vector<double>> LR_coeffs;
    vector<double> temp;
    int total_shapes = m_optimal.size() / m_totalgrams;

    LR_coeffs.resize(m_totalgrams);
    for (int i = 0; i < m_totalgrams; ++i)
    {
        LR_coeffs[i].resize(total_shapes);
        for (int j = 0; j < total_shapes; ++j)
        {
            LR_coeffs[i][j] = exp(log_sigmoid(m_optimal[j*m_totalgrams + i], parameters[i]));
        }
    }

    for (int i = 0; i < m_grammar_labels.size(); ++i)
        temp.push_back(LR_coeffs[m_grammar_labels[i] - 1][i]);

    return LR_coeffs;
}

void Clusters::kmeans(vector<vector<double>> data, int transpose_req)
{
    int total_shapes = m_optimal.size() / m_totalgrams;
    m_clusters.resize(total_shapes);
    vector<vector<double>> cluster_centers;
    cluster_centers.resize(m_numclusters);
    if (transpose_req)
        data = transpose(data);
    for (int i = 0; i < m_clusters.size(); ++i)
    {
        m_clusters[i] = rand() % m_numclusters;
    }

    for (int i = 0; i < m_numclusters; ++i)
        cluster_centers[i].resize(data[0].size());

    bool converged = false;
    while (!converged)
    {
        vector<int> old_clusters = m_clusters;
        vector<vector<vector<double>>> spec_cluster;
        spec_cluster.resize(m_numclusters);
        for (int i = 0; i < data.size(); ++i)
        {
            spec_cluster[m_clusters[i]].push_back(data[i]);
        }
        
        for (int j = 0; j < m_numclusters; ++j)
            cluster_centers[j] = mean(spec_cluster[j]);

        findClusters(data, cluster_centers);

        if (isequal(m_clusters, old_clusters))
            converged = true;
    }

    ofstream fout;
    QString path = m_folder + "Clusters.txt";
    // TODO check if the file exist
    fout.open(path.toStdString());
    for (int i = 0; i < m_clusters.size(); ++i)
        fout << m_clusters[i] << endl;

    fout.close();

    vector<vector<int>> test_matrix;
    test_matrix.resize(m_numclusters);
    for (int i = 0; i < m_numclusters; ++i)
    {
        test_matrix[i].resize(m_totalgrams);
        for (int j = 0; j < m_totalgrams; ++j)
            test_matrix[i][j] = 0;
    }
    for (int i = 0; i < m_grammar_labels.size(); ++i)
    {
        test_matrix[m_clusters[i]][m_grammar_labels[i] - 1] ++;
    }
    
    path = m_folder + "TestClusters.txt";
    fout.open(path.toStdString());
    for (int i = 0; i < m_numclusters; ++i)
    {
        for (int j = 0; j < m_totalgrams; ++j)
            fout << test_matrix[i][j] << " ";
    
        fout << endl;
    }
    
    fout.close();

}

void Clusters::findClusters(vector<vector<double>> data, vector<vector<double>> cluster_centers)
{
    for (int i = 0; i < data.size(); ++i)
    {
        double min_dist_i = distance(data[i], cluster_centers[0]);
        m_clusters[i] = 0;
        for (int j = 0; j < cluster_centers.size(); ++j)
        {
            double d = distance(data[i], cluster_centers[j]);
            if (d < min_dist_i)
            {
                min_dist_i = d;
                m_clusters[i] = j;
            }
        }
    }
}

void Clusters::simpleClustering()
{
    vector<vector<int>> test_matrix;
    QString path_g = m_folder + "NewClassification.txt";
    ofstream gout;
    gout.open(path_g.toStdString());
    test_matrix.resize(m_totalgrams);
    for (int i = 0; i < m_totalgrams; ++i)
    {
        test_matrix[i].resize(m_totalgrams);
        for (int j = 0; j < m_totalgrams; ++j)
            test_matrix[i][j] = 0;
    }

    float thresh = 10.0f;

    for (int i = 0; i < m_errormatrix.size(); ++i)
    {
        int bestIndex = 0;
        for (int j = 0; j < m_totalgrams; ++j)
        {
            //if (m_errormatrix[i][j] < m_errormatrix[i][bestIndex])
            //    bestIndex = j;

            if (abs(m_errormatrix[i][j] - m_errormatrix[i][bestIndex]) < thresh)
            {
                int kj = m_boxes[i*m_totalgrams + j].position().size();
                int kb = m_boxes[i*m_totalgrams + bestIndex].position().size();
                if (kj < kb)
                    bestIndex = j;
            }
        }
        //test_matrix[bestIndex][m_grammar_labels[i] - 1] ++;
        gout << bestIndex << endl;
    }
    gout.close();

    ofstream fout;
    QString path = m_folder + "SimpleClusters.txt";
    // TODO check if the file exist
    fout.open(path.toStdString());
    for (int i = 0; i < m_totalgrams; ++i)
    {
        for (int j = 0; j < m_totalgrams; ++j)
            fout << test_matrix[i][j] << " ";

        fout << endl; 
    }
    fout.close();
}

void Clusters::OrganizeClusters()
{
    QString path = m_folder + "NewClassification.txt";
    ifstream fin;
    fin.open(path.toStdString(), std::ifstream::in);
    int label;
    int i = 0;
    vector<int> new_clusters;
    new_clusters.resize(m_errormatrix.size());
    while (!fin.eof() && i < new_clusters.size())
    {
        fin >> label;
        new_clusters[i] = label;
        i++;
    }
    fin.close();

    QDirIterator it(m_folder, QStringList() << "*.obj");
    int k = 0;

    while (it.hasNext())
    {
        qDebug() << it.next();

        delete m_object;
        m_object = new Object(m_folder+it.fileName(), true, true, true, vec3(0.0f, -0.5f, 0.0f), vec3(1.0, 1.0, 1.0), vec4(180.0f, 0.0f, 1.0f, 0.0f), vec4(1.0f, 1.0f, 1.0f, 1.0f));

        m_scene->saveScene(m_folder + "SimpleClusters/" + QString::number(new_clusters[k]) + "_" + it.fileName());
        k++;
    }
}

double Clusters::log_sigmoid(vector<double> x, vector<double> theta)
{
    double sum = 0.0;
    for (int i = 0; i < x.size(); ++i)
        sum += x[i] * theta[i];

    return -log(1 + exp(-sum));
}

vector<vector<double>> Clusters::transpose(vector<vector<double>> x)
{
    vector<vector<double>> y;
    y.resize(x[0].size());
    for (int i = 0; i < x[0].size(); ++i)
        y[i].resize(x.size());

    for (int i = 0; i < y.size(); ++i)
        for (int j = 0; j < y[i].size(); ++j)
            y[i][j] = x[j][i];

    return y;
}

vector<double> Clusters::mean(vector<vector<double>> a)
{
    vector<double> c;
    c.resize(a[0].size());
    c = a[0];
    for (int i = 1; i < a.size(); ++i)
        c = sum(c, a[i]);

    if (a.size() > 0)
        c = scale(c, 1.0 / a.size());

    return c;
}

double Clusters::distance(vector<double> a, vector<double> b)
{
    vector<double> c = difference(a, b);
    return sqrt(dot_product(c,c));
}

vector<double> Clusters::sum(vector<double> a, vector<double> b)
{
    vector<double> c;
    c.resize(a.size());

    for (int i = 0; i < a.size(); ++i)
        c[i] = a[i] + b[i];

    return c;
}

vector<double> Clusters::difference(vector<double> a, vector<double> b)
{
    vector<double> c = sum(a, scale(b, -1.0));
    return c;
}

double Clusters::dot_product(vector<double> a, vector<double> b)
{
    double dot = 0.0;
    for (int i = 0; i < a.size(); ++i)
        dot += a[i] * b[i];

    return dot;
}

vector<double> Clusters::pointwise_product(vector<double> a, vector<double> b)
{
    vector<double> c;
    c.resize(a.size());
    for (int i = 0; i < a.size(); ++i)
        c[i] = a[i] * b[i];

    return c;
}

vector<double> Clusters::scale(vector<double> a, double b)
{
    vector<double> c;
    c.resize(a.size());
    for (int i = 0; i < a.size(); ++i)
        c[i] = a[i] * b;

    return c;
}

bool Clusters::isequal(vector<int> a, vector<int> b)
{
    int i = 0;
    while (i < a.size())
    {
        if (a[i] != b[i])
            return false;

        i++;
    }

    return true;
}

void Clusters::render(const Transform &trans)
{
    if (m_object)
        m_object->render(trans);
}