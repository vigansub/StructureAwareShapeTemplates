#ifndef CLUSTERING
#define CLUSTERING

#include <iostream>
#include <random>
#include <math.h>
#include "Headers.h"
#include "Proc_Model.h"
#include "Grammar.h"
#include "Object.h"

class Scene;

class Clusters
{
private :
    vector<int> m_clusters;
    vector<vector<double>> m_errormatrix;
    vector<vector<double>> m_optimal;
    vector<Proc_Model> m_boxes;
    vector<QString> m_grammarlist;
    vector<int> m_grammar_labels;

    int m_totalgrams;
    int m_numtraining = 50;
    int m_numclusters = 4;

    Scene *m_scene;
    QString m_folder;
    Grammar *m_grammar;

    Object *m_object;

public:
    Clusters(Scene *scene);
    ~Clusters();

    vector<int> obtain_clusters();

    void load_error();
    void load_optimal();
    void load_boxes();
    void load_labels();
    vector<vector<double>> error();

    vector<vector<int>> create_training_output();
    vector<vector<vector<double>>> create_training_input();
    vector<double> binary_logistic_classifier_train(vector<vector<double>> train_input, vector<int> train_output);
    vector<vector<double>> classify_logistic_grammars(vector<vector<vector<double>>> train_input, vector<vector<int>> train_output);
    vector<vector<double>> create_LR_coeffs();

    void kmeans(vector<vector<double>> data, int transpose_req);
    void findClusters(vector<vector<double>> data, vector<vector<double>> cluster_centers);
    void simpleClustering();

    void render(const Transform &trans);


    void OrganizeClusters();

    double log_sigmoid(vector<double> x, vector<double> theta);
    vector<vector<double>> transpose(vector<vector<double>> x);
    vector<double> mean(vector<vector<double>> a);
    double distance(vector<double> a, vector<double> b);
    vector<double> sum(vector<double> a, vector<double> b);
    vector<double> difference(vector<double> a, vector<double> b);
    double dot_product(vector<double> a, vector<double> b);
    vector<double> pointwise_product(vector<double> a, vector<double> b);
    vector<double> scale(vector<double> a, double b);
    bool isequal(vector<int> a, vector<int> b);
};

#endif