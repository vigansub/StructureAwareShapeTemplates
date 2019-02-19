#ifndef NAIVE
#define NAIVE

#include <iostream>
#include <random>
#include <math.h>
#include "Headers.h"
#include "Proc_Model.h"
#include <time.h>
#include "Grammar.h"

class Scene;

class Genetic
{
private:
    int m_iter, m_dimensions, m_numbits, m_population;
    vector<double> m_mi, m_ma, m_fitness, m_best;
    vector<vector<double>> m_chromosomes;
    double m_crossover, m_mutation;
    
    vector<double> m_parameters;

    vector<vec3> m_samples;

    Grammar *m_grammar;
    Scene *m_scene;

public:
    Genetic(Scene *scene, int iterations, int bits, int population, vector<double> min, vector<double> max, double crossover, double mutation);
    ~Genetic();

    double energy_function(vector<double> position);

    void optimize();

    int selection();
    void swap(vector<double> &a, vector<double> &b, vector<int> p);
    void mutate(vector<double> &a, double mutation);
    int find_best_fitness();

    double unif_rand();
    vector<double> norm_fitness(vector<double> &fitness);
    vector<double> bin2dec(vector<vector<double>> &a);
    vector<vector<double>> dec2bin(vector<double> &a);
    
    void setSamples(vector<vec3> samplevals);
    void setGrammar(Grammar *grammar);

    double errorBB(vector<vec3> &bbox, Proc_Model &p, vector<vec3> &boxEnds);
    double errorDisentangle(Proc_Model &p);
    double errorSkinning(Proc_Model &p);

    vector<double> global_opt();
};

#endif