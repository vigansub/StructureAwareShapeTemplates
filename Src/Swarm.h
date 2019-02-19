#ifndef SWARM
#define SWARM

#include <iostream>
#include <random>
#include <math.h>
#include "Headers.h"
#include "Proc_Model.h"
#include "Grammar.h"
#include <time.h>

class Scene;

class Swarm
{
private:
    int m_particles, m_iter, m_dimensions;
    vector<double> m_global, m_mi, m_ma;
    vector<vector<double>> m_positions;
    vector<vector<double>> m_locals;
    vector<vector<double>> m_velocities;

    vector<double> m_parameters;
    
    vector<vec3> m_samples;

    Grammar *m_grammar;
    Scene *m_scene;

public:
    Swarm(Scene *scene, Grammar* grammar, int particles, int iterations, vector<double> min, vector<double> max, double w, double phi_p, double phi_g);
    ~Swarm();

    double energy_function(vector<double> position);

    void optimize();

    void setSamples(vector<vec3> samplevals);
    void setGrammar(Grammar *grammar);

    vector<double> global_opt();

    double errorBB(vector<vec3> &bbox, Proc_Model &p, vector<vec3> &boxEnds);
    double errorDisentangle(Proc_Model &p);
    double errorSkinning(Proc_Model &p);
};

#endif