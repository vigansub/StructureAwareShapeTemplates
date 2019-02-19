#ifndef	SIMULATEDANNEALING_H
#define SIMULATEDANNEALING_H

#include <iostream>
#include <random>
#include <math.h>
#include "Headers.h"
#include "Proc_Model.h"
#include "Grammar.h"

class Scene;

class Sim_Anneal
{
private:
    vector<double> m_trial;
    vector<double> m_old;
    vector<vec3> m_samples;
    int m_dimension;

    Grammar *m_grammar;
    Scene *m_scene;

public:
	Sim_Anneal(int dimension, Scene *scene);
	~Sim_Anneal();

    vector<double> optimize();
    double acceptance_probability(double old_cost, double new_cost, double T);
    double energy_function();

    double errorBB(Proc_Model &p, vector<vec3> &boxEnds);
    double errorDisentangle(Proc_Model &p);
    double errorSkinning(Proc_Model &p);

    void setSampleVals(const vector<vec3> &samples);
    void setGrammar(Grammar *grammar);
    void generateTrial();
	/*double acceptance_probability(double old_cost, double new_cost, double T);
	vector<double> sim_ann_normal_dist(vector<vector<double>> samples, gaussParams g);
	vector<double> sim_ann_normal_dist_var(vector<vector<double>> samples, gaussParams g);
	vector<double> sim_ann_dist_field(MH_MC m, const float* data);
	vector<double> sim_ann_dist_field_self(MH_MC m, Proc_Model p);
	vector<double> optimization(MH_MC m, Proc_Model proc);
	vector<double> samples_gen(vector<double> mean, Eigen::MatrixXd variance);
	vector<double> samples_gen_var(vector<double> inp);*/
};

#endif
