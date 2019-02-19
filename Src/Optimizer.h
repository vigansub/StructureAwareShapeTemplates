#ifndef OPTIMIZER
#define OPTIMIZER

#include "Headers.h"
#include "desolver.h"
#include "Proc_Model.h"

class Grammar;

class Optimizer : public DESolver
{

public:
    Optimizer(Scene *scene, double scale, double probability, int maxGenerations, int nrDimensions);
    ~Optimizer();

    double error() const { return energy(); }
    const double* solve();

    void setRange(const vector<double> &mi, const vector<double> &ma);
    void setDistanceField(const vector<double> &df);
	void setSampleVals(const vector<vec3> &samples);
	void setGrammar(const vector<double> &grammar);

    double errorBB(Proc_Model &p, vector<vec3> &bbox, vector<vec3> &boxEnds);
	double errorGrammar(vector<vector<int>> &v, vector<vec3> &boxEnds);
	double errorDisentangle(Proc_Model &p);
    double errorSkinning(Proc_Model &p);

	vec3 conventionPoint(vec3 &start, vec3 &end, int num);

	void setGrammar(Grammar *gram);
      
private:
    double energy_function(double trial[], bool &foundSolution);

	double energy_function_proc_model(double trial[], bool &foundSolution);
	double energy_function_grammar(double trial[], bool &foundSolution);

    vector<double> m_minRange;
    vector<double> m_maxRange;
    vector<double> m_df;
	vector<vec3> m_samples;

	//vector<double> m_grammar;

    double m_scale;
    double m_probability;
    int m_maxGenerations;
	int m_nrDimensions;

	Grammar *m_grammar;
};


#endif
