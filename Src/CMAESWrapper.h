#ifndef CMAESWRAPPER_H
#define CMAESWRAPPER_H

#ifdef WIN32
#include "cma-es\cmaes.h"
#else
#include "cma-es/cmaes.h"
#endif
#include "Headers.h"
#include "Proc_Model.h"
#include "Grammar.h"
#include "flann-1.7.1\include\flann\flann.h"

class Scene;

//#include <stdlib.h>
//#include <iostream>

class CMAESWrapper
{
public:

    CMAESWrapper(int dimension, Scene *scene);
    ~CMAESWrapper();

    double fitfun(double const *x);
    //static double fitfun(vector<double> x, int N, vector<double> dim);
    bool is_feasible(double const *v);

    //const double* optimize();
    vector<double> optimize();

    void setSampleVals(const vector<vec3> &samples);
    void setGrammar(Grammar *gram);
    void setVolume(float volume);
    void setGrid(vector<vector<vector<int>>> grid);
    void setGridBB(vector<vec3> m_gridbb);

    double errorBB(vector<vec3> &bbox, Proc_Model &p, vector<vec3> &boxEnds);
    double errorGrammar(vector<vector<int>> &v, vector<vec3> &boxEnds);
    double errorDisentangle(Proc_Model &p);
    double errorSkinning(Proc_Model &p);
    double errorEmpty();
    double errorPenalty(Proc_Model &p);

    vec3 conventionPoint(vec3 &start, vec3 &end, int num);

    double dimension();

private:
    vector<vec3> m_samples;
    double m_dimension;
    float m_volume;
    vector<vector<vector<int>>> m_grid;
    vector<vec3> m_gridbb;

public:
    Grammar *m_grammar;
    Scene *m_scene;

    double m_besterror;
    vector<vector<vector<double>>> all_error_vals;
    vector<vector<double>> local_error;

};

#endif

