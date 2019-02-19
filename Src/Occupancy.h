#ifndef OCUUPANCY_H
#define OCUUPANCY_H

#include "Headers.h"
#include "vector"
#include "Proc_Model.h"

class Occupancy
{
public:
    Occupancy(vector<vec3> points);
    ~Occupancy();

    vector<vector<vector<int>>> return_occupancy();
    vector<vec3> bbox(vector<vec3> points);
    void create_occupancy();
    void write_grid();

private:
    vector<vector<vector<int>>> m_occupancy_grid;
    int m_numx, m_numy, m_numz;
    vector<vec3> m_bbox;
    vector<vec3> m_pc;
};
#endif