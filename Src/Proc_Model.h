#ifndef PROC_H
#define PROC_H

#include <iostream>
#include <random>
#include <math.h>
#include "Headers.h"
#include <time.h>
#include <flann/flann.hpp>

class Proc_Model
{
private:
	vector<vec3> m_pos;
	vector<vec3> m_size;
	vector<double> m_df;
    vector<int> m_colorMode;

	vector<vector<double>> m_grammar;

	double m_range;
	int m_tot, m_numsamples_per_face;
public:
	Proc_Model();
	~Proc_Model();
	void new_proc(int num_boxes);
	void new_proc(const vector<vec3> &positions, const vector<vec3> &sizes);
    void setColor(vector<int> colors) { m_colorMode = colors; };
	vector<double> create_distance_field(vector<vec3> size_vec, vector<vec3> pos_vec);
	void create_box_distance_field(const vector<vec3> &size_vec, const vector<vec3> &pos_vec, vector<double> &df);
	vector<double> create_point_cloud_distance_field(vector<vec3> samples);
	void create_default_distance_field();
	void create_FLANN_default_distance_field();
	void sample_shape_surface(int num_points, float samples[]);
	int optimization(vector<double> solvals);

	void set_grammar(vector<vector<double>> grammar);
	
	double Df_Cube(double x1, double x2, double y1, double y2, double z1, double z2);
	double trilinear_interp(double x, double y, double z, int x1, int x2, int y1, int y2, int z1, int z2);
	double evaluate_df_val(vector<vec3> instance, vector<vec3> size);
	double compare_df(vector<vec3> instance, vector<vec3> size);
    double project_point_to_box(vec3 point, vec3 position, vec3 size);
	double project_point_to_model(vec3 point);
    double project_model_to_pc(vector<vec3> points);
   
    double volume();

    vector<vec3> bounding_box();
    vector<vec3> bounding_box(vector<Vertex> &samples);
	vector<vec3> bounding_box(vector<vec3> &samples);
	vector<vec3> &position();
	vector<vec3> &size();
	vector<double> &df();
	vector<vector<double>> &grammar();
    vector<int> &colorMode() { return m_colorMode; };

};

#endif