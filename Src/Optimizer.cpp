#include "Optimizer.h"
#include <ctime>
#include "Proc_Model.h"
#include "Grammar.h"
#include "Node.h"

Optimizer::Optimizer(Scene *scene, double scale, double probability, int maxGenerations, int nrDimensions)
: DESolver(scene),
  m_scale(scale),
  m_probability(probability), 
  m_maxGenerations(maxGenerations), 
  m_nrDimensions(nrDimensions), 
  m_grammar(NULL)
{
}

Optimizer::~Optimizer() 
{
}

const double* Optimizer::solve() 
{
    const int DIMENSION = m_nrDimensions;
    const int POPULATIONSIZE = DIMENSION * 10;
    const double SCALE = m_scale;
    const double PROBABILITY = m_probability;
    const int MAXGENERATIONS = m_maxGenerations;
    const DESolver::Strategy STRATEGY = (DESolver::Strategy)params::inst()->optStrategy;

    initialize(DIMENSION, POPULATIONSIZE, m_minRange.data(), m_maxRange.data(), STRATEGY, SCALE, PROBABILITY);

    DESolver::solve(MAXGENERATIONS);  

    return solution();
}

double Optimizer::energy_function(double trial[], bool &foundSolution)
{
	double error = 0;

	error = energy_function_grammar(trial, foundSolution);

	return error;
}

double Optimizer::energy_function_proc_model(double trial[], bool &foundSolution)
{
	double error_project = 0.0, error_negative = 0.0, error_df = 0.0;
	Proc_Model p;

	vector<vec3> positions, sizes, boxEnds;
	int num_boxes = 3;

	for (int i = 0; i < num_boxes; i++)
	{
		//positions.push_back(vec3(trial[6 * i], trial[6 * i + 1], trial[6 * i + 2]));
		//sizes.push_back(vec3(trial[6 * i + 3], trial[6 * i + 4], trial[6 * i + 5]));

		if (i == 0)
		{
			positions.push_back(vec3(trial[0], trial[1], trial[2]));
			sizes.push_back(vec3(trial[3], trial[4], trial[5]));
		}
		else if (i == 1)
		{
			positions.push_back(vec3(positions[0].x + sizes[0].x, positions[0].y, positions[0].z));
			sizes.push_back(vec3(trial[6], trial[7], trial[8]));
		}
		else if (i == 2)
		{
			positions.push_back(vec3(positions[1].x + sizes[1].x - trial[9], positions[1].y - trial[10], positions[1].z));
			sizes.push_back(vec3(trial[9], trial[10], trial[11]));
		}

		boxEnds.push_back(positions[i]);
		boxEnds.push_back(positions[i] + sizes[i]);
	}

	p.new_proc(positions, sizes);

	//p.create_default_distance_field();
	//const vector<double>& df = p.df();
	//
	//for (int i = 0; i < df.size(); ++i)
	//{
	//	error_df += (m_df[i] - df[i]) * (m_df[i] - df[i]);
	//}

	//qDebug() << error;


	//Need to rethink if this can hurt the optimization. 
	for (int i = 0; i < m_dimension; i++)
	{
		if (trial[i] < 0)
		{//		return 10000;
			error_negative += 10000 * abs(trial[i]);
			//error_negative = 100000;
		}
	}


	for (int i = 0; i < m_samples.size(); i++)
	{
		error_project += p.project_point_to_model(m_samples[i]);
	}

    double errorBoundingBox = 0.0;//errorBB(p, boxEnds);

	double errorDisEnt = errorDisentangle(p);

	vector<vector<int>> v;
	vector<int> s, t;
	s.push_back(0);
	s.push_back(2);
	s.push_back(1);
	s.push_back(1);
	t.push_back(1);
	t.push_back(2);
	t.push_back(2);
	t.push_back(4);
	v.push_back(s);
	v.push_back(t);

	double error_Grammar = 0.0;// errorGrammar(v, boxEnds);

	double error = (params::inst()->weightProjection * error_project +
		params::inst()->weightNegative * error_negative +
		params::inst()->weightBoundingBox * errorBoundingBox
		+ params::inst()->weightDistanceField * error_df
		+ params::inst()->weightGrammar * error_Grammar
		+ +params::inst()->weightDisentangle * errorDisEnt);

	params::inst()->error = error;

	return error;
}

double Optimizer::energy_function_grammar(double trial[], bool &foundSolution)
{
	double error_project = 0.0, error_negative = 0.0, error_df = 0.0;
	Proc_Model p;

	m_grammar->softCleanUp();
	m_grammar->adjustParameters(trial, m_grammar->m_nrIndependentVariables);
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

	//Need to rethink if this can hurt the optimization. 
	for (int i = 0; i < sizes.size(); i++)
	{
		if (sizes[i].x < 0)
		{
			error_negative += 10000 * abs(sizes[i].x);
		}
		if (sizes[i].y < 0)
		{
			error_negative += 10000 * abs(sizes[i].y);
		}
		if (sizes[i].z < 0)
		{
			error_negative += 10000 * abs(sizes[i].z);
		}
	}

    //if(error_negative > 0)
    //    qDebug() << error_negative;

	for (int i = 0; i < m_samples.size(); i++)
	{
		error_project += p.project_point_to_model(m_samples[i]);
	}

    vector<vec3> bbox = p.bounding_box(m_samples);
	double errorBoundingBox = errorBB(p, bbox, boxEnds);

    double errorSkin = errorSkinning(p);

	//double errorDisEnt = errorDisentangle(p);

    double error = (params::inst()->weightProjection * error_project * error_project +
        params::inst()->weightNegative * error_negative * error_negative +
        params::inst()->weightBoundingBox * errorBoundingBox * errorBoundingBox
        + params::inst()->weightSkin * errorSkin * errorSkin);

	params::inst()->error = error;

	return 1000*error/m_samples.size();
}

double Optimizer::errorBB(Proc_Model &p, vector<vec3> &bbox, vector<vec3> &boxEnds)
{
	double error_bbox_vols = 0.0;

	vector<vec3> bbox_trial = p.bounding_box(boxEnds);
	vec3 outerbox = bbox[1], innerbox = bbox[0], outerbox_trial = bbox_trial[1], innerbox_trial = bbox_trial[0];
	double volume_bbox = (outerbox.x - innerbox.x)*(outerbox.y - innerbox.y)*(outerbox.z - innerbox.z);
	double volume_bbox_trial = (outerbox_trial.x - innerbox_trial.x)*(outerbox_trial.y - innerbox_trial.y)*(outerbox_trial.z - innerbox_trial.z);
	error_bbox_vols = abs(volume_bbox - volume_bbox_trial);
	
    return error_bbox_vols;
}

double Optimizer::errorGrammar(vector<vector<int>> &v, vector<vec3> &boxEnds)
{
	double error = 0.0;
	for (int i = 0; i < v.size(); i++)
	{
		vector<int> order = v[i];
		vec3 point1 = conventionPoint(boxEnds[order[0] * 2], boxEnds[order[0] * 2 + 1], order[1]);
		vec3 point2 = conventionPoint(boxEnds[order[2] * 2], boxEnds[order[2] * 2 + 1], order[3]);

		error += dot(point1 - point2, point1 - point2);
	}

	return error;
}

double Optimizer::errorDisentangle(Proc_Model &p)
{
	double disentanglement = 0.0;
	vector<vec3> proc_pos = p.position(), proc_size = p.size();
	for (int i = 0; i < proc_pos.size(); i++)
	{
		vec3 pos_i = proc_pos[i], size_i = proc_size[i];
		for (int j = i + 1; j < proc_pos.size(); j++)
		{
			vec3 pos_j = proc_pos[j], size_j = proc_size[j];
			vec3 min_box = vec3(max(pos_i.x, pos_j.x), max(pos_i.y, pos_j.y), max(pos_i.z, pos_j.z));
			vec3 max_box = vec3(min(pos_i.x + size_i.x, pos_j.x + size_j.x), min(pos_i.y + size_i.y, pos_j.y + size_j.y), min(pos_i.z + size_i.z, pos_j.z + size_j.z));
			disentanglement += max(0.0f, max_box.x - min_box.x)*max(0.0f, max_box.y - min_box.y)*max(0.0f, max_box.z - min_box.z);
		}
	}

	return disentanglement;
}

double Optimizer::errorSkinning(Proc_Model &p)
{
    double skinerror = 0.0;

    vector<vec3> proc_pos = p.position(), proc_size = p.size();
    for (int i = 0; i < proc_pos.size(); i++)
    {
        vec3 size_i = proc_size[i];
        skinerror += abs(size_i.x*size_i.y*size_i.z);//dot(size_i, size_i);
    }

    return 100.0*skinerror;
}

vec3 Optimizer::conventionPoint(vec3 &start, vec3 &end, int num)
{
	vec3 v;
	switch (num){
	case 1: v = start;                            break;
	case 2: v = vec3(end.x, start.y, start.z);	  break;
	case 3: v = vec3(start.x, end.y, start.z);	  break;
	case 4: v = vec3(end.x, end.y, start.z);		  break;
	case 5: v = vec3(start.x, start.y, end.z);	  break;
	case 6: v = vec3(end.x, start.y, end.z);	  break;
	case 7: v = vec3(start.x, end.y, end.z);	  break;
	case 8: v = end;							  break;
	}

	return v;
}

void Optimizer::setRange(const vector<double> &mi, const vector<double> &ma)
{
    m_minRange = mi;
    m_maxRange = ma;
}

void Optimizer::setDistanceField(const vector<double> &df)
{
    m_df = df;
}

void Optimizer::setSampleVals(const vector<vec3> &samples)
{
	m_samples = samples;
}

void Optimizer::setGrammar(Grammar *gram)
{
	m_grammar = gram;
}
