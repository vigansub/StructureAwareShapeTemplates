#include "SimulatedAnnealing.h"
#include "Proc_Model.h"
#include "Node.h"
#include "Scene.h"


std::default_random_engine gen(time(0));
std::normal_distribution<double> s(0.0, 3.0);

Sim_Anneal::Sim_Anneal(int dimension, Scene* scene)
{
    srand(time(0));
    m_dimension = dimension;
    for (int i = 0; i < m_dimension; i++)
        m_trial.push_back(rand(0.0, 1.0));

    m_old = m_trial;

    m_scene = scene;
}

Sim_Anneal::~Sim_Anneal()
{
}

vector<double> Sim_Anneal::optimize()
{

	vector<double> old_sol, sol;
	old_sol = m_trial;
	sol = old_sol;
	int j = 1;
	double old_cost = energy_function(); //Cost computation
	double T = 1.0f;
	double alpha = 0.999f, T_min = 0.0001f;
    srand(time(0));

	while (T > T_min)
	{
		int i = 0;
		while (i < 100)
		{
            generateTrial();
			vector<double> new_sol = m_trial;
			double new_cost = energy_function();
			double ap = acceptance_probability(old_cost, new_cost, T);
            double k = rand(0.0, 1.0);
            if ((ap > log(k)))
            {
         //       std::cout << "Old Cost: " << old_cost << "New Cost: " << new_cost << std::endl;
                m_old = new_sol;
                old_cost = new_cost;
            }
			i++;
		}
		T *= alpha;
        std::cout << endl;
        std::cout << T << ":" << old_cost << ":";
        for (int i = 0; i < m_old.size(); i++)
            std::cout << m_old[i] << " ";

        if (m_scene)
            m_scene->updateGrammarGeometry_SimAnneal();
	}

	return m_old;
}

double Sim_Anneal::acceptance_probability(double old_cost, double new_cost, double T)
{
    double a = (old_cost-new_cost)/(1.0*T);
    return a;
}

double Sim_Anneal::energy_function()
{
    double error_project = 0.0, error_negative = 0.0, error_df = 0.0;
    Proc_Model p;

    m_grammar->softCleanUp();
    m_grammar->adjustParameters(m_old.data(), m_grammar->m_nrIndependentVariables);
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

    for (int i = 0; i < m_samples.size(); i++)
    {
        error_project += p.project_point_to_model(m_samples[i]);
    }

    double errorBoundingBox = errorBB(p, boxEnds);

    double errorSkin = errorSkinning(p);

    double errorDisEnt = errorDisentangle(p);

    double error = (params::inst()->weightProjection * error_project +
        params::inst()->weightNegative * error_negative +
        params::inst()->weightBoundingBox * errorBoundingBox
        + params::inst()->weightDistanceField * error_df
        + params::inst()->weightDisentangle * errorDisEnt
        + params::inst()->weightSkin * errorSkin);

    params::inst()->error = error;

    //qDebug() << error_project/error << errorBoundingBox/error << error_negative/error << errorSkin/error;

    return error / m_samples.size();

}

double Sim_Anneal::errorBB(Proc_Model &p, vector<vec3> &boxEnds)
{
    double error_bbox_vols = 0.0;

    vector<vec3> bbox = p.bounding_box(m_samples), bbox_trial = p.bounding_box(boxEnds);
    vec3 outerbox = bbox[1], innerbox = bbox[0], outerbox_trial = bbox_trial[1], innerbox_trial = bbox_trial[0];
    double volume_bbox = (outerbox.x - innerbox.x)*(outerbox.y - innerbox.y)*(outerbox.z - innerbox.z);
    double volume_bbox_trial = (outerbox_trial.x - innerbox_trial.x)*(outerbox_trial.y - innerbox_trial.y)*(outerbox_trial.z - innerbox_trial.z);
    error_bbox_vols = abs(volume_bbox - volume_bbox_trial);

    return error_bbox_vols;
}

double Sim_Anneal::errorDisentangle(Proc_Model &p)
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

double Sim_Anneal::errorSkinning(Proc_Model &p)
{
    double skinerror = 0.0;

    vector<vec3> proc_pos = p.position(), proc_size = p.size();
    for (int i = 0; i < proc_pos.size(); i++)
    {
        vec3 size_i = proc_size[i];
        skinerror += dot(size_i, size_i);
        //skinerror += abs(size_i.x*size_i.y*size_i.z);
    }

    return skinerror;
}

void Sim_Anneal::setSampleVals(const vector<vec3> &samples)
{
    m_samples = samples;
}

void Sim_Anneal::generateTrial()
{
    srand(time(0));
    for (int i = 0; i < m_dimension; i++)
    {   // m_trial[i] = rand(0.0, 3.0);
        m_trial[i] = m_old[i] + s(gen);
    }
}

void Sim_Anneal::setGrammar(Grammar *grammar)
{
    m_grammar = grammar; 
}

/*
vector<double> Sim_Anneal::sim_ann_normal_dist(vector<vector<double>> samples, gaussParams g)
{
	vector<double> old_sol, sol;
	Eigen::MatrixXd variance = g.variance();

	Eigen::LLT<Eigen::MatrixXd> lltOfA(variance); // compute the Cholesky decomposition of A
	Eigen::MatrixXd L = lltOfA.matrixL();

	old_sol = samples[0];
	sol = old_sol;
	int j = 1;
	double old_cost;
	int num_tests = 10;
	old_cost = 0.0f;
	for (int i = 0; i < num_tests; i++)
	{
		vector<double> old_sol_sample;
		old_sol_sample = samples_gen(old_sol, L);
		old_cost += g.prob_val(old_sol_sample); //Cost computation
	}
	old_cost /= num_tests;

	double T = 1.0f;
	double alpha = 0.99f, T_min = 0.1f;

	while (T > T_min)
	{
		int i = 0;
		while (i < 100)
		{
			vector<double> new_sol = samples[j];
			j++;
			double new_cost = 0.0f;
			for (int i = 0; i < num_tests; i++)
			{
				vector<double> new_sol_sample;
				new_sol_sample = samples_gen(new_sol, L);
				new_cost += g.prob_val(new_sol_sample); //Cost computation
			}
			new_cost /= num_tests;

			double ap = acceptance_probability(old_cost, new_cost, T);
			if (ap > rand(0, 1))
			{
				sol = new_sol;
				old_cost = new_cost;
			}
			i++;
		}
		T *= alpha;
	}
	return sol;
}

vector<double> Sim_Anneal::sim_ann_normal_dist_var(vector<vector<double>> samples, gaussParams g)
{
	vector<double> old_sol, sol;
	Eigen::MatrixXd variance = g.variance();

	Eigen::LLT<Eigen::MatrixXd> lltOfA(variance); // compute the Cholesky decomposition of A
	Eigen::MatrixXd L = lltOfA.matrixL();

	old_sol = samples[0];
	sol = old_sol;
	int j = 1;
	double old_cost;
	int num_tests = 100;
	old_cost = 0.0f;
	for (int i = 0; i < num_tests; i++)
	{
		vector<double> old_sol_sample;
		old_sol_sample = samples_gen_var(old_sol);
		old_cost += g.prob_val(old_sol_sample); //Cost computation
	}
	old_cost /= num_tests;

	double T = 1.0f;
	double alpha = 0.9f, T_min = 0.0001f;

	while ((T > T_min) && (j<samples.size()))
	{
		int i = 0;
		while (i < 100)
		{
			vector<double> new_sol = samples[j];
			j++;
			double new_cost = 0.0f;
			for (int i = 0; i < num_tests; i++)
			{
				vector<double> new_sol_sample;
				new_sol_sample = samples_gen_var(new_sol);
				new_cost += g.prob_val(new_sol_sample); //Cost computation
			}
			new_cost /= num_tests;

			double ap = acceptance_probability(old_cost, new_cost, T);
			if (ap > rand(0, 1))
			{
				sol = new_sol;
				old_cost = new_cost;
			}
			i++;

			if (j > samples.size())
				return sol;
		}
		T *= alpha;
	}
	return sol;
}

vector<double> Sim_Anneal::sim_ann_dist_field(MH_MC m, const float* data)
{
	vector<double> old_sol, sol;
	mat4 origToVoxel = m.origToVoxel();
	vector<vector<double>> samples = m.samples();

	old_sol = samples[0];
	sol = old_sol;
	int j = 1;
	
	vec4 samples_start = vec4(sol[0], sol[2], sol[4], 1);
	//samples_start = vec4(0.0f, 0.0f, 0.0f, 1);
	vec4 startShift = origToVoxel * samples_start;
	vec4 samples_end = vec4(sol[0] + sol[1], sol[2] + sol[3], sol[4] + sol[5], 1);
	//samples_end = vec4(0.0f, 0.0f, 0.0f, 1);
	vec4 endShift = origToVoxel * samples_end;
	double old_cost = m.Df_Cube(startShift.x / startShift.w, endShift.x / endShift.w, startShift.y / startShift.w, endShift.y / endShift.w, startShift.z / startShift.w, endShift.z / endShift.w, data);

	double T = 1.0f;
	double alpha = 0.9f, T_min = 0.01f;

	while (T > T_min)
	{
		int i = 0;
		while (i < 100)
		{
			vector<double> new_sol = samples[j];
			j++;

			vec4 samples_start = vec4(new_sol[0], new_sol[2], new_sol[4], 1);
			//samples_start = vec4(0.0f, 0.0f, 0.0f, 1);
			vec4 startShift = origToVoxel * samples_start;
			vec4 samples_end = vec4(new_sol[0] + new_sol[1], new_sol[2] + new_sol[3], new_sol[4] + new_sol[5], 1);
			//samples_end = vec4(0.0f, 0.0f, 0.0f, 1);
			vec4 endShift = origToVoxel * samples_end;
			double new_cost = m.Df_Cube(startShift.x / startShift.w, endShift.x / endShift.w, startShift.y / startShift.w, endShift.y / endShift.w, startShift.z / startShift.w, endShift.z / endShift.w, data);

			double ap = acceptance_probability(old_cost, new_cost, T);
			if (ap > rand(0, 1))
			{
				sol = new_sol;
				old_cost = new_cost;
			}
			i++;
		}
		T *= alpha;
	}

	return sol;
}

vector<double> Sim_Anneal::sim_ann_dist_field_self(MH_MC m, Proc_Model proc)
{
	vector<double> old_sol, sol, best_sol;
	mat4 origToVoxel = m.origToVoxel();
	vector<vector<double>> samples = m.samples();
	vector<double> data = proc.df();
	int num_boxes = proc.position().size();
	double best_cost;

	for (int iter = 0; iter < 10; iter++)
	{
		old_sol = samples[iter*10000 + 0];
		sol = old_sol;
		int j = 1;

		// vec4 samples_start = vec4(sol[0], sol[2], sol[4], 1);
		// vec4 samples_end = vec4(sol[0] + sol[1], sol[2] + sol[3], sol[4] + sol[5], 1);

		vector<vec3> sample_instance, sample_size;
		for (int box_count = 0; box_count < num_boxes; box_count++)
		{
			vec3 sample_curr_instance;
			vec3 sample_curr_size;
			int ival = box_count * 6;
			sample_curr_instance = vec3(sol[ival], sol[ival + 2], sol[ival + 4]);
			sample_curr_size = vec3(sol[ival + 1], sol[ival + 3], sol[ival + 5]);

			sample_instance.push_back(sample_curr_instance);
			sample_size.push_back(sample_curr_size);
		}

		//double old_cost = proc.evaluate_df_val(sample_instance, sample_size);
		double old_cost = proc.compare_df(sample_instance, sample_size);
		

		double T = 1.0f;
		double alpha = 0.9f, T_min = 0.01f;

		while (T > T_min)
		{
			int i = 0;
			while (i < 100)
			{
				vector<double> new_sol = samples[iter * 10000 + j];
				j++;

				// vec4 samples_start = vec4(new_sol[0], new_sol[2], new_sol[4], 1);
				// vec4 samples_end = vec4(new_sol[0] + new_sol[1], new_sol[2] + new_sol[3], new_sol[4] + new_sol[5], 1);

				vector<vec3> sample_instance, sample_size;
				for (int box_count = 0; box_count < num_boxes; box_count++)
				{
					vec3 sample_curr_instance;
					vec3 sample_curr_size;
					int ival = box_count * 6;
					sample_curr_instance = vec3(new_sol[ival], new_sol[ival + 2], new_sol[ival + 4]);
					sample_curr_size = vec3(new_sol[ival + 1], new_sol[ival + 3], new_sol[ival + 5]);
					sample_instance.push_back(sample_curr_instance);
					sample_size.push_back(sample_curr_size);
				}

				//double new_cost = proc.evaluate_df_val(sample_instance, sample_size);
				double new_cost = proc.compare_df(sample_instance, sample_size);

				double ap = acceptance_probability(old_cost, new_cost, T);
				if (ap > rand(0, 1))
				{
					sol = new_sol;
					old_cost = new_cost;
				}
				i++;
			}
			T *= alpha;
		}

		if (iter == 0)
		{
			best_sol = sol;
			best_cost = old_cost;
		}
		else if (old_cost < best_cost)
		{
			best_sol = sol;
			best_cost = old_cost;
		}
	}

	return best_sol;
}

vector<double> Sim_Anneal::optimization(MH_MC m, Proc_Model proc)
{
	vector<vector<double>> samples = m.samples();
	
	vector<double> old_sol = samples[0];
	vector<double> sol = old_sol;
	int num_boxes = proc.position().size();

	vector<vec3> sample_instance, sample_size;
	for (int box_count = 0; box_count < num_boxes; box_count++)
	{
		vec3 sample_curr_instance;
		vec3 sample_curr_size;
		int ival = box_count * 6;
		sample_curr_instance = vec3(sol[ival], sol[ival + 2], sol[ival + 4]);
		sample_curr_size = vec3(sol[ival + 1], sol[ival + 3], sol[ival + 5]);

		sample_instance.push_back(sample_curr_instance);
		sample_size.push_back(sample_curr_size);
	}
	double old_cost = proc.compare_df(sample_instance, sample_size);

	for (int i = 1; i < samples.size(); i++)
	{
		
		vector<double> new_sol = samples[i];
		vector<vec3> sample_instance, sample_size;
		for (int box_count = 0; box_count < num_boxes; box_count++)
		{
			vec3 sample_curr_instance;
			vec3 sample_curr_size;
			int ival = box_count * 6;
			sample_curr_instance = vec3(new_sol[ival], new_sol[ival + 2], new_sol[ival + 4]);
			sample_curr_size = vec3(new_sol[ival + 1], new_sol[ival + 3], new_sol[ival + 5]);
			sample_instance.push_back(sample_curr_instance);
			sample_size.push_back(sample_curr_size);
		}

		double new_cost = proc.compare_df(sample_instance, sample_size);

		if (new_cost < old_cost)
		{
			old_cost = new_cost;
			sol = new_sol;
		}
	}

	return sol;
}


double Sim_Anneal::acceptance_probability(double old_cost, double new_cost, double T)
{
	double a = exp(-(new_cost - old_cost) / T);
	return a;
}

vector<double> Sim_Anneal::samples_gen(vector<double> mean, Eigen::MatrixXd L)
{
	Eigen::VectorXd mean1;
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0, 1.0);

	mean1.resize(mean.size(), 1);
	for (int i = 0; i < mean.size(); i++)
	{
		mean1[i] = mean[i];
	}
	Eigen::VectorXd k;
	k.resize(mean.size(), 1);
	for (int i = 0; i < mean.size(); i++)
	{
		k[i] = distribution(generator);
	}
	k = mean1 + L*k;
	vector<double> k1;
	for (int i = 0; i < mean.size(); i++)
		k1.push_back(k[i]);

	return k1;
}

vector<double> Sim_Anneal::samples_gen_var(vector<double> inp)
{
	vector<double> k;
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0, 1.0);
	for (int i = 0; i < inp.size() / 2; i++)
	{
		double mean = inp[2 * i];
		double var = inp[2 * i + 1];
		double std = sqrt(var);
		k.push_back(mean + std*distribution(generator));
	}
	return k;
}*/