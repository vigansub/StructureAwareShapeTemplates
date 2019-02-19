#include "Swarm.h"
#include "Proc_Model.h"
#include "Node.h"
#include "Scene.h"
//#include <ppl.h>

Swarm::Swarm(Scene *scene, Grammar* grammar, int particles, int iterations, vector<double> min, vector<double> max, double w, double phi_p, double phi_g)
{
    srand(time(0));
    m_scene = scene;
    m_grammar = grammar;
    m_particles = particles;
    m_iter = iterations;
    m_mi = min;
    m_ma = max;
    m_dimensions = m_mi.size();
    m_parameters.push_back(w);
    m_parameters.push_back(phi_p);
    m_parameters.push_back(phi_g);

    double best = 0.0;

    for (int i = 0; i < m_particles; i++)
    {
        vector<double> particle_i, velocity_i;
        particle_i.resize(m_dimensions);
        velocity_i.resize(m_dimensions);
        for (int j = 0; j < m_dimensions; j++)
        {
            particle_i[j] = rand(m_mi[j], m_ma[j]);
            double range = abs(m_ma[j] - m_mi[j]);
            velocity_i[j] = 0.01*particle_i[j];// rand(-range, range);
        }
        m_positions.push_back(particle_i);
        m_velocities.push_back(velocity_i);

        if (i == 0)
        {
            m_global = particle_i;
            best = energy_function(m_global);
        }
        else if (energy_function(particle_i) < best)
        {
            m_global = particle_i;
            best = energy_function(m_global);
        }
    }
    m_locals = m_positions;

}

Swarm::~Swarm()
{

}

double Swarm::energy_function(vector<double> position)
{
    double error_project = 0.0, error_negative = 0.0, error_df = 0.0;
    Proc_Model p;

    //Grammar *g = new Grammar;
    //
    //memcpy(g, m_grammar, sizeof(m_grammar));
    //
    //g->softCleanUp();
    //g->adjustParameters(position.data(), g->m_nrIndependentVariables);
    //g->derive(); 
    //
    //QMultiMap<QString, Node*> nodes = g->nodes();

    m_grammar->softCleanUp();
    m_grammar->adjustParameters(position.data(), m_grammar->m_nrIndependentVariables);
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
        error_project += p.project_point_to_model(m_samples[i]);

    vector<vec3> bbox = p.bounding_box(m_samples);
    double errorBoundingBox = errorBB(bbox, p, boxEnds);

    double errorSkin = errorSkinning(p);

    double errorDisEnt = errorDisentangle(p);

    double error = (params::inst()->weightProjection * error_project +
        params::inst()->weightNegative * error_negative +
        params::inst()->weightBoundingBox * errorBoundingBox
        + params::inst()->weightDisentangle * errorDisEnt
        + params::inst()->weightSkin * errorSkin);

    error /= m_samples.size();

    params::inst()->error = error;

    //delete g;

    return 1000.0*error;
}

void Swarm::optimize()
{
/*
    srand(time(0));
    double w = m_parameters[0], phi_p = m_parameters[1], phi_g = m_parameters[2];
    int stall = 0;

    //Concurrency::parallel_for(0, n,
    //    [](int i)
    //{
    //    a[i] += a[i - 1]; // incorrect! 
    //});

    for (int iter = 0; iter < m_iter; iter++)
    {
        if (iter<2500)
           w = 0.4 + 0.5*iter/2500;
        
        //phi_p = -2.5*iter / m_iter + 3.0;
        //phi_g = 2.5*iter / m_iter + 0.5;

        Concurrency::parallel_for(
            0, m_particles-1, [&](int i)
     //       for (int i = 0; i < m_particles; i++)
            {
                vector<double> velocity_i = m_velocities[i], position_i = m_positions[i], local_i = m_locals[i];
                for (int d = 0; d < m_dimensions; d++)
                {
                    double rp = rand(0.0, 1.0), rg = rand(0.0, 1.0);
                    velocity_i[d] = velocity_i[d] * w + (local_i[d] - position_i[d]) * rp * phi_p + (m_global[d] - position_i[d]) * rg * phi_g;
                    position_i[d] += velocity_i[d];
                }

                m_velocities[i] = velocity_i;
                m_positions[i] = position_i;
                
                //if (energy_function(m_positions[i]) < energy_function(m_locals[i]))
                //{
                //    m_locals[i] = m_positions[i];
                //}
            }
        );

        for (int i = 0; i < m_particles; i++)
        {
            if (energy_function(m_positions[i]) < energy_function(m_locals[i]))
            {
                m_locals[i] = m_positions[i];
            }

            if (energy_function(m_locals[i]) < energy_function(m_global))
                m_global = m_locals[i];
        }

        double best_val = energy_function(m_global);

        std::cout << iter << ": Error = " << best_val << ": ";
        for (int i = 0; i < m_dimensions; i++)
            std::cout << m_global[i] << " ";

        std::cout << endl;

        if (m_scene)
            m_scene->updateGrammarGeometry();
    }
*/
}

vector<double> Swarm::global_opt()
{
    return m_global;
}

void Swarm::setSamples(vector<vec3> samplevals)
{
    m_samples = samplevals;
}

void Swarm::setGrammar(Grammar *grammar)
{
    m_grammar = grammar;
}

double Swarm::errorBB(vector<vec3> &bbox, Proc_Model &p, vector<vec3> &boxEnds)
{
    double error_bbox_vols = 0.0;

    vector<vec3> bbox_trial = p.bounding_box(boxEnds);
    vec3 outerbox = bbox[1], innerbox = bbox[0], outerbox_trial = bbox_trial[1], innerbox_trial = bbox_trial[0];
    double volume_bbox = (outerbox.x - innerbox.x)*(outerbox.y - innerbox.y)*(outerbox.z - innerbox.z);
    double volume_bbox_trial = (outerbox_trial.x - innerbox_trial.x)*(outerbox_trial.y - innerbox_trial.y)*(outerbox_trial.z - innerbox_trial.z);
    error_bbox_vols = abs(volume_bbox - volume_bbox_trial);

    return error_bbox_vols;
}

double Swarm::errorDisentangle(Proc_Model &p)
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

double Swarm::errorSkinning(Proc_Model &p)
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
