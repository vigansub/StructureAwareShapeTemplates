#include "Genetic_Algorithm.h"
#include "Node.h"
#include "Scene.h"
//#include <ppl.h>

Genetic::Genetic(Scene *scene, int iterations, int bits, int population, vector<double> min, vector<double> max, double crossover, double mutation)
{
    m_scene = scene;
    m_iter = iterations;
    m_numbits = bits;
    m_population = population;
    m_mi = min, m_ma = max;
    m_dimensions = m_mi.size();
    m_crossover = crossover; 
    m_mutation = mutation;

    m_chromosomes.resize(m_population);
    m_fitness.resize(m_population);
}

Genetic::~Genetic()
{

}

double Genetic::energy_function(vector<double> position)
{
    double error_project = 0.0, error_negative = 0.0, error_df = 0.0;
    Proc_Model p;

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
    
    return 1.0/error;
}

void Genetic::optimize()
{
    srand(time(0));

    vector<double> k;
    k.resize(m_dimensions);
    for (int i = 0; i < m_population; ++i)
    {
        for (int j = 0; j < m_dimensions; ++j)
        {
            k[j] = rand(m_mi[j],m_ma[j]);
        }
        m_chromosomes[i] = k;
        m_fitness[i] = energy_function(k);
    }

    m_fitness = norm_fitness(m_fitness);
    m_best = m_chromosomes[find_best_fitness()];

    for (int i = 0; i < m_iter; i++)
    {
        int new_tot = 0;
        vector<vector<double>> p;
        p.resize(m_population);
        vector<double> p1, p2, new_fitness;
        while (new_tot < m_population)
        {
            int k1 = selection(), k2 = selection();
            p1 = m_chromosomes[k1], p2 = m_chromosomes[k2];
            if (rand(0.0,1.0) < m_crossover)
            {
                vector<int> k;
                k.resize(m_dimensions);
                for (int j = 0; j < m_dimensions; ++j)
                {
                    k[j] = rand() % m_numbits;
                }
                swap(p1, p2, k);
            }
            p[new_tot] = p1;
            new_fitness.push_back(energy_function(p1));
            new_tot++;
            p[new_tot] = p2;
            new_fitness.push_back(energy_function(p2));
            new_tot++;
        }
        m_fitness = norm_fitness(new_fitness);
        m_chromosomes = p;

        if (energy_function(m_chromosomes[find_best_fitness()]) > energy_function(m_best))
            m_best = m_chromosomes[find_best_fitness()];

        double best_val = 1000.0/(energy_function(m_best));

        std::cout << i << ": Error = " << best_val << ": ";
        for (int i = 0; i < m_dimensions; i++)
            std::cout << m_best[i] << " ";

        std::cout << endl;

        if (m_scene)
            m_scene->updateGrammarGeometry();

    }
}

void Genetic::setSamples(vector<vec3> samplevals)
{
    m_samples = samplevals;
}

void Genetic::setGrammar(Grammar *grammar)
{
    m_grammar = grammar;
}

double Genetic::errorBB(vector<vec3> &bbox, Proc_Model &p, vector<vec3> &boxEnds)
{
    double error_bbox_vols = 0.0;

    vector<vec3> bbox_trial = p.bounding_box(boxEnds);
    vec3 outerbox = bbox[1], innerbox = bbox[0], outerbox_trial = bbox_trial[1], innerbox_trial = bbox_trial[0];
    double volume_bbox = (outerbox.x - innerbox.x)*(outerbox.y - innerbox.y)*(outerbox.z - innerbox.z);
    double volume_bbox_trial = (outerbox_trial.x - innerbox_trial.x)*(outerbox_trial.y - innerbox_trial.y)*(outerbox_trial.z - innerbox_trial.z);
    error_bbox_vols = abs(volume_bbox - volume_bbox_trial);

    return error_bbox_vols;
}

double Genetic::errorDisentangle(Proc_Model &p)
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

double Genetic::errorSkinning(Proc_Model &p)
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

int Genetic::selection()
{
    vector<double> fitness_cumsum;
    fitness_cumsum.resize(m_fitness.size());
    fitness_cumsum[0] = m_fitness[0];
    for (int i = 1; i < m_fitness.size(); ++i)
        fitness_cumsum[i] = fitness_cumsum[i - 1] + m_fitness[i];

    double k = rand(0.0, 1.0);
    for (int i = 0; i < m_fitness.size() - 1; ++i)
    {
        if (k <= fitness_cumsum[0])
            return 0;
        else if ((k > fitness_cumsum[i]) && (k <= fitness_cumsum[i + 1]))
            return i + 1;
    }

    return m_fitness.size() - 1;
}

double Genetic::unif_rand()
{
    int s = rand();
    double p = (1.0 * s / RAND_MAX * 1.0);
    return p;
}

vector<double> Genetic::norm_fitness(vector<double> &fitness)
{
    double tot_fitness = 0.0;

    for (int i = 0; i < fitness.size(); ++i)
        tot_fitness += fitness[i];

    for (int i = 0; i < fitness.size(); ++i)
        fitness[i] /= tot_fitness;

    double k = fitness[0];
    for (int i = 0; i < fitness.size(); ++i)
        k = max(k, fitness[i]);

    return fitness;
}

void Genetic::swap(vector<double> &a, vector<double> &b, vector<int> p)
{
    vector<vector<double>> a_bin = dec2bin(a), b_bin = dec2bin(b);
    for (int i = 0; i < a.size(); ++i)
    {
        vector<double> a_i = a_bin[i], b_i = b_bin[i];
        int p_i = p[i];
        for (int j = p_i; j < a_i.size(); ++j)
        {
            double c = a_i[j];
            a_i[j] = b_i[j];
            b_i[j] = c;
        }
        mutate(a_i, m_mutation);
        mutate(b_i, m_mutation);
        a_bin[i] = a_i, b_bin[i] = b_i;
    }
    a = bin2dec(a_bin), b = bin2dec(b_bin);
}

void Genetic::mutate(vector<double> &a, double mutation)
{
    for (int i = 0; i < a.size(); ++i)
    {
        if (rand(0.0, 1.0) < mutation)
        {
            if (a[i] == 0)
                a[i] = 1;
            else a[i] = 0;
        }
    }
}

int Genetic::find_best_fitness()
{
    int best_fitness = 0;
    for (int i = 1; i < m_fitness.size(); i++)
        if (m_fitness[i] > m_fitness[best_fitness])
            best_fitness = i;

    return best_fitness;
}

vector<double> Genetic::bin2dec(vector<vector<double>> &a)
{
    vector<double> op;
    op.resize(m_dimensions);
    for (int i = 0; i < m_dimensions; ++i)
    {
        double sum = 0.0;
        vector<double> k = a[i];
        for (int j = 0; j < m_numbits; ++j)
        {
            sum *= 2.0;
            sum += k[j];
        }
        op[i] = m_mi[i] + (m_ma[i] - m_mi[i])*sum / (pow(2, m_numbits));
    }

    return op;
}

vector<vector<double>> Genetic::dec2bin(vector<double> &a)
{
    vector<vector<double>> op;
    op.resize(m_dimensions);
    vector<double> k;
    k.resize(m_numbits);

    for (int i = 0; i < m_dimensions; ++i)
    {
        int in = (a[i]-m_mi[i])*pow(2, m_numbits)/(m_ma[i]-m_mi[i]);
        for (int j = 0; j < m_numbits; ++j)
        {
            k[m_numbits - j - 1] = (in % 2)*1.0;
            in /= 2;
        }
        op[i] = k;
    }
    return op;
}

vector<double> Genetic::global_opt()
{
    return m_best;
}
