#include "CMAESWrapper.h"
#include <ctime>
#include "Proc_Model.h"
#include "Grammar.h"
#include "Node.h"
#include "Scene.h"

CMAESWrapper::CMAESWrapper(int dimension, Scene* scene)
{
    m_dimension = dimension;
    m_scene = scene;
    m_besterror = 10000.0;
    //m_grammar = NULL;
}

CMAESWrapper::~CMAESWrapper()
{
}

double CMAESWrapper::fitfun(double const *x)
//double CMAESWrapper::fitfun(vector<double> x, int N, vector<double> dim)
{

    double error_project = 0.0, error_negative = 0.0, error_df = 0.0;
    Proc_Model p;

    m_grammar->softCleanUp();
    m_grammar->adjustParameters(x, m_grammar->m_nrIndependentVariables);
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
    int num_boxes = positions.size();

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
    error_project *= (3000.0 / m_samples.size());

    double error_reverse = 0.0;// p.project_model_to_pc(m_samples);
    //double error_reverse = p.project_model_to_pc(m_samples);
    vector<vec3> bbox = p.bounding_box(m_samples);
    double errorBoundingBox = errorBB(bbox, p, boxEnds);
    //double errorBoundingBox = 0.0;// errorBB(bbox, p, boxEnds);

    double errorSkin = errorSkinning(p);
    //double errorSkin = 0.0;// errorSkinning(p);

    //errorBoundingBox = 100.0*abs(p.volume() - m_volume);

    double errorDisEnt = errorDisentangle(p);

    double penalty_error = errorPenalty(p);

    double error = (params::inst()->weightProjection * error_project*error_project +
        params::inst()->weightProjection * error_reverse*error_reverse +
        params::inst()->weightNegative * error_negative*error_negative +
        params::inst()->weightBoundingBox * errorBoundingBox*errorBoundingBox
        + params::inst()->weightSkin * errorSkin*errorSkin
        + params::inst()->weightDisentangle*errorDisEnt*errorDisEnt);

    //if (params::inst()->isPartial)
    //    error += 10000 *(abs(errorBoundingBox));

    //qDebug() << errorBoundingBox;

    vector<double> iter_error_vals;
    iter_error_vals.resize(5);
    iter_error_vals[0] = error;
    iter_error_vals[1] = error_project;
    iter_error_vals[2] = error_negative;
    iter_error_vals[3] = errorBoundingBox;
    iter_error_vals[4] = errorSkin;
    local_error.push_back(iter_error_vals);
    //error += 0.25 * p.position().size();
    params::inst()->error = error;

    //qDebug() << p.volume() << " " << (bbox[1].x - bbox[0].x)*(bbox[1].y - bbox[0].y)*(bbox[1].z - bbox[0].z);

    //qDebug() << errorBoundingBox << m_volume;

    //qDebug() << error_project/error << errorBoundingBox/error << error_negative/error << errorSkin/error;

    error *= (m_grammar->m_nrIndependentVariables/10.0);



    return 1000.0*error / m_samples.size();

}

bool CMAESWrapper::is_feasible(double const *v)
{
    bool result = false;

    if (*v > 0.0 && *v < 10.0)
    {
        result = true;
    }

    return result;
}

// Main place for optimization
vector<double> CMAESWrapper::optimize()
{
    const clock_t begin_time = clock();
    CMAES<double> evo;
    double *arFunvals, *const*pop, *xtemp;
    vector<vector<double>> pops;
    // Initialize everything
    const int dim = m_dimension, iter = 1;
    vector<double> xbest;
    xbest.resize(dim);
    all_error_vals.clear();
    
    all_error_vals.resize(iter);
    int bestIndex = 0;

    //vector<double> dimcheck;
    //for (int i = 0; i < dim; i++)
    //{
    //	//dimcheck.push_back(rand(4.0f, 12.0f));
    //	dimcheck.push_back(5.0f);
    //	std::cout << dimcheck[i]<< std::endl;
    //}

    //for (int i = 0; i < dim/2; i++)
    //{
    //	dimcheck.push_back(5.0f);
    //	dimcheck.push_back(1.0f);
    //	std::cout << dimcheck[2*i] <<" "<< dimcheck[2*i+1]<< std::endl;
    //}

    double bestError = 10000.0;

    ofstream fout;
    fout.open("test/Optimization Error/optimization_error.txt");

    for (int i = 0; i < iter; ++i)
    {
        local_error.clear();

        vector<double> xstart;
        xstart.resize(dim);

        m_grammar->update();
        m_grammar->softCleanUp();
        int set_Val = 0;
        std::cout << "Set custom initialization mode:";
        std::cin >> set_Val;
        if (set_Val == 1)
            m_scene->obtain_gram();
        else if (set_Val == 0)
            m_grammar->derive();
        else
        {
            vector<double> r;
            r.resize(m_grammar->m_nrIndependentVariables);
            for (int nums = 0; nums < r.size(); ++nums)
                r[nums] = rand(0.0, 1.0);

            m_grammar->adjustParameters(r.data(), r.size());
        }
    
        xstart = m_grammar->setParameterVals();
        
        //for (int i = 0; i < dim; i++)
        //    xstart[i] = 0.2f; // params::inst()->optMean;

        vector<double> stddev;
        stddev.resize(dim);

        for (int i = 0; i < dim; i++)
            stddev[i] = 0.005f; // params::inst()->optStd;

        Parameters<double> parameters;
        double avgError = 0.0;
        // TODO Adjust parameters here
        parameters.init(dim, xstart.data(), stddev.data());
        arFunvals = evo.init(parameters);

        std::cout << evo.sayHello() << std::endl;

        // Iterate until stop criterion holds
        while (!evo.testForTermination())
        {
            // Generate lambda new search points, sample population
            pop = evo.samplePopulation(); // Do not change content of pop

            //for (int i = 0; i < 17; i++)
            //{
            //	pops.push_back(dimcheck);
            //}

            /* Here you may resample each solution point pop[i] until it
            becomes feasible, e.g. for box constraints (variable
            boundaries). function is_feasible(...) needs to be
            user-defined.
            Assumptions: the feasible domain is convex, the optimum is
            not on (or very close to) the domain boundary, initialX is
            feasible and initialStandardDeviations are sufficiently small
            to prevent quasi-infinite looping.
            */

            /*for (int i = 0; i < evo.get(CMAES<double>::PopSize); ++i)
            {
            while (!is_feasible(pop[i]))
            evo.reSampleSingle(i);
            }*/

            // evaluate the new search points using fitfun from above

            for (int i = 0; i < evo.get(CMAES<double>::Lambda); ++i)
            {
                arFunvals[i] = fitfun(pop[i]);
                //arFunvals[i] = fitfun(pops[i], (int)evo.get(CMAES<double>::Dimension), dimcheck);
                avgError += arFunvals[i];
            }
            avgError /= evo.get(CMAES<double>::Lambda)*1.0;
            const double *xcurrent;
            xcurrent = evo.getNew(CMAES<double>::XBestEver);
            QString solutionString;
            for (int i = 0; i < m_dimension; ++i)
            {
                solutionString.append(QString::number(xcurrent[i], 'f', 2));
                if (i < m_dimension - 1)
                {
                    solutionString.append(", ");
                }
            }

            //qDebug() << avgError << qPrintable(solutionString);

            if (m_scene) m_scene->updateGrammarGeometry_CMAES();


            // update the search distribution used for sampleDistribution()
            evo.updateDistribution(arFunvals);

            //fout << avgError << endl;
        }

        std::cout << "Stop:" << std::endl << evo.getStopMessage();
        //evo.writeToFile(CMAES<double>::WKResume, "resumeevo1.dat"); // write resumable state of CMA-ES

        // get best estimator for the optimum, xmean
        xtemp = evo.getNew(CMAES<double>::XMean); // "XBestEver" might be used as well

        if (bestError > avgError)
        {
            bestError = avgError;
            bestIndex = i;
            for (int j = 0; j < dim; ++j)
                xbest[j] = xtemp[j];
        }

        all_error_vals[i] = local_error;
    }

    fout.close();

    for (int i = 0; i < dim; i++)
        qDebug() << xbest[i];

    // do something with final solution and finally release memory

    const double *xfinal = xbest.data();

    m_grammar->softCleanUp();
    m_grammar->adjustParameters(xfinal, m_grammar->m_nrIndependentVariables);
    m_grammar->derive();

    qDebug() << "Best Error:" << bestError;
    local_error = all_error_vals[bestIndex];
    m_besterror = bestError;
    std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC;
    return xbest;
    //delete[] xfinal;
}

void CMAESWrapper::setSampleVals(const vector<vec3> &samples)
{
    m_samples = samples;
}

void CMAESWrapper::setVolume(float volume)
{
    m_volume = volume;
}

void CMAESWrapper::setGrid(vector<vector<vector<int>>> grid)
{
    m_grid = grid;
}

void CMAESWrapper::setGridBB(vector<vec3> gridbb)
{
    m_gridbb = gridbb;
}

double CMAESWrapper::errorBB(vector<vec3> &bbox, Proc_Model &p, vector<vec3> &boxEnds)
{
    double error_bbox_vols = 0.0;

    vector<vec3> bbox_trial = p.bounding_box(boxEnds);
    vec3 outerbox = bbox[1], innerbox = bbox[0], outerbox_trial = bbox_trial[1], innerbox_trial = bbox_trial[0];
    double volume_bbox = (outerbox.x - innerbox.x)*(outerbox.y - innerbox.y)*(outerbox.z - innerbox.z);
    double volume_bbox_trial = (outerbox_trial.x - innerbox_trial.x)*(outerbox_trial.y - innerbox_trial.y)*(outerbox_trial.z - innerbox_trial.z);
    //error_bbox_vols = abs(volume_bbox - volume_bbox_trial);
    double common_volume = 0.0, bigger_volume = 0.0;
    vector<vec3> common_box, bigger_box;
    common_box.resize(2);
    bigger_box.resize(2);
    common_box[0] = vec3(max(innerbox.x, innerbox_trial.x),max(innerbox.y, innerbox_trial.y),max(innerbox.z, innerbox_trial.z));
    common_box[1] = vec3(min(outerbox.x, outerbox_trial.x),min(outerbox.y, outerbox_trial.y),min(outerbox.z, outerbox_trial.z));
    common_volume = (max(0.0f,common_box[1].x-common_box[0].x))*(max(0.0f, common_box[1].y-common_box[0].y))*(max(0.0f, common_box[1].z - common_box[0].z));
    bigger_box[0] = vec3(min(innerbox.x, innerbox_trial.x), min(innerbox.y, innerbox_trial.y), min(innerbox.z, innerbox_trial.z));
    bigger_box[1] = vec3(max(outerbox.x, outerbox_trial.x), max(outerbox.y, outerbox_trial.y), max(outerbox.z, outerbox_trial.z));
    bigger_volume = (bigger_box[1].x - bigger_box[0].x)*(bigger_box[1].y - bigger_box[0].y)*(bigger_box[1].z - bigger_box[0].z);
    //error_bbox_vols = volume_bbox_trial - common_volume;
    error_bbox_vols = bigger_volume - common_volume;
    return error_bbox_vols;
}

double CMAESWrapper::errorGrammar(vector<vector<int>> &v, vector<vec3> &boxEnds)
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

double CMAESWrapper::errorDisentangle(Proc_Model &p)
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

double CMAESWrapper::errorSkinning(Proc_Model &p)
{
    double skinerror = 0.0;

    vector<vec3> proc_pos = p.position(), proc_size = p.size();
    for (int i = 0; i < proc_pos.size(); i++)
    {
        vec3 size_i = proc_size[i];
        skinerror += 0.01*dot(size_i, size_i);
        skinerror += abs(size_i.x*size_i.y*size_i.z);
        skinerror += 0.1*abs(size_i.x*size_i.y + size_i.y*size_i.z + size_i.z*size_i.x);
        //if (params::inst()->isPartial)
        //    skinerror += 0.1*(abs(size_i.x) + abs(size_i.y) + abs(size_i.z));
    }

    return 100*skinerror;
}

double CMAESWrapper::errorEmpty()
{
    double error = 0.0;

    return error;
}

double CMAESWrapper::errorPenalty(Proc_Model &p)
{
    double penalty = 0.0;
    vector<vec3> proc_pos = p.position(), proc_size = p.size();
    for (int i = 0; i < proc_pos.size(); i++)
    {
        vec3 size_i = proc_size[i];
        double error_val = abs(size_i.x*size_i.y*size_i.z);
        if (error_val > 0)
            penalty += (0.001 / error_val);
        else penalty += 1000;
    }
    return penalty;
}

vec3 CMAESWrapper::conventionPoint(vec3 &start, vec3 &end, int num)
{
    vec3 v;
    switch (num){
    case 1: v = start;                            break;
    case 2: v = vec3(end.x, start.y, start.z);	  break;
    case 3: v = vec3(start.x, end.y, start.z);	  break;
    case 4: v = vec3(end.x, end.y, start.z);	  break;
    case 5: v = vec3(start.x, start.y, end.z);	  break;
    case 6: v = vec3(end.x, start.y, end.z);	  break;
    case 7: v = vec3(start.x, end.y, end.z);	  break;
    case 8: v = end;							  break;
    }

    return v;
}

void CMAESWrapper::setGrammar(Grammar *gram)
{
    m_grammar = gram;
}

double CMAESWrapper::dimension()
{
    return m_dimension;
}
