#include "Proc_Model.h"

Proc_Model::Proc_Model()
{

}

Proc_Model::~Proc_Model()
{

}

vector<vec3> &Proc_Model::position()
{
	return m_pos;
}

vector<vec3> &Proc_Model::size()
{
	return m_size;
}

vector<double> &Proc_Model::df()
{
	return m_df;
}


vector<vector<double>> &Proc_Model::grammar()
{
	return m_grammar;
}

void Proc_Model::new_proc(const vector<vec3> &positions, const vector<vec3> &sizes)
{
	m_pos = positions;
	m_size = sizes;

	m_range = 16.0;
	m_tot = 64;
	m_numsamples_per_face = 100;
	m_df.resize(m_tot*m_tot*m_tot);
}

void Proc_Model::new_proc(int num_boxes)
{
	srand(time(NULL));
	std::default_random_engine unif_generator;
	std::uniform_real_distribution<double> unif_distribution(0.0, 10.0), unif_distribution1(0.0, 6.0);

	m_range = 16.0;
	m_tot = 32;
	m_numsamples_per_face = 100;
	m_df.resize(m_tot*m_tot*m_tot);

	for (int i = 0; i < num_boxes; i++)
	{
		double posx = unif_distribution(unif_generator);
		double posy = unif_distribution(unif_generator);
		double posz = unif_distribution(unif_generator);
		double sizex = unif_distribution1(unif_generator);
		double sizey = unif_distribution1(unif_generator);
		double sizez = unif_distribution1(unif_generator);
		m_pos.push_back(vec3(posx, posy, posz));
		m_size.push_back(vec3(sizex, sizey, sizez));
	}
}

vector<double> Proc_Model::create_distance_field(vector<vec3> size_vec, vector<vec3> pos_vec)
{
	vector<double> df;
	vector<vec3> size, pos; 
	vector<vec3> start_up, start_down, end_up, end_down;
	const int tot = m_tot;
	vector<vector<double>> neigs;
	double range = m_range;
	int num_samples_per_face = m_numsamples_per_face;

	std::default_random_engine unif_generator;
	std::uniform_real_distribution<double> unif_distribution(0.0, 1.0);

	
	for (int i = 0; i < size_vec.size(); i++)
	{
		vector<vec3> samples;
		vec3 size = size_vec[i];
		vec3 pos = pos_vec[i];

		for (int j = 0; j < num_samples_per_face; j++)
		{
			double d1 = unif_distribution(unif_generator);
			double d2 = unif_distribution(unif_generator);
			vec3 samples_new;

			//Bottom face
			samples_new.z = pos.z;
			samples_new.x = pos.x + d1*size.x;
			samples_new.y = pos.y + d2*size.y;
			samples.push_back(samples_new);

			//Top face
			samples_new.z += size.z;
			samples.push_back(samples_new);

			//Side face x-
			samples_new.x = pos.x;
			samples_new.y = pos.y + d1*size.y;
			samples_new.z = pos.z + d2*size.z;
			samples.push_back(samples_new);

			//Side face x+
			samples_new.x += size.x;
			samples.push_back(samples_new);

			//Side face y-
			samples_new.y = pos.y;
			samples_new.x = pos.x + d2*size.x;
			samples_new.z = pos.z + d1*size.z;
			samples.push_back(samples_new);

			//Side face y+
			samples_new.y += size.y;
			samples.push_back(samples_new);
		}

		for (int xval = 0; xval < tot; xval++)
		{
			for (int yval = 0; yval < tot; yval++)
			{
				for (int zval = 0; zval < tot; zval++)
				{
					vec3 point_val;
					point_val.x = xval*range / (tot - 1.0f);
					point_val.y = yval*range / (tot - 1.0f);
					point_val.z = zval*range / (tot - 1.0f);

					double min_val = tot * 2;
					for (int i = 0; i < samples.size(); i++)
					{
						vec3 diff = point_val - samples[i];
						double d = sqrt(dot(diff, diff));
						if (min_val > d)
							min_val = d;
					}
					if (i == 0)
						df.push_back(min_val);
					else
					{
						int iden = zval + yval*tot + xval*tot*tot;
						df[iden] = min(df[iden], min_val);
					}
				}
			}
		}
	}
	/*for (int i = 0; i < m_size.size(); i++)
	{
		size.push_back(m_size[i]*(tot/range));
		pos.push_back(m_pos[i] * (tot / range));
		start_up.push_back(vec3(ceil(size[i].x), ceil(size[i].y), ceil(size[i].z)));
		start_down.push_back(vec3(floor(size[i].x), floor(size[i].y), floor(size[i].z)));
		end_up.push_back(vec3(ceil(size[i].x + pos[i].x), ceil(size[i].y + pos[i].y), ceil(size[i].z + pos[i].z)));
		end_down.push_back(vec3(floor(size[i].x + pos[i].x), floor(size[i].y + pos[i].y), floor(size[i].z + pos[i].z)));
	}*/

	/*for (int xval = 0; xval < tot; xval++)
	{
		vector<double> xshot;
		xshot.push_back(xval);
		if (xval == 0)
			xshot.push_back(xval + 1);
		else if (xval == tot - 1)
			xshot.push_back(xval - 1);
		else
		{
			xshot.push_back(xval - 1);
			xshot.push_back(xval + 1);
		}

		for (int yval = 0; yval < tot; yval++)
		{
			vector<double> yshot;
			yshot.push_back(yval);
			if (yval == 0)
				yshot.push_back(yval + 1);
			else if (yval == tot - 1)
				yshot.push_back(yval - 1);
			else
			{
				yshot.push_back(yval - 1);
				yshot.push_back(yval + 1);
			}

			for (int zval = 0; zval < tot; zval++)
			{
				vector<double> zshot;
				zshot.push_back(zval);
				if (zval == 0)
					zshot.push_back(zval + 1);
				else if (xval == tot - 1)
					zshot.push_back(zval - 1);
				else
				{
					zshot.push_back(zval - 1);
					zshot.push_back(zval + 1);
				}

				vector<double> spec_neigs;

				for (int i = 0; i < xshot.size(); i++)
				{
					int t1 = xshot[i];
					for (int j = 0; j < yshot.size(); j++)
					{
						int t2 = yshot[j];
						for (int k = 0; k < zshot.size(); k++)
						{
							int t3 = zshot[k];
							if (i + j + k >0)
							{
								spec_neigs.push_back(t1*tot*tot + t2*tot + t3);
							}
						}
					}
				}
				neigs.push_back(spec_neigs);
			}
		}
	}*/

	return df;
}

void Proc_Model::create_box_distance_field(const vector<vec3> &size_vec, const vector<vec3> &pos_vec, vector<double>& df)
{
	double eval = 0.0;
	int tot = m_tot;
	int range = m_range;

	for (int i = 0; i < pos_vec.size(); i++)
	{
		vec3 a = pos_vec[i];
		vec3 b = pos_vec[i] + size_vec[i];

		for (int xval = 0; xval < tot; xval++)
			{
				for (int yval = 0; yval < tot; yval++)
				{
					for (int zval = 0; zval < tot; zval++)
					{
						double val = 10000;
						vec3 point_val;
						point_val.x = xval*range / (tot - 1.0f);
						point_val.y = yval*range / (tot - 1.0f);
						point_val.z = zval*range / (tot - 1.0f);

						double x1, y1, z1;

						z1 = min(abs(point_val.z - a.z), abs(point_val.z - b.z));
						if ((a.x <= point_val.x) && (b.x >= point_val.x))
							x1 = 0;
						else x1 = min(abs(a.x - point_val.x), abs(b.x - point_val.x));
						if ((a.y <= point_val.y) && (b.y >= point_val.y))
							y1 = 0;
						else y1 = min(abs(a.y - point_val.y), abs(b.y - point_val.y));

						val = min(val, sqrt(x1*x1 + y1*y1 + z1*z1));

						x1 = min(abs(point_val.x - a.x), abs(point_val.x - b.x));
						if ((a.z <= point_val.z) && (b.z >= point_val.z))
							z1 = 0;
						else z1 = min(abs(a.z - point_val.z), abs(b.z - point_val.z));
						if ((a.y <= point_val.y) && (b.y >= point_val.y))
							y1 = 0;
						else y1 = min(abs(a.y - point_val.y), abs(b.y - point_val.y));

						val = min(val, sqrt(x1*x1 + y1*y1 + z1*z1));

						y1 = min(abs(point_val.y - a.y), abs(point_val.y - b.y));
						if ((a.z <= point_val.z) && (b.z >= point_val.z))
							z1 = 0;
						else z1 = min(abs(a.z - point_val.z), abs(b.z - point_val.z));
						if ((a.x <= point_val.x) && (b.x >= point_val.x))
							x1 = 0;
						else x1 = min(abs(a.x - point_val.x), abs(b.x - point_val.x));

						val = min(val, sqrt(x1*x1 + y1*y1 + z1*z1));
						
						if (i == 0)
							df[zval + (yval + xval*tot)*tot] = sqrt(x1*x1 + y1*y1 + z1*z1);
						else df[zval + (yval + xval*tot)*tot] = min(df[zval + (yval + xval*tot)*tot], val);
					}
				}
			}
	}
}

vector<double> Proc_Model::create_point_cloud_distance_field(vector<vec3> samples)
{
	vector<double> df;// = create_distance_field(size, instance);
	double eval = 0.0;
	int tot = m_tot;
	df.resize(tot*tot*tot);
	int range = m_range;

	for (int xval = 0; xval < tot; xval++)
	{
		for (int yval = 0; yval < tot; yval++)
		{
			for (int zval = 0; zval < tot; zval++)
			{
				float val;
				vec3 point_val;
				point_val.x = xval*range / (tot - 1.0f);
				point_val.y = yval*range / (tot - 1.0f);
				point_val.z = zval*range / (tot - 1.0f);

				for (int i = 0; i < samples.size(); i++)
				{
					vec3 sample_point = samples[i];
					vec3 point_diff = sample_point - point_val;
					if (i == 0)
						val = dot(point_diff, point_diff);
					else
					{
						val = min(val, dot(point_diff, point_diff));
					}
				}

				df[zval + (yval + xval*tot)*tot] = val;
			}
		}
	}

	return df;

}

void Proc_Model::create_default_distance_field()
{
	create_box_distance_field(m_size, m_pos, m_df);
}

void Proc_Model::create_FLANN_default_distance_field()
{
	//int num_points = 128;
	//const int tot = m_tot;
	//int overall_tot = tot*tot*tot;
	//int overall_numpts = num_points*num_points*6*m_pos.size();

	//vector<float> sourceData;
	//sourceData.resize(overall_numpts*3);
	//vector<float> queryData;
	//queryData.resize(overall_tot * 3);

	//int range = m_range;
	//int current = 0;

	//double xpoint, ypoint, zpoint;

	//for (int xval = 0; xval < tot; xval++)
	//{
	//	xpoint = xval*range / (tot - 1.0f);
	//	for (int yval = 0; yval < tot; yval++)
	//	{
	//		ypoint = yval*range / (tot - 1.0f);
	//		for (int zval = 0; zval < tot; zval++)
	//		{
	//			zpoint = zval*range / (tot - 1.0f);
	//			queryData[current] = xpoint;
	//			++current;
	//			queryData[current] = ypoint;
	//			++current;
	//			queryData[current] = zpoint;
	//			++current;
	//		}
	//	}
	//}

	//sample_shape_surface(num_points, sourceData.data());

	//flann::Matrix<float> dataset(sourceData.data(), sourceData.size()/3, 3, sizeof(float) * 3);
	//flann::Matrix<float> query(queryData.data(), queryData.size()/3, 3, sizeof(float) * 3);

	//flann::Matrix<int> indices(new int[query.rows*3], query.rows, 3);
	//flann::Matrix<float> dists(new float[query.rows*3], query.rows, 3);

	////flann::Index<flann::L2<float> > index(dataset, flann::KDTreeIndexParams(1));
	//flann::KDTreeSingleIndex<flann::L2<float> > index(dataset, flann::KDTreeSingleIndexParams(16));
	//index.buildIndex();

	//index.knnSearch(query, indices, dists, 3, flann::SearchParams(16));

	//for (int i = 0; i < tot*tot*tot; i++)
	//	m_df[i] = dists[i][0];
}

double Proc_Model::compare_df(vector<vec3> instance, vector<vec3> size)
{
	vector<double> df;// = create_distance_field(size, instance);
	double eval = 0.0;
	int tot = m_tot;
	int range = m_range;
	
	for (int i = 0; i < instance.size(); i++)
	{
		vec3 a = instance[i];
		vec3 b = instance[i] + size[i];

		if ((a.x < 0) || (a.y <0) || (a.z < 0) || (a.x > range) || (a.y >range) || (a.z > range) || (b.x < 0) || (b.y < 0) || (b.z < 0) || (b.x > range) || (b.y > range) || (b.z > range) || (size[i].x < 0) || (size[i].y < 0) || (size[i].z < 0))
			return 10000;
	}
	
	create_box_distance_field(size, instance, df);

	for (int i = 0; i < m_df.size(); i++)
		eval += (df[i] - m_df[i])*(df[i] - m_df[i]);
	
	eval /= m_df.size()*1.0;
	
	return eval;
}

double Proc_Model::evaluate_df_val(vector<vec3> instance, vector<vec3> size)
{
	int num_boxes = m_pos.size();
	double eval = 0.0, lambda_vol = 0.0, lambda_sa = 0.0;
	for (int i = 0; i < num_boxes; i++)
	{
		eval += Df_Cube(instance[i].x, instance[i].x+size[i].x, instance[i].y, instance[i].y+size[i].y, instance[i].z, instance[i].z+size[i].z);
		eval += lambda_vol*abs(m_size[i].x*m_size[i].y*m_size[i].z - size[i].x*size[i].y*size[i].z);
		eval += 2 * lambda_sa*((m_size[i].x*m_size[i].y + m_size[i].y*m_size[i].z + m_size[i].z*m_size[i].x) - (size[i].x*size[i].y+size[i].y*size[i].z+size[i].z*size[i].x));
	}
	return eval;
}

double Proc_Model::project_point_to_model(vec3 point)
{
	double val = 100000; 

	for (int i = 0; i < m_pos.size(); i++)
	{
		vec3 a = m_pos[i], size = m_size[i];
		vec3 b = a + size;

		double x1, y1, z1;

		z1 = min(abs(point.z-a.z), abs(point.z-b.z));
		if ((a.x <= point.x) && (b.x >= point.x))
			x1 = 0;
		else x1 = min(abs(a.x - point.x), abs(b.x - point.x));
		if ((a.y <= point.y) && (b.y >= point.y))
			y1 = 0;
		else y1 = min(abs(a.y - point.y), abs(b.y - point.y));

		val = min(val, sqrt(x1*x1 + y1*y1 + z1*z1));

		x1 = min(abs(point.x - a.x), abs(point.x - b.x));
		if ((a.z <= point.z) && (b.z >= point.z))
			z1 = 0;
		else z1 = min(abs(a.z - point.z), abs(b.z - point.z));
		if ((a.y <= point.y) && (b.y >= point.y))
			y1 = 0;
		else y1 = min(abs(a.y - point.y), abs(b.y - point.y));

		val = min(val, sqrt(x1*x1 + y1*y1 + z1*z1));

		y1 = min(abs(point.y - a.y), abs(point.y - b.y));
		if ((a.x <= point.x) && (b.x >= point.x))
			x1 = 0;
		else x1 = min(abs(a.x - point.x), abs(b.x - point.x));
		if ((a.z <= point.z) && (b.z >= point.z))
			z1 = 0;
		else z1 = min(abs(a.z - point.z), abs(b.z - point.z));

		val = min(val, sqrt(x1*x1 + y1*y1 + z1*z1));
	}

	return val;
}

double Proc_Model::project_point_to_box(vec3 point, vec3 position, vec3 size)
{
    vec3 a = position;
    vec3 b = a + size;

    double x1, y1, z1, val = 100000;

    z1 = min(abs(point.z - a.z), abs(point.z - b.z));
    if ((a.x <= point.x) && (b.x >= point.x))
        x1 = 0;
    else x1 = min(abs(a.x - point.x), abs(b.x - point.x));
    if ((a.y <= point.y) && (b.y >= point.y))
        y1 = 0;
    else y1 = min(abs(a.y - point.y), abs(b.y - point.y));

    val = min(val, sqrt(x1*x1 + y1*y1 + z1*z1));

    x1 = min(abs(point.x - a.x), abs(point.x - b.x));
    if ((a.z <= point.z) && (b.z >= point.z))
        z1 = 0;
    else z1 = min(abs(a.z - point.z), abs(b.z - point.z));
    if ((a.y <= point.y) && (b.y >= point.y))
        y1 = 0;
    else y1 = min(abs(a.y - point.y), abs(b.y - point.y));

    val = min(val, sqrt(x1*x1 + y1*y1 + z1*z1));

    y1 = min(abs(point.y - a.y), abs(point.y - b.y));
    if ((a.x <= point.x) && (b.x >= point.x))
        x1 = 0;
    else x1 = min(abs(a.x - point.x), abs(b.x - point.x));
    if ((a.z <= point.z) && (b.z >= point.z))
        z1 = 0;
    else z1 = min(abs(a.z - point.z), abs(b.z - point.z));

    val = min(val, sqrt(x1*x1 + y1*y1 + z1*z1));

    return val;
}

double Proc_Model::project_model_to_pc(vector<vec3> points)
{
    int num_points = 16;
    int overall_numpts = num_points*num_points * 6 * m_pos.size();

    vector<float> queryData;
    queryData.resize(overall_numpts * 3);
    sample_shape_surface(num_points, queryData.data());

    vector<float> sourceData;
    sourceData.resize(3 * points.size());
    for (int i = 0; i < points.size(); ++i)
    {
        sourceData[3 * i] = points[i].x;
        sourceData[3 * i + 1] = points[i].y;
        sourceData[3 * i + 2] = points[i].z;
    }

    flann::Matrix<float> dataset(sourceData.data(), sourceData.size() / 3, 3, sizeof(float) * 3);
    flann::Matrix<float> query(queryData.data(), queryData.size() / 3, 3, sizeof(float) * 3);

    flann::Matrix<int> indices(new int[query.rows], query.rows, 1);
    flann::Matrix<float> dists(new float[query.rows], query.rows, 1);

    flann::KDTreeSingleIndex<flann::L2<float> > index(dataset, flann::KDTreeSingleIndexParams(16));
    index.buildIndex();

    index.knnSearch(query, indices, dists, 1, flann::SearchParams(128));

    double project_error = 0.0;
    for (int i = 0; i < queryData.size() / 3; ++i)
        project_error += dists[i][0];

    delete[] indices.ptr();
    delete[] dists.ptr();

    return project_error / overall_numpts;
}

void Proc_Model::sample_shape_surface(int num_points, float samples[])
{
	int current = 0;
	for (int i = 0; i < m_pos.size(); i++)
	{

		float xpoint, ypoint, zpoint;

		//Positive x-side
		xpoint = m_pos[i].x + m_size[i].x;
		for (int yval = 0; yval < num_points; yval++)
		{
			ypoint = m_pos[i].y + yval*m_size[i].y / ((num_points - 1)*1.0f);
			for (int zval = 0; zval < num_points; zval++)
			{
				zpoint = m_pos[i].z + zval*m_size[i].z / ((num_points - 1)*1.0f);
				samples[current] = xpoint;
				++current;
				samples[current] = ypoint;
				++current;
				samples[current] = zpoint;
				++current;
			}
		}

		//Negative x-side
		xpoint = m_pos[i].x;
		for (int yval = 0; yval < num_points; yval++)
		{
			ypoint = m_pos[i].y + yval*m_size[i].y / ((num_points - 1)*1.0f);
			for (int zval = 0; zval < num_points; zval++)
			{
				zpoint = m_pos[i].z + zval*m_size[i].z / ((num_points - 1)*1.0f);
				samples[current] = xpoint;
				++current;
				samples[current] = ypoint;
				++current;
				samples[current] = zpoint;
				++current;
			}
		}

		//Positive y-side
		ypoint = m_pos[i].y + m_size[i].y;
		for (int xval = 0; xval < num_points; xval++)
		{
			xpoint = m_pos[i].x + xval*m_size[i].x / ((num_points - 1)*1.0f);
			for (int zval = 0; zval < num_points; zval++)
			{
				zpoint = m_pos[i].z + zval*m_size[i].z / ((num_points - 1)*1.0f);
				samples[current] = xpoint;
				++current;
				samples[current] = ypoint;
				++current;
				samples[current] = zpoint;
				++current;
			}
		}


		//Negative y-side
		ypoint = m_pos[i].y;
		for (int xval = 0; xval < num_points; xval++)
		{
			xpoint = m_pos[i].x + xval*m_size[i].x / ((num_points - 1)*1.0f);
			for (int zval = 0; zval < num_points; zval++)
			{
				zpoint = m_pos[i].z + zval*m_size[i].z / ((num_points - 1)*1.0f);
				samples[current] = xpoint;
				++current;
				samples[current] = ypoint;
				++current;
				samples[current] = zpoint;
				++current;
			}
		}

		//Positive z-side
		zpoint = m_pos[i].z + m_size[i].z;
		for (int xval = 0; xval < num_points; xval++)
		{
			xpoint = m_pos[i].x + xval*m_size[i].x / ((num_points - 1)*1.0f);
			for (int yval = 0; yval < num_points; yval++)
			{
				ypoint = m_pos[i].y + yval*m_size[i].y / ((num_points - 1)*1.0f);
				samples[current] = xpoint;
				++current;
				samples[current] = ypoint;
				++current;
				samples[current] = zpoint;
				++current;
			}
		}

		//Negative z-side
		zpoint = m_pos[i].z;
		for (int xval = 0; xval < num_points; xval++)
		{
			xpoint = m_pos[i].x + xval*m_size[i].x / ((num_points - 1)*1.0f);
			for (int yval = 0; yval < num_points; yval++)
			{
				ypoint = m_pos[i].y + yval*m_size[i].y / ((num_points - 1)*1.0f);
				samples[current] = xpoint;
				++current;
				samples[current] = ypoint;
				++current;
				samples[current] = zpoint;
				++current;
			}
		}
	}
}

vector<vec3> Proc_Model::bounding_box()
{
    float min_x, max_x, min_y, max_y, min_z, max_z;
    vector<vec3> samples;
    for (int i = 0; i < m_pos.size(); ++i)
    {
        samples.push_back(m_pos[i]);
        samples.push_back(m_pos[i] + m_size[i]);
    }

    return bounding_box(samples);
}

vector<vec3> Proc_Model::bounding_box(vector<Vertex> &samples)
{
    vector<vec3> m_samples;
    m_samples.resize(samples.size());
    for (int i = 0; i < samples.size(); ++i)
        m_samples[i] = samples[i].position;

    return bounding_box(m_samples);
}

vector<vec3> Proc_Model::bounding_box(vector<vec3> &samples)
{
	float min_x, max_x, min_y, max_y, min_z, max_z;
	for (int i = 0; i < samples.size(); i++)
	{
		vec3 p = samples[i];
		if (i == 0)
		{
			min_x = p.x;
			max_x = min_x;
			min_y = p.y;
			max_y = min_y;
			min_z = p.z;
			max_z = min_z;
		}
		else
		{
			min_x = min(min_x, p.x);
			min_y = min(min_y, p.y);
			min_z = min(min_z, p.z);
			max_x = max(max_x, p.x);
			max_y = max(max_y, p.y);
			max_z = max(max_z, p.z);
		}
	}
	vector<vec3> bbox;
	bbox.push_back(vec3(min_x, min_y, min_z));
	bbox.push_back(vec3(max_x, max_y, max_z));

	return bbox;
}

int Proc_Model::optimization(vector<double> solvals)
{
	int minval = 0;
	for (int i = 0; i < solvals.size(); i++)
	{
		if (solvals[i] < solvals[minval])
			minval = i;
	}

	return minval;
}

double Proc_Model::Df_Cube(double x1, double x2, double y1, double y2, double z1, double z2)
{
	int tot = 10;
	int max_grid = 128, max_val = 16;
	vector<double> data = m_df;
	
	x1 *= ((max_grid-1)/max_val);
	x2 *= ((max_grid - 1) / max_val);
	y1 *= ((max_grid - 1) / max_val);
	y2 *= ((max_grid - 1) / max_val);
	z1 *= ((max_grid - 1) / max_val);
	z2 *= ((max_grid - 1) / max_val);

	float stepx = (x2 - x1)*1.0f / tot;
	float stepy = (y2 - y1)*1.0f / tot;
	float stepz = (z2 - z1)*1.0f / tot;
	float df_val = 0;
	int xplus, xminus, yplus, yminus, zplus, zminus;
	float z = z1;
	zminus = floor(z);
	zplus = floor(z) + 1;

	if (z < 0)
		z = 0;
	if (z > (max_grid-1))
		z = (max_grid - 1);
	if (zminus < 0)
		zminus = 0;
	if (zplus < 0)
		zplus = 0;
	if (zminus >(max_grid - 1))
		zminus = (max_grid - 1);
	if (zplus >(max_grid - 1))
		zplus = (max_grid - 1);
	for (int i = 0; i < tot + 1; i++)
	{
		float x = x1 + i*stepx;
		xminus = floor(x);
		xplus = floor(x) + 1;
		if (x < 0)
			x = 0;
		if (x > 127)
			x = 127;
		if (xminus < 0)
			xminus = 0;
		if (xplus < 0)
			xplus = 0;
		if (xminus >(max_grid - 1))
			xminus = (max_grid - 1);
		if (xplus >(max_grid - 1))
			xplus = (max_grid - 1);
		for (int j = 0; j < tot + 1; j++)
		{
			float y = y1 + j*stepy;
			yminus = floor(y);
			yplus = floor(y) + 1;
			if (y < 0)
				y = 0;
			if (y > (max_grid-1))
				y = (max_grid-1);
			if (yminus < 0)
				yminus = 0;
			if (yplus < 0)
				yplus = 0;
			if (yminus > (max_grid-1))
				yminus = (max_grid-1);
			if (yplus > (max_grid-1))
				yplus = (max_grid-1);
			df_val += trilinear_interp(x, y, z, xminus, xplus, yminus, yplus, zminus, zplus);
		}
	}
	z = z2;
	zminus = floor(z);
	zplus = floor(z) + 1;
	if (z < 0)
		z = 0;
	if (z > (max_grid-1))
		z = (max_grid-1);
	if (zminus < 0)
		zminus = 0;
	if (zplus < 0)
		zplus = 0;
	if (zminus > (max_grid-1))
		zminus = (max_grid-1);
	if (zplus > (max_grid-1))
		zplus = (max_grid-1);
	for (int i = 0; i < tot + 1; i++)
	{
		float x = x1 + i*stepx;
		xminus = floor(x);
		xplus = floor(x) + 1;
		if (x < 0)
			x = 0;
		if (x > (max_grid-1))
			x = (max_grid-1);
		if (xminus < 0)
			xminus = 0;
		if (xplus < 0)
			xplus = 0;
		if (xminus > (max_grid-1))
			xminus = (max_grid-1);
		if (xplus > (max_grid-1))
			xplus = (max_grid-1);
		for (int j = 0; j < tot + 1; j++)
		{
			float y = y1 + j*stepy;
			yminus = floor(y);
			yplus = floor(y) + 1;
			if (y < 0)
				y = 0;
			if (y > (max_grid-1))
				y = (max_grid-1);
			if (yminus < 0)
				yminus = 0;
			if (yplus < 0)
				yplus = 0;
			if (yminus > (max_grid-1))
				yminus = (max_grid-1);
			if (yplus > (max_grid-1))
				yplus = (max_grid-1);
			df_val += trilinear_interp(x, y, z, xminus, xplus, yminus, yplus, zminus, zplus);
		}
	}
	float x = x1;
	xminus = floor(x);
	xplus = floor(x) + 1;
	if (x < 0)
		x = 0;
	if (x > (max_grid-1))
		x = (max_grid-1);
	if (xminus < 0)
		xminus = 0;
	if (xplus < 0)
		xplus = 0;
	if (xminus > (max_grid-1))
		xminus = (max_grid-1);
	if (xplus > (max_grid-1))
		xplus = (max_grid-1);
	for (int i = 0; i < tot + 1; i++)
	{
		float z = z1 + i*stepz;
		zminus = floor(z);
		zplus = floor(z) + 1;
		if (z < 0)
			z = 0;
		if (z > (max_grid-1))
			z = (max_grid-1);
		if (zminus < 0)
			zminus = 0;
		if (zplus < 0)
			zplus = 0;
		if (zminus > (max_grid-1))
			zminus = (max_grid-1);
		if (zplus > (max_grid-1))
			zplus = (max_grid-1);
		for (int j = 0; j < tot + 1; j++)
		{
			float y = y1 + j*stepy;
			yminus = floor(y);
			yplus = floor(y) + 1;
			if (y < 0)
				y = 0;
			if (y > 127)
				y = 127;
			if (yminus < 0)
				yminus = 0;
			if (yplus < 0)
				yplus = 0;
			if (yminus > (max_grid-1))
				yminus = (max_grid-1);
			if (yplus > (max_grid-1))
				yplus = (max_grid-1);
			df_val += trilinear_interp(x, y, z, xminus, xplus, yminus, yplus, zminus, zplus);
		}
	}
	x = x2;
	xminus = floor(x);
	xplus = floor(x) + 1;
	if (x < 0)
		x = 0;
	if (x > (max_grid-1))
		x = (max_grid-1);
	if (xminus < 0)
		xminus = 0;
	if (xplus < 0)
		xplus = 0;
	if (xminus > (max_grid-1))
		xminus = (max_grid-1);
	if (xplus > (max_grid-1))
		xplus = (max_grid-1);
	for (int i = 0; i < tot + 1; i++)
	{
		float z = z1 + i*stepz;
		zminus = floor(z);
		zplus = floor(z) + 1;
		if (z < 0)
			z = 0;
		if (z > 127)
			z = 127;
		if (zminus < 0)
			zminus = 0;
		if (zplus < 0)
			zplus = 0;
		if (zminus > (max_grid-1))
			zminus = (max_grid-1);
		if (zplus > (max_grid-1))
			zplus = (max_grid-1);
		for (int j = 0; j < tot + 1; j++)
		{
			float y = y1 + j*stepy;
			yminus = floor(y);
			yplus = floor(y) + 1;
			if (y < 0)
				y = 0;
			if (y > 127)
				y = 127;
			if (yminus < 0)
				yminus = 0;
			if (yplus < 0)
				yplus = 0;
			if (yminus > (max_grid-1))
				yminus = (max_grid-1);
			if (yplus > (max_grid-1))
				yplus = (max_grid-1);
			df_val += trilinear_interp(x, y, z, xminus, xplus, yminus, yplus, zminus, zplus);
		}
	}
	float y = y1;
	yminus = floor(y);
	yplus = floor(y) + 1;
	if (y < 0)
		y = 0;
	if (y > (max_grid-1))
		y = (max_grid-1);
	if (yminus < 0)
		yminus = 0;
	if (yplus < 0)
		yplus = 0;
	if (yminus > (max_grid-1))
		yminus = (max_grid-1);
	if (yplus > (max_grid-1))
		yplus = (max_grid-1);
	for (int i = 0; i < tot + 1; i++)
	{
		float z = z1 + i*stepz;
		zminus = floor(z);
		zplus = floor(z) + 1;
		if (z < 0)
			z = 0;
		if (z > (max_grid-1))
			z = (max_grid-1);
		if (zminus < 0)
			zminus = 0;
		if (zplus < 0)
			zplus = 0;
		if (zminus > (max_grid-1))
			zminus = (max_grid-1);
		if (zplus > (max_grid-1))
			zplus = (max_grid-1);
		for (int j = 0; j < tot + 1; j++)
		{
			float x = x1 + j*stepx;
			xminus = floor(x);
			xplus = floor(x) + 1;
			if (x < 0)
				x = 0;
			if (x > (max_grid-1))
				x = (max_grid-1);
			if (xminus < 0)
				xminus = 0;
			if (xplus < 0)
				xplus = 0;
			if (xminus > (max_grid-1))
				xminus = (max_grid-1);
			if (xplus > (max_grid-1))
				xplus = (max_grid-1);
			df_val += trilinear_interp(x, y, z, xminus, xplus, yminus, yplus, zminus, zplus);
		}
	}
	y = y2;
	yminus = floor(y);
	yplus = floor(y) + 1;
	if (y < 0)
		y = 0;
	if (y > (max_grid-1))
		y = (max_grid-1);
	if (yminus < 0)
		yminus = 0;
	if (yplus < 0)
		yplus = 0;
	if (yminus > (max_grid-1))
		yminus = (max_grid-1);
	if (yplus > (max_grid-1))
		yplus = (max_grid-1);
	for (int i = 0; i < tot + 1; i++)
	{
		float z = z1 + i*stepz;
		zminus = floor(z);
		zplus = floor(z) + 1;
		if (z < 0)
			z = 0;
		if (z > (max_grid-1))
			z = (max_grid-1);
		if (zminus < 0)
			zminus = 0;
		if (zplus < 0)
			zplus = 0;
		if (zminus > (max_grid-1))
			zminus = (max_grid-1);
		if (zplus > (max_grid-1))
			zplus = (max_grid-1);
		for (int j = 0; j < tot + 1; j++)
		{
			float x = x1 + j*stepx;
			xminus = floor(x);
			xplus = floor(x) + 1;
			if (x < 0)
				x = 0;
			if (x > (max_grid-1))
				x = (max_grid-1);
			if (xminus < 0)
				xminus = 0;
			if (xplus < 0)
				xplus = 0;
			if (xminus > (max_grid-1))
				xminus = (max_grid-1);
			if (xplus > (max_grid-1))
				xplus = (max_grid-1);
			df_val += trilinear_interp(x, y, z, xminus, xplus, yminus, yplus, zminus, zplus);
		}
	}
	df_val /= (6.0f * pow(tot + 1, 2));
	return df_val;
}

double Proc_Model::trilinear_interp(double x, double y, double z, int x1, int x2, int y1, int y2, int z1, int z2)
{
	//vector<double> df = m_df;
	double p000 = m_df[x1 * 128 * 128 + y1 * 128 + z1];
	double p001 = m_df[x1 * 128 * 128 + y1 * 128 + z2];
	double p010 = m_df[x1 * 128 * 128 + y2 * 128 + z1];
	double p011 = m_df[x1 * 128 * 128 + y2 * 128 + z2];
	double p100 = m_df[x2 * 128 * 128 + y1 * 128 + z1];
	double p101 = m_df[x2 * 128 * 128 + y1 * 128 + z2];
	double p110 = m_df[x2 * 128 * 128 + y2 * 128 + z1];
	double p111 = m_df[x2 * 128 * 128 + y2 * 128 + z2];
	double al_x, al_y, al_z;
	(x2 > x1) ? (al_x = (x - x1) / ((x2 - x1) * 1.0)) : (al_x = 1);
	(y2 > y1) ? (al_y = (y - y1) / ((y2 - y1) * 1.0)) : (al_y = 1);
	(z2 > z1) ? (al_z = (z - z1) / ((z2 - z1) * 1.0)) : (al_z = 1);
	float val = (1.0 - al_x)*((1.0 - al_y)*((1.0 - al_z)*p000 + al_z*p001) + al_y*((1.0 - al_z)*p010 + al_z*p011)) + al_x*((1.0 - al_y)*((1.0 - al_z)*p100 + al_z*p101) + al_y*((1.0 - al_z)*p110 + al_z*p111));
	return val;
}

double Proc_Model::volume()
{
    double volume = 0.0;
    for (int i = 0; i < m_size.size(); ++i)
        volume += (m_size[i].x * m_size[i].y * m_size[i].z);

    return volume;
}

void Proc_Model::set_grammar(vector<vector<double>> grammar)
{
	m_grammar = grammar;
}
