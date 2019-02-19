#include "pointcloudio.h"
#include <iostream>
#include <fstream>
#include <QDebug>
using namespace std;

std::vector<vec3> PointCloudIO::load_point_cloud(QString path)
{
    qDebug() << "PointCloudIO::load_point_cloud";
    ifstream fin;
    fin.open (path.toStdString(), std::ifstream::in);
    int NV; // number of vertices
    fin >> NV;
    double x, y, z;
    std::vector<vec3> points;
    points.reserve(NV);
    for(int i=0; i<NV; ++i)
    {
        fin >> x >> y >> z;
        points.push_back(vec3(x,y,z));
    }
    qDebug() << "load_point_cloud finished";
    qDebug() << "NV: " << NV;
    qDebug() << "Number of Points Loaded: " << points.size();
    return points;
}


std::vector<std::vector<vec3>> PointCloudIO::load_point_cloud_ply(QString path)
{
    qDebug() << "PointCloudIO::load_point_cloud";
    ifstream fin;
    string line;
    fin.open(path.toStdString(), std::ifstream::in);
    for (int i = 0; i < 12; ++i)
        getline(fin, line, '\n');
    
    std::vector<vec3> points, normals;
    std::vector<std::vector<vec3>> pc;
    float x, y, z, n_x, n_y, n_z;
    int i = 0;
    while (!fin.eof())
    {
        fin >> x >> y >> z >> n_x >> n_y >> n_z;
        points.push_back(vec3(x, z, -y));
        normals.push_back(vec3(n_x, n_z, -n_y));
        ++i;
    }
    pc.push_back(points);
    pc.push_back(normals);
    qDebug() << "Number of points loaded: " << points.size();
    return pc;

}

void PointCloudIO::save_point_cloud(std::vector<vec3> points, QString path)
{
    qDebug() << "Saving Point Cloud to " << path;
    ofstream fout;
    // TODO check if the file exist
    fout.open(path.toStdString());
    fout << points.size() << endl; // number of vertices
    for(auto v : points)
    {
        fout << v.x <<" "<< v.y <<" "<< v.z <<endl;
    }
    fout.close();
    qDebug() << "Point Cloud Saved to " << path;
}
