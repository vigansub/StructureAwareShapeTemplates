#ifndef POINTCLOUDIO_H
#define POINTCLOUDIO_H
#include <vector>
#include "Math.h"
#include <QString>

namespace PointCloudIO
{
std::vector<vec3> load_point_cloud(QString path);
std::vector<std::vector<vec3>> load_point_cloud_ply(QString path);
void save_point_cloud(std::vector<vec3> points, QString path);
}

#endif // POINTCLOUDIO_H
