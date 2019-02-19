#include "Occupancy.h"
#include "Scene.h"

Occupancy::Occupancy(vector<vec3> points)
{

    m_pc = points;
    m_bbox = bbox(m_pc);
    m_numx = 16;
    m_numy = 16;
    m_numz = 16;
    create_occupancy();
}

Occupancy::~Occupancy()
{

}

vector<vector<vector<int>>> Occupancy::return_occupancy()
{
    return m_occupancy_grid;
}

vector<vec3> Occupancy::bbox(vector<vec3> points)
{
    Proc_Model p;
    return p.bounding_box(points);
}

void Occupancy::create_occupancy()
{
    m_occupancy_grid.resize(m_numy);
    
    for (int i = 0; i < m_numy; ++i)
    {
        m_occupancy_grid[i].resize(m_numx);
        for (int j = 0; j < m_numx; ++j)
        {
            m_occupancy_grid[i][j].resize(m_numz);
            for (int k = 0; k < m_numz; ++k)
                m_occupancy_grid[i][j][k] = 0;
        }
    }

    float min_x = m_bbox[0].x, min_y = m_bbox[0].y, min_z = m_bbox[0].z, max_x = m_bbox[1].x, max_y = m_bbox[1].y, max_z = m_bbox[1].z;

    for (int i = 0; i < m_pc.size(); ++i)
    {
        vec3 point = m_pc[i];
        int idx_x, idx_y, idx_z;
        idx_x = ((point.x - min_x)*m_numx / (max_x - min_x));
        idx_y = ((point.y - min_y)*m_numy / (max_y - min_y));
        idx_z = ((point.z - min_z)*m_numz / (max_z - min_z));
        if (idx_x == m_numx) idx_x--;
        if (idx_y == m_numy) idx_y--;
        if (idx_z == m_numz) idx_z--;

        if (m_occupancy_grid[idx_y][idx_x][idx_z] == 0)
            m_occupancy_grid[idx_y][idx_x][idx_z] = 1;
    }

    ofstream fout;
    QString folder = "Data/Objs/SpecChairs/simple/test/";
    fout.open((folder + "OccupancyGrid.txt").toStdString());

    for (int i = 0; i < m_numy; ++i)
    {
        vector<vector<int>> plane = m_occupancy_grid[i];
        for (int j = 0; j < m_numx; ++j)
        {
            vector<int> line = plane[j];
            fout << line[0];
            for (int k = 1; k < line.size(); ++k)
                fout << " " << line[k];
            fout << endl;
        }
   }
}
