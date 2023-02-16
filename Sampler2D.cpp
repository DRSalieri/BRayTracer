#include "geometry.hpp"
#include "Sampler2D.hpp"
#include <cmath>

Sampler2D::Sampler2D()
{
}
Sampler2D::~Sampler2D()
{
}

/////////////////////////////////////////////////////////////////////////////////
RandomSampler2D::RandomSampler2D()
{
    width = 1;
    height = 1;
}
RandomSampler2D::RandomSampler2D(float _w, float _h)
{
    width = _w;
    height = _h;
}
RandomSampler2D::~RandomSampler2D()
{
}
bool RandomSampler2D::sample2D(Vec2f &output)
{
    output.x = random_float(0, width);
    output.y = random_float(0, height);
    return true;
}
/////////////////////////////////////////////////////////////////////////////////
UniformSampler2D::UniformSampler2D()
{
    width = height = 1;
    row = column = 10;
}
UniformSampler2D::UniformSampler2D(float _w, float _h, int _row, int _column)
{
    width = _w;
    height = _h;
    row = _row;
    column = _column;
}
UniformSampler2D::~UniformSampler2D()
{
}
bool UniformSampler2D::sample2D(Vec2f &output)
{
    int x = random_int(0, column);
    int y = random_int(0, row);
    output.x = ((float)x + 0.5) * width / column;
    output.y = ((float)y + 0.5) * height / row;
    return true;
}
/////////////////////////////////////////////////////////////////////////////////
BridsonSampler2D::BridsonSampler2D()
{
    width = height = 1;
    r = 0.05;
    r2 = r * r;
    k = 30;
    cell = 0.5 / sqrt(2);
    column = (int)(width / cell);
    row = (int)(height / cell);

    BridsonGenerate();
}
BridsonSampler2D::BridsonSampler2D(float _w, float _h, float _r)
{
    width = _w;
    height = _h;
    r = _r;
    r2 = _r * _r;
    k = 30;
    cell = _r / sqrt(2);
    column = (int)(_w / cell);
    row = (int)(_h / cell);

    BridsonGenerate();
}
bool BridsonSampler2D::sample2D(Vec2f &output)
{
    output = results[random_int(0, results.size())];
    return true;
}
void BridsonSampler2D::BridsonGenerate()
{
    clock_t endTime, startTime;

    startTime = clock();

    std::cerr << "BridsonSampler: Start generate bridson sampler." << std::endl;

    BridsonPoint *infoArray = new BridsonPoint[row * column];
    std::vector<Vec2f> activeArray;

    // 随机生成一个点作为初始点
    Vec2f p = Vec2f(random_float(0, width), random_float(0, height));
    AddPoint(p, activeArray, infoArray);

    while (!activeArray.empty())
    {
        Vec2f activePoint = activeArray[0];
        
        bool k_trial_success = false;
        
        for(int i = 0; i < k ;i++)
        {
            Vec2f genP = random_point_around(activePoint);
            if(Judge(genP, activeArray, infoArray)){
                k_trial_success = true;
                AddPoint(genP, activeArray, infoArray);
            }
        }

        if(k_trial_success == false)
        {
            activeArray.erase(activeArray.begin());
        }
    }

    results.clear();

    for(int i = 0; i < row ; i++)
    {
        for(int j = 0;j < column;j++)
        {
            int ind = index(Vec2i(i, j));

            if(infoArray[ind].mask == true)
            {
                results.push_back(infoArray[ind].p);
            }
        }
    }

    delete[] infoArray;

    endTime = clock();

    std::cerr << "BridsonSampler: Complete generate [" << results.size() << "] points in [" <<  (double)(endTime - startTime) / CLOCKS_PER_SEC << "] seconds." << std::endl;

}
int BridsonSampler2D::index(Vec2i v)
{
    return v.x * column + v.y;
}
int BridsonSampler2D::index(Vec2f v)
{
    return index(coor(v));
}
Vec2i BridsonSampler2D::coor(int i)
{
    return Vec2i(i / column, i % column);
}
Vec2i BridsonSampler2D::coor(Vec2f v)
{
    return Vec2i((int)(v.x / cell), (int)(v.y / cell));
}
void BridsonSampler2D::AddPoint(Vec2f p, std::vector<Vec2f> &activeArray, BridsonPoint *infoArray)
{
    activeArray.push_back(p);
    infoArray[index(p)].mask = true;
    infoArray[index(p)].p = p;
}
bool BridsonSampler2D::Judge(Vec2f p, std::vector<Vec2f> &activeArray, BridsonPoint *infoArray)
{
    // 判断p是否在范围内
    if (p.x < 0 || p.x >= width)
        return false;
    if (p.y < 0 || p.y >= height)
        return false;
    Vec2i co = coor(p);
    if(co.x < 0 || co.x >= column)
        return false;
    if(co.y < 0 || co.y >= row)
        return false;
    if(infoArray[index(co)].mask == true)
        return false;
    // 判断周围的5x5范围内是否满足情况
    int minx = std::max(0, co.x - 2);
    int maxx = std::min(column, co.x + 3);
    int miny = std::max(0, co.y - 2);
    int maxy = std::min(row, co.y + 3);

    for (int i = minx; i < maxx; i++)
    {
        for (int j = miny; j < maxy; j++)
        {
            if(i == co.x && j == co.y) continue;
            // 四个角也可以不用判断
            int ind = index(Vec2i(i,j));
            if(infoArray[ind].mask == false) continue;
            if((infoArray[ind].p - p).norm2() <= r2) return false;
        }
    }

    return true;
}

Vec2f BridsonSampler2D::random_point_around(Vec2f p)
{
    float x1 = random_float(r, r * 2);
    float x2 = random_float(0, 2 * M_PI);
    return Vec2f(p.x + x1 * sin(x2), p.y + x1 * cos(x2));
}