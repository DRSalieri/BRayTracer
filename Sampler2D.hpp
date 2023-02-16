#ifndef __SAMPLER2D_H__
#define __SAMPLER2D_H__

#include "geometry.hpp"
#include <vector>
#include <deque>

class Sampler2D
{
public:
    Sampler2D();
    ~Sampler2D();
    virtual bool sample2D(Vec2f& output) = 0;
public:
    float width;
    float height;
};

class RandomSampler2D : public Sampler2D
{
public:
    RandomSampler2D();
    RandomSampler2D(float _w, float _h);
    ~RandomSampler2D();
    virtual bool sample2D(Vec2f& output);
};

class UniformSampler2D : public Sampler2D
{
public:
    UniformSampler2D();
    UniformSampler2D(float _w, float _h, int _row,int _column);
    ~UniformSampler2D();
    virtual bool sample2D(Vec2f& output);
public:
    int row;
    int column;
};

class BridsonSampler2D : public Sampler2D
{
public:
    struct BridsonPoint{
        bool mask = false;
        Vec2f p;
    };
public:
    BridsonSampler2D();
    BridsonSampler2D(float _w, float _h, float _r);
    ~BridsonSampler2D();
    virtual bool sample2D(Vec2f& output);    
private:
    void BridsonGenerate();
    Vec2f random_point_around(Vec2f p);
    int index(Vec2i v);
    int index(Vec2f v);
    Vec2i coor(int i);
    Vec2i coor(Vec2f v);
    bool Judge(Vec2f p, std::vector<Vec2f>& activeArray, BridsonPoint* infoArray);      // 判断一个点是否合理
    void AddPoint(Vec2f p, std::vector<Vec2f>& activeArray, BridsonPoint* infoArray);
public:
    float r, r2;            
    float cell;         // 网格边长
    int k;
    int row;
    int column;
    std::vector<Vec2f> results;
};

#endif