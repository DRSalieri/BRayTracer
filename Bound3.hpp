#ifndef __BOUNDS3_H__
#define __BOUNDS3_H__

#include "ray.hpp"
#include "geometry.hpp"
#include <array>

class Bound3
{
public:
    Bound3();
    Bound3(const Vec3lf p) : pMin(p), pMax(p) {}
    Bound3(const Vec3lf p1, const Vec3lf p2);
    Vec3lf Diagonal() const;
    int maxExtent() const;
    float surfaceArea() const;
    Vec3lf Centroid();                                  // rertun the centroid coordinates
    Bound3 Intersect(const Bound3 &b);
    Vec3lf Offset(const Vec3lf &p);
    bool Overlaps(const Bound3 &b1, const Bound3 &b2);
    bool Inside(const Vec3lf &p, const Bound3 &b);

public:
    inline const Vec3lf &operator[](int i) const;
    inline bool IntersectP(const Ray &ray, const Vec3lf &invDir,
                           const std::array<int, 3> &dirisNeg) const;
    inline Bound3 Union(const Bound3 &b);
    inline Bound3 Union(const Vec3lf &v);

public:
    Vec3lf pMin, pMax;
};

inline const Vec3lf &Bound3::operator[](int i) const
{
    return (i == 0) ? pMin : pMax;
}

inline bool Bound3::IntersectP(const Ray &ray, const Vec3lf &invDir,
                               const std::array<int, 3> &dirIsNeg) const
{
    float t_min_x = (pMin.x - ray.orig.x) * invDir[0];
    float t_min_y = (pMin.y - ray.orig.y) * invDir[1];
    float t_min_z = (pMin.z - ray.orig.z) * invDir[2];
    float t_max_x = (pMax.x - ray.orig.x) * invDir[0];
    float t_max_y = (pMax.y - ray.orig.y) * invDir[1];
    float t_max_z = (pMax.z - ray.orig.z) * invDir[2];

    // 如果方向为负方向，交换max与min
    if (!dirIsNeg[0])
    {
        float t = t_min_x;
        t_min_x = t_max_x;
        t_max_x = t;
    }
    if (!dirIsNeg[1])
    {
        float t = t_min_y;
        t_min_y = t_max_y;
        t_max_y = t;
    }
    if (!dirIsNeg[2])
    {
        float t = t_min_z;
        t_min_z = t_max_z;
        t_max_z = t;
    }

    float t_enter = std::max(t_min_x, std::max(t_min_y, t_min_z));
    float t_exit = std::min(t_max_x, std::min(t_max_y, t_max_z));
    if (t_enter <= t_exit && t_exit >= 0)
    {
        return true;
    }
    else
        return false;
}

inline Bound3 Bound3::Union(const Bound3 &b)
{
    this->pMin = Vec3lf::Min(this->pMin, b.pMin);
    this->pMax = Vec3lf::Max(this->pMax, b.pMax);
    return *this;
}
inline Bound3 Bound3::Union(const Vec3lf &p)
{
    this->pMin = Vec3lf::Min(this->pMin, p);
    this->pMax = Vec3lf::Max(this->pMax, p);
    return *this;
}
Bound3 Union(const Bound3& b1, const Bound3& b2);

#endif