#include "ray.hpp"
#include "geometry.hpp"
#include "Bound3.hpp"
#include <limits>

Bound3::Bound3()
{
    float minNum = std::numeric_limits<float>::lowest();
    float maxNum = std::numeric_limits<float>::max();
    pMax = Vec3lf(minNum, minNum, minNum);
    pMin = Vec3lf(maxNum, maxNum, maxNum);
}

Bound3::Bound3(const Vec3lf p1, const Vec3lf p2)
{
    pMin = Vec3lf(fmin(p1.x, p2.x), fmin(p1.y, p2.y), fmin(p1.z, p2.z));
    pMax = Vec3lf(fmax(p1.x, p2.x), fmax(p1.y, p2.y), fmax(p1.z, p2.z));
}

Vec3lf Bound3::Diagonal() const
{
    return pMax - pMin;
}
/***
 * maxtent: the max extent of x,y,z
 * 0 - x
 * 1 - y
 * 2 - z
 **/
int Bound3::maxExtent() const
{
    Vec3lf d = Diagonal();
    if (d.x > d.y && d.x > d.z)
        return 0;
    else if (d.y > d.z)
        return 1;
    else
        return 2;
}
float Bound3::surfaceArea() const
{
    Vec3lf d = Diagonal();
    return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
}
Vec3lf Bound3::Centroid() { return 0.5 * pMin + 0.5 * pMax; }
Bound3 Bound3::Intersect(const Bound3 &b)
{
    return Bound3(Vec3lf(fmax(pMin.x, b.pMin.x), fmax(pMin.y, b.pMin.y),
                         fmax(pMin.z, b.pMin.z)),
                  Vec3lf(fmin(pMax.x, b.pMax.x), fmin(pMax.y, b.pMax.y),
                         fmin(pMax.z, b.pMax.z)));
}

Vec3lf Bound3::Offset(const Vec3lf &p)
{
    Vec3lf o = p - pMin;
    if (pMax.x > pMin.x)
        o.x /= pMax.x - pMin.x;
    if (pMax.y > pMin.y)
        o.y /= pMax.y - pMin.y;
    if (pMax.z > pMin.z)
        o.z /= pMax.z - pMin.z;
    return o;
}
bool Bound3::Overlaps(const Bound3 &b1, const Bound3 &b2)
{
    bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
    bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
    bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
    return (x && y && z);
}
bool Bound3::Inside(const Vec3lf &p, const Bound3 &b)
{
    return (p.x >= b.pMin.x && p.x <= b.pMax.x && p.y >= b.pMin.y &&
            p.y <= b.pMax.y && p.z >= b.pMin.z && p.z <= b.pMax.z);
}
Bound3 Union(const Bound3 &b1, const Bound3 &b2)
{
    Bound3 ret;
    ret.pMin = Vec3lf::Min(b1.pMin, b2.pMin);
    ret.pMax = Vec3lf::Max(b1.pMax, b2.pMax);
    return ret;
}