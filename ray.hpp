#ifndef __RAY_H__
#define __RAY_H__

#include "geometry.hpp"

class Ray{
public:
    Ray() {}
    Ray(const Point3lf& origin,const Vec3lf& direction)
        : orig(origin), dir(direction){}

    Point3lf origin() const { return orig; }
    Vec3lf direction() const { return dir; } 

    Point3lf at(float t) const{
        return orig + t * dir;
    }

    Ray& epsilon(const Vec3lf& n) {
        EpsilonRay(orig, dir, n);
        return *this;
    }
    
public:
    bool canReflect = true;
    Point3lf orig;
    Vec3lf dir;
};

#endif