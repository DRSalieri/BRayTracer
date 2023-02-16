#ifndef __HIT_H__
#define __HIT_H__

class Hittable;
class Material;

#include "geometry.hpp"
#include "ray.hpp"
#include <memory>


class Hit_record{
public:
    Hit_record(){ init(); }
    ~Hit_record(){}
    void init(){
        happened = false;
        mat_ptr = NULL;
        obj = NULL;
    }

    // make sure dot(r.dir, normal) < 0
    inline void set_face_normal(const Ray& r, const Vec3lf& outward_normal) {
        front_face = (dot(r.dir, outward_normal) < 0);
        if(front_face == false)
        {
            int x = 1;
            x = 2;
        }
        normal = front_face ? outward_normal : -outward_normal;
    }
public:
    bool happened = false;
    bool uv = false;
    Point3lf p;                 // Point
    Vec3lf normal;              // normal
    float t,u,v;                   // t for ray
    bool front_face;
    Material* mat_ptr;
    Hittable* obj;
};

class HittableSample{
public:
    HittableSample(){
        normal = Vec3lf(0,0,0);
        pos = Vec3lf(0,0,0);
        emit = Vec3lf(0,0,0);
        isEmit = false;
    }
public:
    Vec3lf normal;
    Point3lf pos;
    Vec3lf emit;
    bool isEmit;
};

#endif