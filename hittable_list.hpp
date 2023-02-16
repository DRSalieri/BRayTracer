#ifndef __HITTABLE_LIST_H__
#define __HITTABLE_LIST_H__

#include "BVH.hpp"

#include <memory>
#include <vector>
#include <limits>

using std::shared_ptr;
using std::make_shared;

class Hittable;

class Hittable_list {
    public:
        Hittable_list() {}
        Hittable_list(Hittable* object) { add(object); }

        void clear() { objects.clear(); }
        void add(Hittable* object) { objects.push_back(object); }

        bool hit(const Ray& r, Hit_record& rec, uint32_t& index);
        void BuildBVH();
    public:
        BVHTree* bvh;
        std::vector<Hittable*> objects;
};

#endif