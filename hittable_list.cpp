#include "hittable_list.hpp"

bool Hittable_list::hit(const Ray &r, Hit_record &rec, uint32_t &index)
{

    Hit_record temp_rec;
    bool hit_anything = false;
    float t_min = std::numeric_limits<float>().infinity();

    return this->bvh->hit(r, rec);
    /*
        for (const auto& object : objects) {
            if (object->hit(r, temp_rec) && temp_rec.t < t_min) {
                // object->outputString();
                hit_anything |= true;
                t_min = temp_rec.t;
                rec = temp_rec;
            }
        }
        */

    // return hit_anything;
}

void Hittable_list::BuildBVH()
{
    bvh = new BVHTree(objects);
}
