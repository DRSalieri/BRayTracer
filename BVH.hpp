#ifndef __BVH_H__
#define __BVH_H__

#include "Bound3.hpp"
#include "hittable.hpp"
#include "geometry.hpp"
#include "hit.hpp"

struct BVHNode;
class Hittable;
class Hit_record;

class BVHTree
{
public:
    enum class SplitMethod
    {
        NAIVE,
        SAH
    };

    BVHTree(std::vector<Hittable *> p, int maxPrimsInNode = 1, SplitMethod SplitMethod = SplitMethod::NAIVE);
    ~BVHTree();

    Bound3 WorldBound() const;
    bool hit(const Ray& ray, Hit_record& rec) const;
    bool getHitRec(BVHNode* node, const Ray& ray, Hit_record& rec)const;
    bool IntersectP(const Ray& ray) const;

    BVHNode* recursiveBuild(std::vector<Hittable *> objects);

public:
    BVHNode *root;
    const int maxPrimsInNode;
    const SplitMethod splitMethod;
    std::vector<Hittable *> primitives;
};

struct BVHNode
{
    Bound3 bounds;
    BVHNode *left;
    BVHNode *right;
    Hittable *object;

public:
    int splitAxis = 0;
    int firstPrimOffset = 0;
    int nPrimitives = 0;
    BVHNode()
    {
        bounds = Bound3();
        left = nullptr;
        right = nullptr;
        object = nullptr;
    }
};

#endif