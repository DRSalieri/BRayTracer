#include "BVH.hpp"
#include <algorithm>

BVHTree::BVHTree(std::vector<Hittable *> p, int maxPrimsInNode,
                 SplitMethod SplitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    time_t start, stop;
    time(&start);
    if (primitives.empty())
        return;

    root = recursiveBuild(primitives);

    time(&stop);
    float diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    /*
    printf(
        "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
        hrs, mins, secs);
        */
}

BVHNode *BVHTree::recursiveBuild(std::vector<Hittable *> objects)
{
    // root
    BVHNode *node = new BVHNode();

    // scene Bounds
    Bound3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds.Union(objects[i]->getBounds());

    // 3.when to terminate
    if (objects.size() == 1)
    {
        // leaf _BVHNode_
        node->bounds = objects[0]->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        return node;
    }

    else if (objects.size() == 2)
    {
        // just 2 object, 1 left, 1 right
        node->left = recursiveBuild(std::vector<Hittable *>{objects[0]});
        node->right = recursiveBuild(std::vector<Hittable *>{objects[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }

    else
    {
        // 1.how to split: max Extent of x,y,z
        // compute the maxExtent of centroidBounds
        Bound3 centroidBounds;
        for (int i = 0; i < objects.size(); ++i)
            centroidBounds.Union(objects[i]->getBounds().Centroid());
        int dim = centroidBounds.maxExtent();
        switch (dim)
        {
        case 0:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2)
                      { return f1->getBounds().Centroid().x <
                               f2->getBounds().Centroid().x; });
            break;
        case 1:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2)
                      { return f1->getBounds().Centroid().y <
                               f2->getBounds().Centroid().y; });
            break;
        case 2:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2)
                      { return f1->getBounds().Centroid().z <
                               f2->getBounds().Centroid().z; });
            break;
        }

        // 2.where to cut: the middle number of objects

        auto beginning = objects.begin();
        auto middling = objects.begin() + (objects.size() / 2);
        auto ending = objects.end();

        auto leftshapes = std::vector<Hittable *>(beginning, middling);
        auto rightshapes = std::vector<Hittable *>(middling, ending);

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild(leftshapes);
        node->right = recursiveBuild(rightshapes);

        node->bounds = Union(node->left->bounds, node->right->bounds);
    }

    return node;
}

bool BVHTree::hit(const Ray &ray, Hit_record &rec) const
{
    if (root == NULL){
        return false;
    }
    return BVHTree::getHitRec(root, ray, rec);
}

bool BVHTree::getHitRec(BVHNode *node, const Ray &ray, Hit_record &rec) const
{
    if (node == NULL)
        return false;
    Vec3lf invDir(1. / ray.dir.x, 1. / ray.dir.y, 1. / ray.dir.z);
    std::array<int, 3> dirIsNeg;
    dirIsNeg[0] = ray.dir.x > 0;
    dirIsNeg[1] = ray.dir.y > 0;
    dirIsNeg[2] = ray.dir.z > 0;

    // if not hit parent bounds, return;
    if (!node->bounds.IntersectP(ray, invDir, dirIsNeg))
    {
        return false;
    }
    // leaf
    if (node->left == nullptr && node->right == nullptr)
    {
        return node->object->hit(ray, rec);
    }
    else
    {
        // recursive
        Hit_record recL;
        bool interL = getHitRec(node->left, ray, recL);
        Hit_record recR;
        bool interR = getHitRec(node->right, ray, recR);

        if (interL == true)
        {
            if (interR == true)
            {
                rec = recL.t < recR.t ? recL : recR;
                return true;
            }
            else
            {
                rec = recL;
                return true;
            }
        }
        else
        {
            if (interR == true)
            {
                rec = recR;
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    return false;
}
