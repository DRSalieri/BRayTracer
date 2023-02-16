#ifndef __HITTABLE_H__
#define __HITTABLE_H__

class BVHTree;

#include "Sampler2D.hpp"
#include "hit.hpp"
#include "ray.hpp"
#include "material.hpp"
#include "Bound3.hpp"

bool rayTriangleIntersect(const Vec3lf &v0, const Vec3lf &v1, const Vec3lf &v2, const Vec3lf &orig, const Vec3lf &dir, float &tnear, float &u, float &v);

#pragma region Hittable

class Hittable
{
public:
    virtual bool get_uv(const Hit_record& rec, Vec2lf& res);
    virtual bool sample(Hittable_list &scene, const Hit_record &rec, color &output);
    virtual bool hit(const Ray &ray, Hit_record &rec) = 0;
    // virtual bool hit(const Ray &ray, Hit_record &rec, float &, uint32_t &) = 0;
    virtual Bound3 getBounds() = 0;
    virtual void outputString() = 0;
    virtual Material *getMaterial();

    //path tracing
    Material* mat_ptr;
    float area;
    virtual float getArea();
    virtual bool PTSample(HittableSample& sam, float& pdf);
};

#pragma endregion

#pragma region Sphere

class Sphere : public Hittable
{
public:
    Sphere(){};
    Sphere(Point3lf cen, float r, Material *m)
    {
        center = cen;
        radius = r;
        mat_ptr = m;
        radius2 = radius * radius;
        area = 4 * M_PI * r * r;
    };

public:
    virtual bool hit(const Ray &ray, Hit_record &rec) override;
    // virtual bool hit(const Ray &ray, Hit_record &rec, float &, uint32_t &) override;
    virtual void outputString() override;

    virtual Bound3 getBounds() override;

public:
    Point3lf center;
    float radius, radius2;
};
#pragma endregion

#pragma region Triangle
class Triangle : public Hittable
{
public:
    Triangle(Vec3lf _v0, Vec3lf _v1, Vec3lf _v2, Material *_m, Sampler2D* _sampler = NULL);
    ~Triangle() {}

public:
    void setNormal(Vec3lf _n0, Vec3lf _n1, Vec3lf _n2);
    void setTex(Vec2lf _t0, Vec2lf _t1, Vec2lf _t2);
    virtual bool get_uv(const Hit_record& rec, Vec2lf& res) override;
    virtual bool hit(const Ray &ray, Hit_record &rec) override;
    // virtual bool hit(const Ray &ray, Hit_record &rec, float &, uint32_t &) override;
    bool hit(const Ray &ray, Hit_record &rec, float &t_min);
    virtual void outputString() override;

    virtual bool PTSample(HittableSample& sam, float& pdf) override;

    virtual Bound3 getBounds() override;

public:
    Vec3lf vertice[3];
    Vec3lf normal[3];
    Vec2lf tex[3];
    Vec3lf e1, e2;
    Vec3lf face_normal;
    float area;
    Sampler2D* sampler;
};

#pragma endregion

#pragma region Polygon

class Polygon: public Hittable
{
public:
    Polygon(const std::vector<Point3lf>& p_list, Material* _m, Sampler2D* _sampler = NULL);
    ~Polygon(){};

public:
    virtual bool hit(const Ray& ray, Hit_record& rec) override;
    virtual void outputString() override;
    virtual bool PTSample(HittableSample& sam, float& pdf) override;
    virtual Bound3 getBounds() override;
    virtual float getArea() override;
public:
    Bound3 bounding_box;
    std::vector<Triangle> triangles;
    Sampler2D* sampler;
    float area;
};

#pragma endregion

#pragma region MeshTriangle

class MeshTriangle : public Hittable
{
public:
    /*
    MeshTriangle(
        const Vec3lf *verts, const uint32_t *vertsIndex,
        const uint32_t &numTris, const Vec2lf *st);
    */
    MeshTriangle(
    const std::string &filename,
    const std::string &mtlpath,
    Material *m,
    Matrix modelTransform,
    bool tex);

public:
    virtual bool hit(const Ray &ray, Hit_record &rec) override;
    // virtual bool hit(const Ray &ray, Hit_record &rec, float &, uint32_t &) override;
    virtual void outputString() override;

    virtual Bound3 getBounds() override;

public:
    Bound3 bounding_box;
    BVHTree *bvh;
    uint32_t numTriangles;
    std::vector<Triangle> triangles;

};

#pragma endregion

#pragma region Cube

class Cube : public Hittable
{
public:
    Cube();
    Cube(Vec3lf p1, Vec3lf p2, Material *m);

public:
    virtual bool hit(const Ray &ray, Hit_record &rec) override;
    // virtual bool hit(const Ray &ray, Hit_record &rec, float &, uint32_t &) override;
    virtual void outputString() override;

    virtual Bound3 getBounds() override;
    virtual bool PTSample(HittableSample &sam, float &pdf);
public:
    Bound3 bounding_box;
    Vec3lf pMin;
    Vec3lf pMax;
};

#pragma endregion

#pragma region Rect

class Rect : public Hittable
{
public:
    Rect(Vec3lf _v0, Vec3lf _v1, Vec3lf _v2, Vec3lf _v3, Material *_m, float _step = 0.02)
    {
        v[0] = _v0;
        v[1] = _v1;
        v[2] = _v2;
        v[3] = _v3;
        step = _step;
        normal = cross(_v1 - _v0, _v3 - _v0);
        mat_ptr = _m;
    }

public:
    virtual bool hit(const Ray &ray, Hit_record &rec) override;
    virtual bool sample(Hittable_list &scene, const Hit_record &rec, color &output) override;
    virtual void outputString() override;
    virtual Bound3 getBounds() override;
    virtual Material *getMaterial() override;

public:
    Vec3lf v[4]; // 0 1 2 3 v01 v03 sample
    Vec3lf normal;
    float step;
    Material *mat_ptr;
};

#pragma endregion

#pragma region Disk

class Disk : public Hittable
{
public:
    Disk(Point3lf _p, Vec3lf _x, Vec3lf _n, float _r, Material *_m, float _stepr, float _steptheta, Sampler2D* _sampler = NULL)
    {
        p = _p;
        normal = _n.normalize();
        r = _r;
        mat_ptr = _m;
        x = _x.normalize() * _r;
        y = cross(_x, _n).normalize() * _r;
        r2 = _r * _r;
        area = 2 * M_PI * _r * _r;
        sampler = _sampler;

        step_r = _stepr;
        step_theta = _steptheta;
    }

public:
    virtual bool hit(const Ray &ray, Hit_record &rec) override;
    virtual bool sample(Hittable_list &scene, const Hit_record &rec, color &output) override;
    virtual void outputString() override;
    virtual Bound3 getBounds() override;

    virtual bool PTSample(HittableSample& sam, float& pdf) override;

public:
    Point3lf p;
    Vec3lf x, y;
    Vec3lf normal;
    float r, r2;
    float step_r, step_theta;
    Sampler2D* sampler;
};

#pragma endregion

#pragma region Cylinder

class Cylinder : public Hittable
{
public:
    Cylinder(Point3lf _p, float _r, float _h, Material *_m)
    {
        hmin = _p.y;
        hmax = _p.y + _h;
        h = _h;
        r = _r;
        mat_ptr = _m;
        p = _p;
        r2 = _r * _r;
    }

public:
    virtual bool hit(const Ray &r_in, Hit_record &rec) override;
    virtual void outputString() override;
    virtual Bound3 getBounds() override;

public:
    float r, r2, h;
    Point3lf p;
    float hmin, hmax;
};

#pragma endregion

#endif