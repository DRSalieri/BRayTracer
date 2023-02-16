#define TINYOBJLOADER_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "ray.hpp"
#include "hit.hpp"
#include "hittable.hpp"
#include "material.hpp"
#include "geometry.hpp"
#include "Transform.hpp"
#include "tiny_obj_loader.h"
#include "BVH.hpp"
#include <cstring>
#include <cfloat>
#include <limits>
#include "hittable_list.hpp"
#include "stb_image.h"

#pragma region helper

void VecTransform(const Matrix &mat, Vec3lf &vert)
{
    vec<4, float> temp_vert = embed<4, 3, float>(vert, 1);
    temp_vert = mat * temp_vert;
    temp_vert = temp_vert / temp_vert[3];
    vert = proj<3, 4, float>(temp_vert);
}

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
    float discr = b * b - 4 * a * c;
    if (discr < 0)
        return false;
    else if (discr == 0)
        x0 = x1 = -0.5 * b / a;
    else
    {
        float q = (b > 0) ? -0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr));
        x0 = q / a;
        x1 = c / q;
    }
    if (x0 > x1)
        std::swap(x0, x1);
    return true;
}

bool rayTriangleIntersect(
    const Vec3lf &v0, const Vec3lf &v1, const Vec3lf &v2, const Vec3lf &orig, const Vec3lf &dir, float &tnear, float &u, float &v)
{

    Vec3lf e1 = v1 - v0;
    Vec3lf e2 = v2 - v0;
    Vec3lf s = orig - v0;
    Vec3lf s1 = cross(dir, e2);
    Vec3lf s2 = cross(s, e1);
    float s1e1 = dot(s1, e1);

    tnear = dot(s2, e2) / s1e1;
    u = dot(s1, s) / s1e1;
    v = dot(s2, dir) / s1e1;

    if (u >= 0 && v >= 0 && (1 - u - v) >= 0 && tnear >= 0)
        return true;
    else
        return false;
}

bool rayFlatIntersect3D(
    const Vec3lf &orig, const Vec3lf &dir, const Point3lf &p, const Vec3lf &n, float &tnear, bool &inFlat)
{
    inFlat = false;
    float div = dot(n, dir);
    if (div == 0)
    {
        if (dot(orig - p, n) == 0)
        {
            tnear = 0;
            inFlat = true;
            return true;
        }
        else
            return false;
    }
    else
        div = 1.0 / div;

    tnear = (dot(p, n) - dot(orig, n)) * div;

    if (tnear >= 0)
        return true;
    else
        return false;
}

bool rayLineDist(
    const Vec3lf &orig, const Vec3lf &dir, const Point3lf &p, float &t, float &dist2)
{
    dist2 = cross(p - orig, dir.normalize()).norm2();
    float op = (orig - p).norm2();
    t = (op - dist2) / dir.norm();
    if (dot(dir, p - orig) < 0)
        t = -t;
    return true;
}

#pragma endregion
/////////////////////////////////////////////////////////////////////////////////
#pragma region Hittable

bool Hittable::sample(Hittable_list &scene, const Hit_record &rec, color &output)
{
    return false;
}

Material *Hittable::getMaterial()
{
    return mat_ptr;
}
bool Hittable::get_uv(const Hit_record &rec, Vec2lf &res)
{
    return false;
}
float Hittable::getArea()
{
    return area;
}
bool Hittable::PTSample(HittableSample &sam, float &pdf)
{
    return false;
}

#pragma endregion
/////////////////////////////////////////////////////////////////////////////////
#pragma region Sphere
bool Sphere::hit(const Ray &ray, Hit_record &rec)
{
    Vec3lf L = ray.orig - center;

    float a = dot(ray.dir, ray.dir);
    float b = 2 * dot(ray.dir, L);
    float c = dot(L, L) - radius2;
    float t0, t1;

    if (!solveQuadratic(a, b, c, t0, t1))
        return false;
    if (t0 == t1)
        return false;
    if (t0 < EPSILON)
        t0 = t1;
    if (t0 < EPSILON)
        return false;

    // rec
    rec.happened = true;
    rec.p = Vec3lf(ray.orig + ray.dir * t0);
    rec.t = t0;
    rec.obj = this;
    rec.mat_ptr = mat_ptr;
    rec.set_face_normal(ray, (rec.p - center) / radius);

    return true;
}

Bound3 Sphere::getBounds()
{
    return Bound3(Vec3lf(center.x - radius, center.y - radius, center.z - radius),
                  Vec3lf(center.x + radius, center.y + radius, center.z + radius));
}

void Sphere::outputString()
{
    std::cerr << "Sphere" << std::endl;
}

#pragma endregion
/////////////////////////////////////////////////////////////////////////////////
#pragma region Triangle
Triangle::Triangle(Vec3lf _v0, Vec3lf _v1, Vec3lf _v2, Material *_m, Sampler2D *_sampler)
{
    mat_ptr = _m;
    vertice[0] = _v0;
    vertice[1] = _v1;
    vertice[2] = _v2;
    e1 = _v1 - _v0;
    e2 = _v2 - _v0;
    face_normal = cross(e1, e2).normalize();
    area = cross(e1, e2).norm() * 0.5f;
    sampler = _sampler;
}

bool Triangle::hit(const Ray &r, Hit_record &rec)
{
    float t, u, v;
    if (rayTriangleIntersect(vertice[0], vertice[1], vertice[2], r.orig, r.dir, t, u, v))
    {
        rec.happened = true;
        rec.uv = true;
        rec.t = t;
        rec.u = u;
        rec.v = v;
        rec.normal = dot(face_normal, r.dir) < 0 ? face_normal : -face_normal;
        rec.p = r.at(rec.t);
        rec.mat_ptr = mat_ptr;
        rec.obj = this;
        return true;
    }
    return false;
}

bool Triangle::hit(const Ray &r, Hit_record &rec, float &t_min)
{
    float t, u, v;
    if (rayTriangleIntersect(vertice[0], vertice[1], vertice[2], r.orig, r.dir, t, u, v) && t < t_min)
    {
        t_min = t;

        rec.happened = true;
        rec.uv = true;
        rec.t = t;
        rec.u = u;
        rec.v = v;
        rec.normal = dot(face_normal, r.dir) < 0 ? face_normal : -face_normal;
        rec.p = r.at(rec.t);
        rec.mat_ptr = mat_ptr;
        rec.obj = this;

        return true;
    }
    return false;
}

bool Triangle::get_uv(const Hit_record &rec, Vec2lf &res)
{
    if (rec.uv == false)
        return false;
    res = interpolate(rec.u, rec.v, 1 - rec.u - rec.v, tex[0], tex[1], tex[2], 1);
    return true;
}

Bound3 Triangle::getBounds()
{
    return Bound3(vertice[0], vertice[1]).Union(vertice[2]);
}

void Triangle::outputString()
{
    std::cerr << "Triangle" << std::endl;
}

bool Triangle::PTSample(HittableSample &sam, float &pdf)
{
    Vec2f ret;
    sampler->sample2D(ret);
    // https://stackoverflow.com/questions/4778147/sample-random-point-in-triangle
    float sqr1 = sqrt(ret[0]);
    float sqr2 = sqrt(ret[1]);
    Vec3f pos = vertice[0] * (1.f - sqr1) + vertice[1] * (sqr1 * (1 - ret[1])) + vertice[2] * (sqr1 * ret[1]);

    sam.pos = pos;
    sam.normal = face_normal;
    sam.isEmit = mat_ptr->hasEmission();
    sam.emit = mat_ptr->getEmission();

    pdf = 1.0f / area;
    return true;
}

#pragma endregion
/////////////////////////////////////////////////////////////////////////////////
#pragma region MeshTriangle

MeshTriangle::MeshTriangle(
    const std::string &filename,
    const std::string &mtlpath,
    Material *m,
    Matrix modelTransform,
    bool tex)
{

    // Init
    numTriangles = 0;
    mat_ptr = m;

    // OBJ loader
    tinyobj::ObjReaderConfig reader_config;
    reader_config.mtl_search_path = mtlpath;
    tinyobj::ObjReader reader;

    if (!reader.ParseFromFile(filename, reader_config))
    {
        if (!reader.Error().empty())
        {
            std::cerr << "TinyObjReader: " << reader.Error();
        }
        exit(1);
    }

    if (!reader.Warning().empty())
    {
        std::cerr << "TinyObjReader: " << reader.Warning();
    }

    auto &attrib = reader.GetAttrib();
    auto &shapes = reader.GetShapes();
    auto &materials = reader.GetMaterials();

    Vec3lf min_vert = Vec3lf{std::numeric_limits<float>::infinity(),
                             std::numeric_limits<float>::infinity(),
                             std::numeric_limits<float>::infinity()};
    Vec3lf max_vert = Vec3lf{-std::numeric_limits<float>::infinity(),
                             -std::numeric_limits<float>::infinity(),
                             -std::numeric_limits<float>::infinity()};

    std::array<Vec3lf, 3> shape_vertices;
    std::array<Vec3lf, 3> shape_normals;
    std::array<Vec2lf, 3> shape_texs;

    std::vector<Diffuse *> mats;

    // loop over materials
    for (size_t m = 0; m < materials.size(); m++)
    {
        mats.emplace_back(new Diffuse());
        if (!materials[m].diffuse_texname.empty())
        {
            Texture *texture = new Texture();
            texture->data = stbi_load((mtlpath + materials[m].diffuse_texname).c_str(), &(texture->w), &(texture->h), &(texture->n), 3);
            mats[m]->texture = texture;
        }
    }

    // loop over shapes
    for (size_t s = 0; s < shapes.size(); s++)
    {
        std::cerr << "TinyObjReader[shape]: " << shapes[s].name
                  << "( face: " << shapes[s].mesh.num_face_vertices.size() << ")"
                  << "( vertex: " << shapes[s].mesh.indices.size() << ")"
                  << std::endl;
        numTriangles += shapes[s].mesh.num_face_vertices.size();
        size_t index_offset = 0;
        // loop over faces
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++)
        {
            bool haveNormal = false;
            bool haveTex = false;
            size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);
            // check: must be a triangle face
            if (fv != 3)
            {
                std::cerr << "TinyObjReader: "
                          << "Not a Triangle mesh.\n";
                exit(1);
            }
            // loop over vertices
            for (size_t v = 0; v < fv; v++)
            {
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
                tinyobj::real_t vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
                tinyobj::real_t vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];
                Vec3lf vert = Vec3lf(vx, vy, vz);
                VecTransform(modelTransform, vert);

                min_vert = Vec3lf(std::min(min_vert.x, vert.x),
                                  std::min(min_vert.y, vert.y),
                                  std::min(min_vert.z, vert.z));
                max_vert = Vec3lf(std::max(max_vert.x, vert.x),
                                  std::max(max_vert.y, vert.y),
                                  std::max(max_vert.z, vert.z));

                shape_vertices[v] = vert;

                if (idx.normal_index >= 0)
                {
                    tinyobj::real_t nx = attrib.normals[3 * size_t(idx.normal_index) + 0];
                    tinyobj::real_t ny = attrib.normals[3 * size_t(idx.normal_index) + 1];
                    tinyobj::real_t nz = attrib.normals[3 * size_t(idx.normal_index) + 2];
                    Vec3lf norm = Vec3lf(nx, ny, nz);
                    shape_normals[v] = norm;
                    haveNormal = true;
                }

                // Check if `texcoord_index` is zero or positive. negative = no texcoord data
                if (idx.texcoord_index >= 0)
                {
                    tinyobj::real_t tx = attrib.texcoords[2 * size_t(idx.texcoord_index) + 0];
                    tinyobj::real_t ty = attrib.texcoords[2 * size_t(idx.texcoord_index) + 1];
                    Vec2lf tex = Vec2lf(tx, ty);
                    shape_texs[v] = tex;
                    haveTex = false;
                }
            }
            // std::cerr << shape_vertices[0] << shape_vertices[1] << shape_vertices[2] << std::endl;
            index_offset += fv;
            Triangle tri = Triangle(shape_vertices[0], shape_vertices[1], shape_vertices[2], m);
            if (haveNormal)
                tri.setNormal(shape_normals[0], shape_normals[1], shape_normals[2]);
            if (haveTex)
                tri.setTex(shape_texs[0], shape_texs[1], shape_texs[2]);

            if (tex)
                tri.mat_ptr = mats[shapes[s].mesh.material_ids[f]];
            else
                tri.mat_ptr = m;

            triangles.emplace_back(tri);
        }
    }

    bounding_box = Bound3(min_vert, max_vert);

    // Build BVH Tree
    std::vector<Hittable *> ptrs;
    for (auto &tri : triangles)
        ptrs.push_back(&tri);
    bvh = new BVHTree(ptrs);
}

void Triangle::setNormal(Vec3lf _n0, Vec3lf _n1, Vec3lf _n2)
{
    normal[0] = _n0;
    normal[1] = _n1;
    normal[2] = _n2;
}
void Triangle::setTex(Vec2lf _t0, Vec2lf _t1, Vec2lf _t2)
{
    tex[0] = _t0;
    tex[1] = _t1;
    tex[2] = _t2;
}

bool MeshTriangle::hit(const Ray &ray, Hit_record &rec)
{
    return bvh->hit(ray, rec);

    bool intersect = false;
    float t_min = std::numeric_limits<float>().infinity();

    for (Triangle &tri : triangles)
    {
        if (tri.hit(ray, rec, t_min))
        {
            intersect |= true;
        }
    }
    return intersect;
}

Bound3 MeshTriangle::getBounds()
{
    return bounding_box;
}

void MeshTriangle::outputString()
{
    std::cerr << "MeshTriangle" << std::endl;
}
#pragma endregion
/////////////////////////////////////////////////////////////////////////////////
#pragma region Cube
Cube::Cube()
{
    pMin = Vec3lf(0, 0, 0);
    pMax = Vec3lf(0, 0, 0);
    area = 0;
    mat_ptr = nullptr;
}

Cube::Cube(Vec3lf p1, Vec3lf p2, Material *m)
{
    pMin = Vec3lf::Min(p1, p2);
    pMax = Vec3lf::Max(p1, p2);
    Vec3lf del = pMax - pMin;
    area = del.x * del.y * del.z;
    if (area < 0)
        area = -area;
    mat_ptr = m;
}

bool Cube::hit(const Ray &ray, Hit_record &rec)
{
    Vec3lf invDir(1. / ray.dir.x, 1. / ray.dir.y, 1. / ray.dir.z);
    std::array<int, 3> dirIsNeg;
    dirIsNeg[0] = ray.dir.x > 0;
    dirIsNeg[1] = ray.dir.y > 0;
    dirIsNeg[2] = ray.dir.z > 0;

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

    /**
     * 1 - x
     * 2 - y
     * 3 - z
     */
    int normal_face = 0;
    float t_enter;

    if (t_min_x > t_min_y)
    {
        if (t_min_x > t_min_z)
        {
            normal_face = 1;
            t_enter = t_min_x;
        }
        else
        {
            normal_face = 3;
            t_enter = t_min_z;
        }
    }
    else
    {
        if (t_min_y > t_min_z)
        {
            normal_face = 2;
            t_enter = t_min_y;
        }
        else
        {
            normal_face = 3;
            t_enter = t_min_z;
        }
    }

    float t_exit = std::min(t_max_x, std::min(t_max_y, t_max_z));

    if (t_enter <= t_exit && t_exit >= 0)
    {
        rec.happened = true;

        switch (normal_face)
        {
        case 1:
            rec.normal = dirIsNeg[0] ? Vec3lf(-1, 0, 0) : Vec3lf(1, 0, 0);
            break;
        case 2:
            rec.normal = dirIsNeg[1] ? Vec3lf(0, -1, 0) : Vec3lf(0, 1, 0);
            break;
        case 3:
            rec.normal = dirIsNeg[2] ? Vec3lf(0, 0, -1) : Vec3lf(0, 0, 1);
            break;
        default:
            break;
        }

        rec.t = t_enter > 0 ? t_enter : t_exit;
        rec.p = Vec3lf(ray.orig + ray.dir * rec.t);
        rec.obj = this;
        rec.mat_ptr = mat_ptr;
        return true;
    }
    else
        return false;
}

void Cube::outputString()
{
    std::cerr << "Cube" << std::endl;
}

Bound3 Cube::getBounds()
{
    return Bound3(pMin, pMax);
}

bool Cube::PTSample(HittableSample &sam, float &pdf)
{
    // TODO:Cube Sample
    return false;
}
#pragma endregion
/////////////////////////////////////////////////////////////////////////////////
#pragma region Rect
bool Rect::hit(const Ray &ray, Hit_record &rec)
{
    float t, a, b;
    if (rayTriangleIntersect(v[0], v[1], v[2], ray.orig, ray.dir, t, a, b))
    {
        rec.happened = true;
        rec.t = t;
        rec.normal = dot(normal, ray.dir) < 0 ? normal : -normal;
        rec.p = ray.at(rec.t);
        rec.mat_ptr = mat_ptr;
        rec.obj = this;
        return true;
    }
    else if (rayTriangleIntersect(v[0], v[2], v[3], ray.orig, ray.dir, t, a, b))
    {
        rec.happened = true;
        rec.t = t;
        rec.normal = dot(normal, ray.dir) < 0 ? normal : -normal;
        rec.p = ray.at(rec.t);
        rec.mat_ptr = mat_ptr;
        rec.obj = this;
        return true;
    }
    else
        return false;
}
bool Rect::sample(Hittable_list &scene, const Hit_record &rec, color &output)
{
    Vec3lf v0v1 = v[1] - v[0];
    Vec3lf v0v3 = v[3] - v[0];
    Vec3lf center = (v[0] + v[1] + v[2] + v[3]) / 4;

    Vec3lf pCenter = center - rec.p;
    Vec3lf pNormal = dot(normal, pCenter) < 0 ? normal : -normal;

    Vec3lf samplePoint(0, 0, 0);
    Ray shadowRay;
    uint32_t index;

    int samples = 0;
    Hit_record scene_rec;
    for (float x = 0; x <= 1; x += step)
    {
        for (float y = 0; y <= 1; y += step)
        {
            samples++;
            samplePoint = v[0] + v0v1 * x + v0v3 * y;
            shadowRay.orig = samplePoint;
            shadowRay.dir = rec.p - samplePoint;
            shadowRay.epsilon(pNormal);

            if (scene.hit(shadowRay, scene_rec, index))
            {
                if (scene_rec.obj != rec.obj)
                    continue;

                // face normal
                if (scene_rec.normal.x - rec.normal.x > EPSILON || scene_rec.normal.y - rec.normal.y > EPSILON || scene_rec.normal.z - rec.normal.z > EPSILON)
                    continue;

                float dist = (scene_rec.p - rec.p).norm();
                if (abs(dist) > EPSILON)
                    continue;
            }

            output += mat_ptr->Emission(shadowRay, rec);
        }
    }
    output = output / samples;
    return true;
}
void Rect::outputString()
{
    std::cerr << "Rect" << std::endl;
}
Bound3 Rect::getBounds()
{
    Vec3lf pMin = Vec3lf::Min(Vec3lf::Min(Vec3lf::Min(v[0], v[1]), v[2]), v[3]);
    Vec3lf pMax = Vec3lf::Max(Vec3lf::Max(Vec3lf::Max(v[0], v[1]), v[2]), v[3]);
    return Bound3(pMin, pMax);
}
Material *Rect::getMaterial()
{
    return mat_ptr;
}
#pragma endregion
/////////////////////////////////////////////////////////////////////////////////
#pragma region Disk

bool Disk::hit(const Ray &ray, Hit_record &rec)
{
    float t;
    bool inFlat;
    if (rayFlatIntersect3D(ray.orig, ray.dir, p, normal, t, inFlat))
    {
        if (inFlat)
        {
            float dist2;
            rayLineDist(ray.orig, ray.dir, p, t, dist2);
            if (dist2 <= r2)
            {
                float dt = sqrt(r2 - dist2) / ray.dir.norm();
                if (t + dt < 0)
                    return false;
                if (t - dt >= 0)
                    t = t - dt;
                else
                    t = t + dt;

                rec.happened = true;
                rec.t = t;
                rec.p = ray.at(t);
                rec.normal = dot(normal, ray.dir) < 0 ? normal : -normal;
                rec.obj = this;
                rec.mat_ptr = mat_ptr;
                return true;
            }
            else
            {
                return false;
            }
        }

        if ((ray.at(t) - p).norm2() <= r2)
        {

            rec.happened = true;
            rec.t = t;
            rec.p = ray.at(t);
            rec.normal = dot(normal, ray.dir) < 0 ? normal : -normal;
            rec.obj = this;
            rec.mat_ptr = mat_ptr;

            return true;
        }
        else
        {
            return false;
        }
    }
    return false;
}
bool Disk::sample(Hittable_list &scene, const Hit_record &rec, color &output)
{
    Point3lf samplePoint;
    Ray shadowRay;
    int samples = 0;

    Vec3lf pCenter = p - rec.p;
    Vec3lf pNormal = dot(normal, pCenter) < 0 ? normal : -normal;

    Hit_record scene_rec;
    uint32_t index;

    for (float i = 0; i <= r; i += step_r)
    {
        for (float j = 0.0; j <= 360.0; j += step_theta)
        {
            samples++;
            double rad = deg2rad(j);
            samplePoint = p + i * x * cos(rad) + i * y * sin(rad);
            shadowRay.orig = samplePoint;
            shadowRay.epsilon(pNormal);
            shadowRay.dir = rec.p - shadowRay.orig;

            if (scene.hit(shadowRay, scene_rec, index))
            {
                if (scene_rec.obj != rec.obj)
                    continue;

                // face normal
                if ((scene_rec.normal.normalize() - rec.normal.normalize()).norm2() > EPSILON)
                    continue;

                float dist = (scene_rec.p - rec.p).norm2();
                if (abs(dist) > EPSILON)
                    continue;
            }

            output += mat_ptr->Emission(shadowRay, rec);
        }
    }

    output = output / samples;
    return true;
}
void Disk::outputString()
{
    std::cerr << "Disk" << std::endl;
}
Bound3 Disk::getBounds()
{
    Vec3lf p1 = p + x + y;
    Vec3lf p2 = p - x - y;
    Vec3lf p3 = p + x - y;
    Vec3lf p4 = p - x + y;
    Vec3lf pMin = Vec3lf::Min(p1, p2);
    pMin = Vec3lf::Min(pMin, p3);
    pMin = Vec3lf::Min(pMin, p4);
    Vec3lf pMax = Vec3lf::Max(p1, p2);
    pMax = Vec3lf::Max(pMax, p3);
    pMax = Vec3lf::Max(pMax, p4);
    return Bound3(pMin, pMax);
}
bool Disk::PTSample(HittableSample &sam, float &pdf)
{
    /*
    float rho = r * random_float(0, 1);
    float the = 2.f * M_PI * random_float(0, 1);

    sam.pos = p + rho * x.normalize() * cos(the) + rho * y.normalize() * sin(the);
    sam.normal = normal;
    sam.isEmit = mat_ptr->hasEmission();
    sam.emit = mat_ptr->getEmission();

    pdf = 1.0f / area;
    */
    Vec2f ret;
    do
    {
        sampler->sample2D(ret);
        ret = ret - Vec2f(r, r);
    } while (ret.x * ret.x + ret.y * ret.y >= r2);

    sam.pos = p + x.normalize() * ret.x + y.normalize() * ret.y;
    sam.normal = normal;
    sam.isEmit = mat_ptr->hasEmission();
    sam.emit = mat_ptr->getEmission();

    pdf = 1.0f / area;
    return true;
}
#pragma endregion
/////////////////////////////////////////////////////////////////////////////////
#pragma region Cylinder

bool Cylinder::hit(const Ray &r_in, Hit_record &rec)
{
    Ray r_mov;
    r_mov.orig = r_in.orig - Vec3lf(p.x, 0, p.z);
    r_mov.dir = r_in.dir;
    Point3lf temp_p;

    /*
        float a = r_in.dir.x * r_in.dir.x + r_in.dir.z * r_in.dir.z;
        float b = 2 * (r_in.dir.x * r_in.orig.x + r_in.dir.z * r_in.orig.z - p.x * r_in.dir.x - p.z * r_in.dir.z);
        float c = r_in.orig.x * r_in.orig.x + r_in.orig.z * r_in.orig.z + p.x * p.x + p.z * p.z - 2 * r_in.orig.x * p.x - 2 * r_in.orig.z * p.z - r2;
    */

    float a = r_mov.dir.x * r_mov.dir.x + r_mov.dir.z * r_mov.dir.z;
    float b = 2 * (r_mov.dir.x * r_mov.orig.x + r_mov.dir.z * r_mov.orig.z);
    float c = r_mov.orig.x * r_mov.orig.x + r_mov.orig.z * r_mov.orig.z - r2;

    float t0, t1;

    if (!solveQuadratic(a, b, c, t0, t1))
        return false;
    if (t0 == t1)
        return false;

    temp_p = r_in.at(t0);
    if ((t0 < 0) || (temp_p.y > hmax) || (temp_p.y < hmin))
    {
        t0 = t1;
        temp_p = r_in.at(t0);
        if ((t0 < 0) || (temp_p.y > hmax) || (temp_p.y < hmin))
            return false;
    }

    // rec
    rec.happened = true;
    rec.t = t0;
    rec.p = r_in.at(t0);
    rec.obj = this;
    rec.mat_ptr = mat_ptr;
    rec.set_face_normal(r_in, (rec.p - Vec3lf(p.x, rec.p.y, p.z)).normalize());

    return true;
}

void Cylinder::outputString()
{
    std::cerr << "Cylinder" << std::endl;
}
Bound3 Cylinder::getBounds()
{
    Vec3lf p1 = p + Vec3lf(-r, hmin, -r);
    Vec3lf p2 = p + Vec3lf(r, hmax, r);
    return Bound3(p1, p2);
}
#pragma endregion

#pragma region Polygon

Polygon::Polygon(const std::vector<Point3lf> &p_list, Material *_m, Sampler2D *_sampler)
{
    Vec3lf min_vert = Vec3lf{std::numeric_limits<float>::infinity(),
                             std::numeric_limits<float>::infinity(),
                             std::numeric_limits<float>::infinity()};
    Vec3lf max_vert = Vec3lf{-std::numeric_limits<float>::infinity(),
                             -std::numeric_limits<float>::infinity(),
                             -std::numeric_limits<float>::infinity()};

    for (Point3lf p : p_list)
    {
        min_vert = Vec3lf(std::min(min_vert.x, p[0]),
                          std::min(min_vert.y, p[1]),
                          std::min(min_vert.z, p[2]));
        max_vert = Vec3lf(std::max(max_vert.x, p[0]),
                          std::max(max_vert.y, p[1]),
                          std::max(max_vert.z, p[2]));
    }
    bounding_box = Bound3(min_vert, max_vert);
    mat_ptr = _m;
    sampler = _sampler;
    area = 0;
    triangles.clear();
    // 划分三角形012,023,034...
    for (int i = 1; i < p_list.size() - 1; i++)
    {
        Triangle tri(p_list[0], p_list[i], p_list[i + 1], _m, sampler);
        area += tri.area;
        triangles.push_back(tri);
    }
}

bool Polygon::hit(const Ray &ray, Hit_record &rec)
{
    for(Triangle tri: triangles)
    {
        if(tri.hit(ray, rec) == true)
            return true;
    }
    return false;
}
void Polygon::outputString()
{
    std::cerr << "Polygon" << std::endl;
}
bool Polygon::PTSample(HittableSample &sam, float &pdf)
{
    float p = random_float() * area;
    float emit_area_sum = 0;
    for(Triangle tri: triangles)
    {
        emit_area_sum += tri.area;
        if(p <= emit_area_sum)
        {
            tri.PTSample(sam, pdf);
            return true;
        }
    }
    triangles[triangles.size() - 1].PTSample(sam, pdf);
    return true;
}

Bound3 Polygon::getBounds()
{
    return bounding_box;
}

float Polygon::getArea()
{
    return area;
}
#pragma endregion