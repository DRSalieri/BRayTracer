#include "geometry.hpp"

Vec3lf toWorld(const Vec3lf& a, const Vec3lf& N){
        Vec3lf B, C;
        if (std::fabs(N.x) > std::fabs(N.y)){
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            C = Vec3lf(N.z * invLen, 0.0f, -N.x *invLen);
        }
        else {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C = Vec3lf(0.0f, N.z * invLen, -N.y *invLen);
        }
        B = cross(C, N);
        return a.x * B + a.y * C + a.z * N;
}

// interpolate
Vec3lf interpolate(float alpha, float gamma, float beta, Vec3lf &vert1, Vec3lf &vert2, Vec3lf &vert3, float weight = 1.0)
{
    return (alpha * vert1 + beta * vert2 + gamma * vert3) / weight;
}

Vec2lf interpolate(float alpha, float gamma, float beta, Vec2lf &vert1, Vec2lf &vert2, Vec2lf &vert3, float weight = 1.0)
{
    return (alpha * vert1 + beta * vert2 + gamma * vert3) / weight;
}

std::tuple<float, float, float> computeBarycentric2D(float x, float y, Vec2lf &vert1, Vec2lf &vert2, Vec2lf &vert3)
{
    float alpha = (x * (vert2.y - vert3.y) + y * (vert3.x - vert2.x) + vert2.x * vert3.y - vert2.y * vert3.x) / (vert1.x * (vert2.y - vert3.y) + vert1.y * (vert3.x - vert2.x) + vert2.x * vert3.y - vert2.y * vert3.x);

    float beta = (x * (vert3.y - vert1.y) + y * (vert1.x - vert3.x) + vert3.x * vert1.y - vert3.y * vert1.x) / (vert2.x * (vert3.y - vert1.y) + vert2.y * (vert1.x - vert3.x) + vert3.x * vert1.y - vert3.y * vert1.x);

    float gamma = (x * (vert1.y - vert2.y) + y * (vert2.x - vert1.x) + vert1.x * vert2.y - vert1.y * vert2.x) / (vert3.x * (vert1.y - vert2.y) + vert3.y * (vert2.x - vert1.x) + vert1.x * vert2.y - vert1.y * vert2.x);

    return std::tuple<float, float, float>(alpha, beta, gamma);
}

Vec3lf random_Vec3lf() {
    Vec3lf ret;
    for (size_t i=3; i--; ret[i] = random_float());
    return ret;
}

Vec3lf random_in_unit_sphere() {
    while (true) {
        Vec3lf p = random_Vec3lf();
        if (length_squared(p) >= 1) continue;
        return p;
    }
}

Vec3lf random_unit_vector() {
    return random_in_unit_sphere().normalize();
}

Vec3lf random_in_unit_disk() {
    while (true) {
        Vec3lf p = Vec3lf(random_float(-1,1), random_float(-1,1), 0);
        if (length_squared(p) >= 1) continue;
        return p;
    }
}

Vec3lf vec_refract(const Vec3lf &I, const Vec3lf &N, const float &ior)
{
    float cosi = clamp(-1, 1, dot(I, N));
    float etai = 1, etat = ior;
    Vec3lf n = N;
    if (cosi < 0)
    {
        cosi = -cosi;
    }
    else
    {
        std::swap(etai, etat);
        n = -N;
    }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? Vec3lf(0, 0, 0) : eta * I + (eta * cosi - sqrtf(k)) * n;
}

void fresnel(const Vec3lf &I, const Vec3lf &N, const float &ior, float &kr)
{
    float cosi = clamp(-1, 1, dot(I, N));
    float etai = 1, etat = ior;
    if (cosi > 0)
    {
        std::swap(etai, etat);
    }
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1)
    {
        kr = 1;
    }
    else
    {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        kr = (Rs * Rs + Rp * Rp) / 2;
    }
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
}

float lerp(float a, float b, float t)
{
    return a + t * (b - a);
}

Vec3f lerp(const Vec3f &a, const Vec3f &b, float t)
{
    return a + t * (b - a);
}