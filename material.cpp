#include "material.hpp"
#include "Texture.hpp"
#include "hittable_list.hpp"
#include <algorithm>
#include <tuple>
#include <cmath>

#pragma region Material

bool Material::reflect(
    const Ray &r_in, const Hit_record &rec, Ray &scattered) const
{
    return false;
}
bool Material::refract(
    const Ray &r_in, const Hit_record &rec, Ray &scattered) const
{
    return false;
}
bool Material::get_color(
    const Ray &r_in, const Hit_record &rec, Hittable_list &scene, color &c) const
{
    return false;
}
float Material::get_ior()
{
    return 0;
}
color Material::Emission(const Ray &r, const Hit_record &rec) const
{
    return color(0, 0, 0);
}
bool Material::hasEmission() const
{
    return isEmission;
}
color Material::getEmission() const
{
    return m_Emission;
}
color Material::BRDF(const Vec3lf &wi, const Vec3lf &wo, const Hit_record &hit_rec)
{
    return color(0, 0, 0);
}
bool Material::sampleScatter(const Ray &r_in, const Hit_record &rec, Ray &r_out,float& pdf)
{
    pdf = 0;
    return false;
}
float Material::pdf(const Vec3lf &wi, const Vec3lf &wo, const Vec3lf &N)
{
    return 0;
}
#pragma endregion

///////////////////////////////////////////////////////////////////////////////////////////////////

#pragma region Light

bool Light::get_color(
    const Ray &r_in, const Hit_record &rec, Hittable_list &scene, color &c) const
{
    c = intensity;
    return true;
}

color Light::Emission(const Ray &r, const Hit_record &rec) const
{
    Vec3lf v = rec.p - r.orig;
    float dist2 = v.norm2();

    return dist2 <= 1 ? intensity : intensity / dist2;
}

#pragma endregion

///////////////////////////////////////////////////////////////////////////////////////////////////

#pragma region Diffuse

bool Diffuse::get_color(
    const Ray &r_in, const Hit_record &rec, Hittable_list &scene, color &c) const
{
    int count = 0;
    color temp_c;
    if (texture != nullptr && rec.uv == true)
    {
        Vec2lf uv;
        if (rec.obj->get_uv(rec, uv))
        {
            unsigned char r, g, b;
            texture->load_xy(uv.x, uv.y, r, g, b);
            c = Vec3lf((float)r / 256.0, (float)g / 256.0, (float)b / 256.0);
            // [r,g,b] = texture->load_xy(uv.x, uv.y);
        }
        else
            c = diffuse_color;
    }
    else
        c = diffuse_color;

    Material *mat;
    for (Hittable *object : scene.objects)
    {
        mat = object->getMaterial();
        if (mat)
        {
            if (mat == nullptr)
            {
                std::cerr << "Error: Can't find the Material of ";
                object->outputString();
            }
            if (mat->isEmission == true)
            {
                if (object->sample(scene, rec, temp_c))
                {
                    c *= temp_c;
                }
            }
        }
    }

    return true;
}

#pragma endregion

///////////////////////////////////////////////////////////////////////////////////////////////////

#pragma region Specular

bool Specular::reflect(
    const Ray &r_in, const Hit_record &rec, Ray &scattered) const
{
    if (canReflect == false)
        return false;
    Vec3lf dir = vec_reflect(r_in.direction().normalize(), rec.normal);
    scattered = Ray(rec.p, dir);
    scattered.epsilon(rec.normal);
    return true;
}

bool Specular::refract(
    const Ray &r_in, const Hit_record &rec, Ray &scattered) const
{
    if (canRefract == false)
        return false;
    Vec3lf dir = vec_refract(r_in.direction().normalize(), rec.normal, ior);
    scattered = Ray(rec.p, dir);
    scattered.canReflect = false;
    scattered.epsilon(rec.normal);
    return true;
}

float Specular::get_ior()
{
    return ior;
}

#pragma endregion

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma region Lambertian

color Lambertian::BRDF(const Vec3lf &wi, const Vec3lf &wo, const Hit_record &hit_rec)
{
    float cosalpha = dot(wo, hit_rec.normal);
    if (cosalpha > 0.0f)
    {
        Vec3lf diffuse = kd / M_PI;
        return diffuse;
    }
    else
        return Vec3lf(0.f, 0.f, 0.f);
}
bool Lambertian::sampleScatter(const Ray &r_in, const Hit_record &rec, Ray &r_out,float& pdf)
{
    /*
    Vec3lf randVec = random_in_unit_sphere();
    if(randVec.y < 0) randVec.y = -randVec.y;

    r_out.orig = rec.p;
    r_out.dir = toWorld(randVec, rec.normal);
    r_out.epsilon(rec.normal);
    return true;
    */
    float x_1 = random_float(), x_2 = random_float();
    float z = std::fabs(1.0f - 2.0f * x_1);
    float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
    Vec3lf localRay(r * std::cos(phi), r * std::sin(phi), z);
    r_out.orig = rec.p;
    r_out.dir = toWorld(localRay, rec.normal).normalize();
    r_out.epsilon(rec.normal);

    pdf = this->pdf(r_in.dir, r_out.dir, rec.normal);
    return true;
}
float Lambertian::pdf(const Vec3lf &wi, const Vec3lf &wo, const Vec3lf &N)
{
    if (dot(wo, N) > 0.f)
        return 0.5f / M_PI;
    else
        return 0.f;
}

#pragma endregion

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma region CoookTorrance

color CookTorrance::BRDF(const Vec3lf &wi, const Vec3lf &wo, const Hit_record &hit_rec)
{
    Vec3f h = (wi + wo).normalize();
    float LoH = fabs(dot(wi, h));
    float VoH = fabs(dot(wo, h));

    Vec3f F0 = lerp(baseF0, c, metallic);
    Vec3f F = FresnelSchlick(F0, VoH);

    // Vec3f kd = Vec3f(1, 1, 1);
    Vec3f kd = (1.0f - metallic) * (Vec3f(1, 1, 1) - F);
    Vec3f diffuse = dotMul(c, kd) / M_PI;

    float D = NDF(hit_rec.normal, h, roughness);
    float G = GF(hit_rec.normal, wi, wo);
    float LoN = fabs(dot(wi, hit_rec.normal));
    float VoN = fabs(dot(wo, hit_rec.normal));
    Vec3f ks = Vec3f(1, 1, 1);// - kd;
    Vec3f specular = dotMul(ks, F) * D * G / (4 * LoN * VoN);

    return diffuse + specular;
}

bool CookTorrance::sampleScatter(const Ray &r_in, const Hit_record &rec, Ray &r_out, float& pdf)
{
    float x_1 = random_float(), x_2 = random_float();

    float phi = 2 * M_PI * x_1;
    float a2 = roughness * roughness;
    float cosTheta = sqrt((1 - x_2) / (1 + (a2 - 1) * x_2));
    float sinTheta = sqrt(1 - cosTheta * cosTheta);

    float x = sinTheta * cos(phi);
    float y = sinTheta * sin(phi);
    float z = cosTheta;

    Vec3f h(x, y, z);
    h = toWorld(h.normalize(), rec.normal);

    r_out.orig = rec.p;
    r_out.dir = vec_reflect(-r_in.orig.normalize(), h);
    r_out.epsilon(rec.normal);

    float d = (cosTheta * a2 - cosTheta) * cosTheta + 1;
    float D = a2 / (M_PI * d * d);
    pdf = D * cosTheta;

    return true;
}

float CookTorrance::pdf(const Vec3lf &wi, const Vec3lf &wo, const Vec3lf &N)
{

    Vec3f h = (wi + wo).normalize();
    float NoH = fabs(dot(N, h));
    float VoH = fabs(dot(wo, h));

    float a2 = roughness * roughness;
    float x = (a2 - 1.0f) * NoH * NoH + 1.0f;
    // D * NoH / (4 * VoH)
    return a2 * NoH / (M_PI * x * x) / (4 * VoH);
}

Vec3f CookTorrance::FresnelSchlick(Vec3f F0, float VoH) const
{
    return lerp(F0, Vec3f(1, 1, 1), pow(1 - VoH, 5));
}
float CookTorrance::NDF(const Vec3f &n, const Vec3f &h, float alpha) const
{
    float NoH = fabs(dot(n, h));
    float a2 = alpha * alpha;
    float x = NoH * NoH * (a2 - 1.0f) + 1.0f;
    return a2 / (M_PI * x * x);
}
float CookTorrance::GF(const Vec3f &n, const Vec3f &wi, const Vec3f &wo) const
{
    float NoL = fabs(dot(n, wi));
    float NoV = fabs(dot(n, wo));
    float k = (roughness + 1.0f) * (roughness + 1.0f) / 8.0f;
    return GGXschlick(NoL, k) * GGXschlick(NoV, k);
}

float CookTorrance::GGXschlick(float NoV, float k) const
{
    return NoV / (lerp(k, 1, NoV));
}

#pragma endregion