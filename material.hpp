#ifndef __MATERIAL_H__
#define __MATERIAL_H__

class Hit_record;
class Hittable_list;
class Light_list;

#include "Texture.hpp"
#include "ray.hpp"
#include "hittable.hpp"
#include "material.hpp"

class Material
{
public:
    /***
     *  in this function, you need to:
     * 1. (&scattered)  produce scattered ray
     * 2. (&color)      compute attenuation
     * 3. (return)      judge is scattered?
     ***/
    virtual bool reflect(
        const Ray &r_in, const Hit_record &rec, Ray &scattered) const;
    virtual bool refract(
        const Ray &r_in, const Hit_record &rec, Ray &scattered) const;
    virtual bool get_color(
        const Ray &r_in, const Hit_record &rec, Hittable_list &scene, color &color) const;
    virtual float get_ior();
    virtual color Emission(const Ray &r, const Hit_record &rec) const;

    // Path Tracing
    virtual bool hasEmission() const;
    virtual color getEmission() const;
    virtual color BRDF(const Vec3lf &wi, const Vec3lf &wo, const Hit_record &hit_rec);
    virtual bool sampleScatter(const Ray &r_in, const Hit_record &rec, Ray &r_out, float& pdf);
    virtual float pdf(const Vec3lf &wi, const Vec3lf &wo, const Vec3lf &N);
    // pdf?
public:
    bool isRecursive;
    bool isColor;

    // Path Tracing
    bool isEmission;
    color m_Emission = color(0, 0, 0); // emission color
};

class Light : public Material
{
public:
    Light(color _c) : intensity(_c)
    {
        isEmission = true;
        isRecursive = false;
        isColor = true;
    };

public:
    virtual bool get_color(
        const Ray &r_in, const Hit_record &rec, Hittable_list &scene, color &color) const override;
    virtual color Emission(const Ray &r, const Hit_record &rec) const override;

public:
    color intensity;
};

// Blinn-Phong
class Diffuse : public Material
{
public:
    Diffuse()
    {
        texture = nullptr;
        diffuse_color = Vec3lf(0, 0, 0);
        isEmission = false;
        isRecursive = false;
        isColor = true;
    }
    Diffuse(color _c) : diffuse_color(_c)
    {
        texture = nullptr;
        isEmission = false;
        isRecursive = false;
        isColor = true;
    };

public:
    virtual bool get_color(
        const Ray &r_in, const Hit_record &rec, Hittable_list &scene, color &color) const;

public:
    Texture *texture;
    color diffuse_color;
};

class Specular : public Material
{
public:
    Specular(float _ior, bool _refr, bool _refl) : ior(_ior), canReflect(_refl), canRefract(_refr)
    {
        isEmission = false;
        isRecursive = true;
        isColor = false;
    };
    virtual bool reflect(
        const Ray &r_in, const Hit_record &rec, Ray &scattered) const override;
    virtual bool refract(
        const Ray &r_in, const Hit_record &rec, Ray &scattered) const override;
    virtual float get_ior() override;

public:
    float ior;
    bool canRefract;
    bool canReflect;
};

// Path Tracing
class Lambertian : public Material
{
public:
    Lambertian(Vec3lf _kd, bool _isEmission, color _m_Emission)
    {
        kd = _kd;
        isEmission = _isEmission;
        m_Emission = _m_Emission;
    }

public:
    // Path Tracing
    virtual color BRDF(const Vec3lf &wi, const Vec3lf &wo,  const Hit_record &hit_rec) override;
    virtual bool sampleScatter(const Ray &r_in, const Hit_record &rec, Ray &r_out,float& pdf) override;
    virtual float pdf(const Vec3lf &wi, const Vec3lf &wo, const Vec3lf &N);

public:
    Vec3lf kd;
};

class CookTorrance : public Material
{
public:
    CookTorrance(color _c, Vec3f _baseF0, float _metallic, float _roughness, bool _isEmission, color _m_Emission):
        c(_c), baseF0(_baseF0), metallic(_metallic), roughness(_roughness){
            isEmission = _isEmission;
            m_Emission = _m_Emission;
        };
    ~CookTorrance(){};
public:
    virtual color BRDF(const Vec3lf &wi, const Vec3lf &wo, const Hit_record &hit_rec) override;
    virtual bool sampleScatter(const Ray &r_in, const Hit_record &rec, Ray &r_out,float& pdf) override;
    virtual float pdf(const Vec3lf &wi, const Vec3lf &wo, const Vec3lf &N);
private:
    Vec3f FresnelSchlick(Vec3f F0, float VoH) const;
    float NDF(const Vec3f& n, const Vec3f& h, float alpha) const;
    float GF(const Vec3f& n, const Vec3f& wi, const Vec3f& wo) const;
    float GGXschlick(float NoV, float k) const;
public:
    color c;
    Vec3f baseF0;
    float metallic, roughness;
};

#endif