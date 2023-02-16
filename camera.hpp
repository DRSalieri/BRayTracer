#ifndef __CAMERA_H__
#define __CAMERA_H__

#include "geometry.hpp"
#include "ray.hpp"

class Camera
{
public:
    Camera(Vec3lf lookfrom, Vec3lf lookat, Vec3lf vup, float vfov, float aspect)
    {
        Vec3lf u, v, w;

        float theta = deg2rad(vfov);
        float half_height = tan(theta / 2);
        float half_width = aspect * half_height;

        origin = lookfrom;
        w = (lookfrom - lookat).normalize();
        u = cross(vup, w).normalize();
        v = cross(w, u);

        lower_left_corner = origin - half_width * u - half_height * v - w;
        horizontal = 2 * half_width * u;
        Vertical = 2 * half_height * v;
    }

public:
    Ray get_ray (float s, float t) const
    {
        return Ray(origin, lower_left_corner + s * horizontal + t * Vertical - origin);
    }

public:
    Point3lf origin;
    Point3lf lower_left_corner;
    Vec3lf horizontal;
    Vec3lf Vertical;
};

/*
class Camera
{
public:
    Camera(
        Point3lf lookfrom,
        Point3lf lookat,
        Vec3lf vup,
        float vfov, // vertical field-of-view in degrees
        float aspect_ratio,
        float aperture,
        float focus_dist)
    {
        auto theta = degrees_to_radians(vfov);
        auto h = tan(theta / 2);
        auto viewport_height = 2.0 * h;
        auto viewport_width = aspect_ratio * viewport_height;

        w = (lookfrom - lookat).normalize();
        u = cross(vup, w).normalize();
        v = cross(w, u);

        origin = lookfrom;
        horizontal = focus_dist * viewport_width * u;
        vertical = focus_dist * viewport_height * v;
        lower_left_corner = origin - horizontal / 2 - vertical / 2 - focus_dist * w;

        lens_radius = aperture / 2;
    }

    Ray get_ray(float s, float t) const
    {
        // Vec3lf rd = lens_radius * random_in_unit_disk();
        // Vec3lf offset = u * rd.x + v * rd.y;
        Vec3lf offset = u + v;
        return Ray(
            origin + offset,
            lower_left_corner + s * horizontal + t * vertical - origin - offset);
    }

private:
    Point3lf origin;
    Point3lf lower_left_corner;
    Vec3lf horizontal;
    Vec3lf vertical;
    Vec3lf u, v, w;
    float lens_radius;
};
*/
#endif