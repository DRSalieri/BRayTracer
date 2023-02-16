#ifndef __RENDERER_H__
#define __RENDERER_H__

#include "camera.hpp"
#include "hittable_list.hpp"

// Image
const float aspect_ratio = 3.0 / 2.0;
const int image_width = 300;
const int image_height = static_cast<int>(image_width / aspect_ratio);
const int samples_per_pixel = 20;
const float spp_dist = 1.f / (float) samples_per_pixel;
const int max_depth = 20;
const float RR_prob = 0.8;
/*
const float ks = 0.7;
const float kr = 0.3;
*/
class Renderer
{
    public:
    void Render(Hittable_list& scene, const Camera& cam);

};

#endif