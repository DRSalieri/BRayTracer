#include "Renderer.hpp"
#include "hittable_list.hpp"
#include "output.hpp"

#define __PATH_TRACING__

int times = 0;
void sampleLight(Hittable_list& objects, HittableSample &sam, float &pdf)
{
    // sum of surface area
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.objects.size(); ++k)
    {
        if (objects.objects[k]->mat_ptr->hasEmission())
        {
            emit_area_sum += objects.objects[k]->getArea();
        }
    }

    float p = random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.objects.size(); ++k)
    {
        if (objects.objects[k]->mat_ptr->hasEmission())
        {
            emit_area_sum += objects.objects[k]->getArea();
            if (p <= emit_area_sum)
            {
                objects.objects[k]->PTSample(sam, pdf);
                break;
            }
        }
    }
}
/***
 * r - (eye -> viewport)
 * scene - ( hittable_list )
 * return the color
 **/
color RayCast(const Ray &r, Hittable_list &scene, int depth)
{
    #ifdef __PATH_TRACING__

    Vec3lf color_dir = Vec3lf(0,0,0);
    Vec3lf color_indir = Vec3lf(0,0,0);

    Hit_record rec;
    uint32_t index;
    scene.hit(r, rec, index);
    if(!rec.happened){
        return Vec3lf(0,0,0);
    }

    // emission
    if(rec.mat_ptr->hasEmission())
        return rec.mat_ptr->getEmission();

    // direct light
    HittableSample hit_sam;
    float light_pdf = 0.0f;
    sampleLight(scene, hit_sam, light_pdf);

    Vec3lf shadowDir = hit_sam.pos - rec.p;
    float dis = shadowDir.norm2();
    shadowDir = shadowDir.normalize();
    Ray shadowRay(rec.p, shadowDir);
    shadowRay.epsilon(rec.normal);
    
    hit_sam.normal = dot(hit_sam.normal, shadowRay.dir) < 0 ? hit_sam.normal : - hit_sam.normal; 

    Hit_record shadow_rec;
    scene.hit(shadowRay, shadow_rec, index);
    // judge the distance
    if( abs( (shadow_rec.p - hit_sam.pos).norm() ) < EPSILON )
    {
       Vec3lf f_r = rec.mat_ptr->BRDF( -r.dir.normalize(), shadowRay.dir.normalize(), rec);
        color_dir = dotMul(hit_sam.emit, f_r);
        // 需要判断当dis、pdf过小的情况?
        float cosA = std::max(.0f, dot(shadowRay.dir.normalize(), rec.normal));
        float cosB = std::max(.0f, dot(-shadowRay.dir.normalize(), hit_sam.normal));
        //float brdf_pdf = rec.mat_ptr->pdf(-r.dir.normalize(), shadowRay.dir.normalize(), rec.normal);
        // if(pdf2 >= EPSILON)

        if( light_pdf > 0)
            color_dir = color_dir * cosA * cosB / dis / light_pdf;
        else
            color_dir = color(0,0,0);
    }

    if( depth <= 0 || random_float() > RR_prob )
        return color_dir;

    // indirect light
    Ray scatterRay;
    float scatterPdf;
    rec.mat_ptr->sampleScatter(r, rec, scatterRay, scatterPdf);
    Hit_record scatterRec;
    scene.hit(scatterRay, scatterRec, index);
    if(scatterRec.happened)
    {
        // if(scatterRec.mat_ptr->hasEmission())
        //   return color_dir;
        Vec3lf f_r = rec.mat_ptr->BRDF(-r.dir.normalize(), scatterRay.dir.normalize(),rec);
        float cos = std::max(.0f, dot(scatterRay.dir.normalize(), rec.normal));
        float brdf_pdf = rec.mat_ptr->pdf(-r.dir.normalize(), scatterRay.dir.normalize(), rec.normal);
        color_indir = dotMul( RayCast(scatterRay, scene, depth - 1) , f_r);

        if(brdf_pdf > 0)
            color_indir = color_indir * cos / brdf_pdf / RR_prob;
        else
            color_indir = color(0,0,0);
    }

    color ret = color_indir * 0.5 + color_dir * 0.5;

    return ret;


    #else
    Hit_record rec;
    if (depth <= 0)
        return color(0, 0, 0);
    uint32_t index;
    Vec2lf uv;

    if (scene.hit(r, rec, index))
    {
        color c(0, 0, 0);

        if (rec.mat_ptr->isColor)
        {
            rec.mat_ptr->get_color(r, rec, scene, c);

            return c;
        }

        if (rec.mat_ptr->isRecursive)
        {
            float kr;
            fresnel(r.dir, rec.normal, rec.mat_ptr->get_ior(), kr);

            Ray reflectRay, refractRay;
            if (rec.mat_ptr->refract(r, rec, refractRay))
            {
                c += RayCast(refractRay, scene, depth - 1) * (1 - kr);
            } else {
                kr = 1;
            }
            if (r.canReflect && rec.mat_ptr->reflect(r, rec, reflectRay))
            {
                c += RayCast(reflectRay, scene, depth - 1) * kr;
            }
            return c;
        }
    }

    return Vec3lf(0, 0, 0);
    #endif
}

void Renderer::Render(Hittable_list &scene, const Camera &cam)
{

    std::cout << "P3\n"
              << image_width << ' ' << image_height << "\n255\n";

    for (float j = image_height - 1; j >= 0; --j)
    {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (float i = 0; i < image_width; ++i)
        {

            color pixel_color(0, 0, 0);

            for (int s = 0; s < samples_per_pixel; ++s)
            {
                for(int t = 0;t < samples_per_pixel; ++t){
                    float u = ( i + spp_dist * s + random_float(0, spp_dist) ) ;
                    u /= image_width;
                    float v = ( j + spp_dist * t + random_float(0, spp_dist) ) ;
                    v /= image_height;
                    Ray r = cam.get_ray(u, v);
                    pixel_color += RayCast(r, scene, max_depth);
                }
            }
            write_color(std::cout, pixel_color, samples_per_pixel * samples_per_pixel);

            /*
            if(pixel_color.x > EPSILON)
            {
                printf("1");
            }
            */
        }
    }

    // rintf("\n%d\n",times);
}
