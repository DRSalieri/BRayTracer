#include "geometry.hpp"
#include "material.hpp"
#include "camera.hpp"
#include "hittable_list.hpp"
#include "Renderer.hpp"
#include "Transform.hpp"
#include <iostream>
#include <cstring>

// #define BRIDSON_OUTPUT

Hittable_list scene_4()
{
    Hittable_list scene;

    // floor
    Material *ground_material1 = new Diffuse(color(0.5, 0.5, 0.5));
    scene.add(new Cube(Point3lf(-100, -1, -100), Point3lf(100, 0, 100), ground_material1));

    // back wall
    Material *ground_material2 = new Diffuse(color(0.5, 0.5, 0.5));
    scene.add(new Cube(Point3lf(4, 0, -2), Point3lf(5, 2, 2), ground_material2));

    // left wall
    Material *material1 = new Diffuse(color(0.9, 0.1, 0.1));
    scene.add(new Cube(Point3lf(0, 0, -100), Point3lf(4, 2, -1), material1));

    // right wall
    Material *material2 = new Diffuse(color(0.1, 0.9, 0.1));
    scene.add(new Cube(Point3lf(0, 0, 1), Point3lf(4, 2, 100), material2));

    // ceil
    Material *material3 = new Diffuse(color(0.1, 0.1, 0.9));
    Point3lf v0 = Point3lf(0, 2, -1);
    Point3lf v1 = Point3lf(0, 2, 1);
    Point3lf v2 = Point3lf(4, 2, 1);
    Point3lf v3 = Point3lf(4, 2, -1);
    scene.add(new Rect(v0, v1, v2, v3, material3));

    Material *light = new Light(color(1, 1, 1));

    /*
    v0 = Point3lf(2.8, 1.9999, -0.2);
    v1 = Point3lf(2.8, 1.9999, 0.2);
    v2 = Point3lf(3.2, 1.9999, 0.2);
    v3 = Point3lf(3.2, 1.9999, -0.2);
    scene.add(new Rect(v0, v1, v2, v3, light, 0.04));
    */

    Point3lf p = Point3lf(1.95, 1.98, 0.4);
    float r = 0.2;
    Vec3lf pNormal = Vec3lf(0, -1, 0);
    Vec3lf px = Vec3lf(1, 0, 0);
    scene.add(new Disk(p, px, pNormal, r, light, 0.02, 30));

    // scene.add(new Cylinder(Point3lf(2, 0, -0.6), 0.3, 0.5, material3));
    //  scene.add(new Disk(Point3lf(2,0.5,-0.6),Vec3lf(1,0,0),Vec3lf(0,1,0), 0.3, material3));
    // scene.add(new Cube(Point3lf(1.4, 0, -0.6), Point3lf(1.8, 0.6, -0.2), material3));

    Matrix mat = ModelTransform(Vec3lf(3, 0, -0.4), Vec3lf(0, 90, 0), Vec3lf(10, 10, 10));
    // scene.add(new MeshTriangle("./obj/B/B.obj",material1,mat));

    Material *ball = new Specular(1.8, false, true);
    scene.add(new Sphere(Point3lf(1, 0.4, 0.5), 0.4, ball));

    return scene;
}

Hittable_list scene_1()
{
    Hittable_list scene;

    // floor
    Material *ground_material1 = new Diffuse(color(0.5, 0.5, 0.5));
    scene.add(new Cube(Point3lf(-400, -1, -400), Point3lf(400, 0, 400), ground_material1));
    // back wall
    // scene.add(new Cube(Point3lf(200, 0, -200), Point3lf(210, 200, 200), ground_material1));
    // ceil
    scene.add(new Cube(Point3lf(0, 200, -200), Point3lf(200, 210, 200), ground_material1));

    // left wall
    Material *material1 = new Diffuse(color(0.9, 0.1, 0.1));
    scene.add(new Cube(Point3lf(0, 0, -150), Point3lf(200, 200, -140), material1));

    // right wall
    Material *material2 = new Diffuse(color(0.1, 0.9, 0.1));
    scene.add(new Cube(Point3lf(0, 0, 140), Point3lf(200, 200, 150), material2));

    Material *light = new Light(color(1000, 1000, 1000));

    Point3lf p = Point3lf(100, 199, 0);
    float r = 10;
    Vec3lf pNormal = Vec3lf(0, -1, 0);
    Vec3lf px = Vec3lf(1, 0, 0);
    scene.add(new Disk(p, px, pNormal, r, light, 1, 30));

    /*
    Matrix mat = ModelTransform(Vec3lf(130, 0, -50), Vec3lf(0, 90, 0), Vec3lf(60, 60, 60));
    Material *specular = new Specular(1.8, false, true);
    scene.add(new MeshTriangle("./obj/bunny.obj", "", specular, mat, false));
    */

    // Material *ball = new Specular(1.8, false, true);
    // scene.add(new Sphere(Point3lf(1, 0.4, 0.5), 0.4, ball));

    return scene;
}

Hittable_list scene_2()
{
    Hittable_list scene;
    
    Material *m1 = new CookTorrance(color(0.65f, 0.05f, 0.05f), Vec3f(0.6f, 0.01f, 0.01f) , 0.1, 0.9, false, color(0,0,0));
    Material *m2 = new CookTorrance(color(0.65f, 0.05f, 0.05f), Vec3f(0.6f, 0.01f, 0.01f) , 0.5, 0.9, false, color(0,0,0));
    Material *m3 = new CookTorrance(color(0.65f, 0.05f, 0.05f), Vec3f(0.6f, 0.01f, 0.01f) , 0.9, 0.9, false, color(0,0,0));
    Material *m4 = new CookTorrance(color(0.65f, 0.05f, 0.05f), Vec3f(0.6f, 0.01f, 0.01f) , 0.1, 0.5, false, color(0,0,0));
    Material *m5 = new CookTorrance(color(0.65f, 0.05f, 0.05f), Vec3f(0.6f, 0.01f, 0.01f) , 0.5, 0.5, false, color(0,0,0));
    Material *m6 = new CookTorrance(color(0.65f, 0.05f, 0.05f), Vec3f(0.6f, 0.01f, 0.01f) , 0.9, 0.5, false, color(0,0,0));
    Material *m7 = new CookTorrance(color(0.65f, 0.05f, 0.05f), Vec3f(0.6f, 0.01f, 0.01f) , 0.1, 0.1, false, color(0,0,0));
    Material *m8 = new CookTorrance(color(0.65f, 0.05f, 0.05f), Vec3f(0.6f, 0.01f, 0.01f) , 0.5, 0.1, false, color(0,0,0));
    Material *m9 = new CookTorrance(color(0.65f, 0.05f, 0.05f), Vec3f(0.6f, 0.01f, 0.01f) , 0.9, 0.1, false, color(0,0,0));
    // floor
    Material *ground_material1 = new Lambertian(color(0.5, 0.5, 0.5), false, color(0, 0, 0));
    scene.add(new Cube(Point3lf(-400, -1, -400), Point3lf(400, 0, 400), ground_material1));

    // back wall
    scene.add(new Cube(Point3lf(300, 0, -200), Point3lf(310, 200, 200), ground_material1));
    // ceil
    scene.add(new Cube(Point3lf(0, 200, -200), Point3lf(300, 210, 200), ground_material1));

    // left wall
    
    Material *material1 = new Lambertian(color(0.9, 0.9, 0.1), false, color(0, 0, 0));
    scene.add(new Cube(Point3lf(0, 0, -150), Point3lf(300, 200, -140), material1));
    
    // right wall
    Material *material2 = new Lambertian(color(0.1, 0.9, 0.1), false, color(0, 0, 0));
    scene.add(new Cube(Point3lf(0, 0, 140), Point3lf(300, 200, 150), material2));

    Material *light = new Lambertian(color(1, 1, 1), true, color(5, 5, 5));

    Sampler2D *sampler = new RandomSampler2D(100, 100);
    Point3lf p = Point3lf(150, 199, 0);
    float r = 50;
    Vec3lf pNormal = Vec3lf(0, -1, 0);
    Vec3lf px = Vec3lf(1, 0, 0);
    // scene.add(new Disk(p, px, pNormal, r, light, 1, 30, sampler));

    Sampler2D *triSampler = new RandomSampler2D(1, 1);
    Point3lf p0 = Point3lf(200, 199, -40);
    Point3lf p1 = Point3lf(100, 199, 0);
    Point3lf p2 = Point3lf(100, 199, 80);
    Point3lf p3 = Point3lf(200, 199, 40);
    std::vector<Point3lf> p_list;
    p_list.push_back(p0);
    p_list.push_back(p1);
    p_list.push_back(p2);
    p_list.push_back(p3);
    // scene.add(new Polygon(p_list, light, triSampler));

    Point3lf p10 = Point3lf(200, 199, -20);
    Point3lf p11 = Point3lf(200, 199, 20);
    Point3lf p12 = Point3lf(150, 199, 40);
    Point3lf p13 = Point3lf(100, 199, 20);
    Point3lf p14 = Point3lf(100, 199, -20);
    Point3lf p15 = Point3lf(150, 199, -40);
    std::vector<Point3lf> p_list2;
    p_list2.push_back(p10);
    p_list2.push_back(p11);
    p_list2.push_back(p12);
    p_list2.push_back(p13);
    p_list2.push_back(p14);
    p_list2.push_back(p15);
    scene.add(new Polygon(p_list2, light, triSampler));
    /*
    // Cylinder
    scene.add(new Cylinder(Vec3lf(300, 0, -70), 20, 80, ground_material1));
    scene.add(new Disk(Vec3lf(300, 80, -70), Vec3lf(1, 0, 0), Vec3lf(0, 1, 0), 20, ground_material1, 0.1, 10));
    */
    // sphere
    
    scene.add(new Sphere(Vec3lf(170, 30, 80), 25, m1));
    scene.add(new Sphere(Vec3lf(170, 85, 80), 25, m2));
    scene.add(new Sphere(Vec3lf(170, 140, 80), 25, m3));
    scene.add(new Sphere(Vec3lf(170, 30, 0), 25, m4));
    scene.add(new Sphere(Vec3lf(170, 85, 0), 25, m5));
    scene.add(new Sphere(Vec3lf(170, 140, 0), 25, m6));
    scene.add(new Sphere(Vec3lf(170, 30, -80), 25, m7));
    scene.add(new Sphere(Vec3lf(170, 85, -80), 25, m8));
    scene.add(new Sphere(Vec3lf(170, 140, -80), 25, m9));
    

   // scene.add(new Sphere(Vec3lf(170, 30, 0), 25, m4));

    /*
    Material *material3 = new Lambertian(color(0.1, 0.1, 0.9), false, color(0,0,0));
    Matrix mat = ModelTransform(Vec3lf(50, 0, 0), Vec3lf(0, 90, 0), Vec3lf(60, 60, 60));
    scene.add(new MeshTriangle("./obj/bunny.obj", "", material3, mat, false));
    */
    // Material *material4 = new Lambertian(color(0.9, 0.9, 0.1), false, color(0, 0, 0));
    // Matrix mat = ModelTransform(Vec3lf(130, 0, -50), Vec3lf(0, 90, 0), Vec3lf(7.5, 7.5, 7.5));
    // scene.add(new MeshTriangle("./obj/bunny.obj", "./obj/bunny.obj", material4, mat, false));
    // Material *ball = new Specular(1.8, false, true);
    // scene.add(new Sphere(Point3lf(1, 0.4, 0.5), 0.4, ball));

    return scene;
}

Hittable_list scene_3()
{
    Hittable_list scene;

    Material *lightr = new Lambertian(color(1, 1, 1), true, color(5, 0, 0));
    Material *lightg = new Lambertian(color(1, 1, 1), true, color(0, 5, 0));
    Material *lightb = new Lambertian(color(1, 1, 1), true, color(0, 0, 5));

    Sampler2D *sampler = new RandomSampler2D(100, 100);
    Point3lf p1 = Point3lf(200, 150, 0);
    float r1 = 20;
    Vec3lf pNormal1 = Vec3lf(0, -1, 0);
    Vec3lf px1 = Vec3lf(1, 0, 0);
    scene.add(new Disk(p1, px1, pNormal1, r1, lightr, 1, 30, sampler));

    Point3lf p2 = Point3lf(200, 150, -80);
    float r2 = 5;
    Vec3lf pNormal2 = Vec3lf(0, -1, 0);
    Vec3lf px2 = Vec3lf(1, 0, 0);
    scene.add(new Disk(p2, px2, pNormal2, r2, lightg, 1, 30, sampler));

    Point3lf p3 = Point3lf(200, 150, 80);
    float r3 = 40;
    Vec3lf pNormal3 = Vec3lf(0, -1, 0);
    Vec3lf px3 = Vec3lf(1, 0, 0);
    scene.add(new Disk(p3, px3, pNormal3, r3, lightb, 1, 30, sampler));

    // floor
    Material *m1 = new CookTorrance(Vec3f(0.91f, 0.92f, 0.92f), Vec3f(0.96f, 0.96f, 0.97f) , 1, 0.2, false, color(0,0,0));
    scene.add(new Cube(Point3lf(-400, -1, -400), Point3lf(400, 50, 400), m1));

    return scene;

}


/////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
#ifdef BRIDSON_OUTPUT
    {
        BridsonSampler2D *sampler = new BridsonSampler2D(100, 100, 1);
        for (auto v : sampler->results)
        {
            std::cout << v.x << " " << v.y << std::endl;
        }
    }
#else
    {
        clock_t startTime, endTime;
        startTime = clock();
        // scene
        auto scene = scene_2();
        scene.BuildBVH();
        // Camera
        Point3lf lookfrom(-400, 110, 0);
        Point3lf lookat(0, 90, 0);
        Vec3lf vup = cross(lookat - lookfrom, Vec3lf(0, 0, -1));

        Camera cam(lookfrom, lookat, vup, 30, aspect_ratio);

        Renderer r;

        r.Render(scene, cam);

        endTime = clock();
        std::cerr << "The run time is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
    }
#endif
}