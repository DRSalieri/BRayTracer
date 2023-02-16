#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include "geometry.hpp"
#include <iostream>

// samples_per_pixel
void write_color(std::ostream &out, color pixel_color, int samples_per_pixel) {
    auto r = pixel_color.x;
    auto g = pixel_color.y;
    auto b = pixel_color.z;

    // Divide the color by the number of samples and gamma-correct for gamma=2.0.
    auto scale = 1.0 / samples_per_pixel;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    /*
    int abc =  static_cast<int>(255.99 * clamp(r, 0.0, 0.999)) ;
    if(abc < 0)
    {
        printf("1");
    }
    */
   
    // Write the translated [0,255] value of each color component.
    out << static_cast<int>(255.99 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(255.99 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(255.99 * clamp(b, 0.0, 0.999)) << '\n';
}

#endif // !COLOR_H