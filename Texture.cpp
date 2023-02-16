#include "Texture.hpp"
#include <cmath>

bool Texture::load_xy(float x, float y, unsigned char &r, unsigned char &g, unsigned char &b)
{
    x *= w;
    y *= h;
    int ix = x - floor(x) < 0.5 ? floor(x) : ceil(x);
    int iy = y - floor(y) < 0.5 ? floor(y) : ceil(y);

    r = data[(ix * w + iy) * 3 + 0];
    g = data[(ix * w + iy) * 3 + 1];
    b = data[(ix * w + iy) * 3 + 2];
    return true;
}