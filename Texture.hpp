#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#include <string>

class Texture{
public:
    Texture(){};
public:
    bool load_xy(float x, float y, unsigned char& r, unsigned char& g,unsigned char& b);
public:
    unsigned char* data;
    int w;
    int h;
    int n;
};

#endif