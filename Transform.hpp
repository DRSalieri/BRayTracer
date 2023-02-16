#ifndef __TRANSFORM_H__
#define __TRANSFORM_H__

#include "geometry.hpp"

Matrix ModelTransform(
    Vec3lf translation = Vec3lf(0, 0, 0),
    Vec3lf rotation = Vec3lf(0, 0, 0) ,
    Vec3lf scale = Vec3lf(1,1,1));

#endif