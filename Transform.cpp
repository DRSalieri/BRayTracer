#include "Transform.hpp"
#include "geometry.hpp"

// rotation with degree(not rad)
Matrix ModelTransform(Vec3lf translation,
                      Vec3lf rotation,
                      Vec3lf scale)
{
    Matrix translationM = Matrix::identity();
    translationM[0][3] = translation[0];
    translationM[1][3] = translation[1];
    translationM[2][3] = translation[2];

    Matrix rotationMx = Matrix::identity();
    Matrix rotationMy = Matrix::identity();
    Matrix rotationMz = Matrix::identity();
    float rad, sin, cos;
    if (rotation[0] != 0)
    {
        rad = deg2rad(rotation[0]);
        cos = cosf(rad);
        sin = sinf(rad);
        rotationMx[1][1] = cos;
        rotationMx[1][2] = -sin;
        rotationMx[2][1] = sin;
        rotationMx[2][2] = cos;
    }
    if (rotation[1] != 0)
    {
        rad = deg2rad(rotation[1]);
        cos = cosf(rad);
        sin = sinf(rad);
        rotationMx[0][0] = cos;
        rotationMx[0][2] = sin;
        rotationMx[2][0] = -sin;
        rotationMx[2][2] = cos;
    }
    if (rotation[2] != 0)
    {
        rad = deg2rad(rotation[2]);
        cos = cosf(rad);
        sin = sinf(rad);
        rotationMx[0][0] = cos;
        rotationMx[0][1] = -sin;
        rotationMx[1][0] = sin;
        rotationMx[1][1] = cos;
    }

    Matrix scaleM = Matrix::identity();
    scaleM[0][0] = scale[0];
    scaleM[1][1] = scale[1];
    scaleM[2][2] = scale[2];

    return translationM * rotationMx * rotationMy * rotationMz * scaleM;
}