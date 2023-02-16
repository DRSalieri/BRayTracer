
#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <cmath>
#include <limits>
#include <vector>
#include <cassert>
#include <iostream>
#include <tuple>
template <size_t DimCols, size_t DimRows, typename T>
class mat;

template <size_t DIM, typename T>
struct vec
{
    vec()
    {
        for (size_t i = DIM; i--; data_[i] = T())
            ;
    }
    T &operator[](const size_t i)
    {
        assert(i < DIM);
        return data_[i];
    }
    const T &operator[](const size_t i) const
    {
        assert(i < DIM);
        return data_[i];
    }

private:
    T data_[DIM];
};

/////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct vec<2, T>
{
    vec() : x(T()), y(T()) {}
    vec(T X, T Y) : x(X), y(Y) {}
    template <class U>
    vec<2, T>(const vec<2, U> &v);
    T &operator[](const size_t i)
    {
        assert(i < 2);
        return i <= 0 ? x : y;
    }
    const T &operator[](const size_t i) const
    {
        assert(i < 2);
        return i <= 0 ? x : y;
    }
    float norm() const { return std::sqrt(x * x + y * y); }
    float norm2() const { return x * x + y * y; }

    T x, y;
};

/////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct vec<3, T>
{
    vec() : x(T()), y(T()), z(T()) {}
    vec(T X, T Y, T Z) : x(X), y(Y), z(Z) {}
    template <class U>
    vec<3, T>(const vec<3, U> &v);
    T &operator[](const size_t i)
    {
        assert(i < 3);
        return i <= 0 ? x : (1 == i ? y : z);
    }
    const T &operator[](const size_t i) const
    {
        assert(i < 3);
        return i <= 0 ? x : (1 == i ? y : z);
    }
    float norm() const { return std::sqrt(x * x + y * y + z * z); }
    float norm2() const { return x * x + y * y + z * z; }
    vec<3, T> normalize(T l = 1) const
    {
        return (*this) * (l / norm());
    }
    vec<3, T> &operator+=(const vec<3, T> &v)
    {
        x += v[0];
        y += v[1];
        z += v[2];
        return *this;
    }
    vec<3, T> &operator-=(const vec<3, T> &v)
    {
        x -= v[0];
        y -= v[1];
        z -= v[2];
        return *this;
    }
    vec<3, T> &operator*=(const vec<3, T> &v)
    {
        x *= v[0];
        y *= v[1];
        z *= v[2];
        return *this;
    }
    static vec<3, T> Min(const vec<3, T> &v1, const vec<3, T> &v2)
    {
        vec<3, T> v;
        v.x = std::min(v1.x, v2.x);
        v.y = std::min(v1.y, v2.y);
        v.z = std::min(v1.z, v2.z);
        return v;
    }
    static vec<3, T> Max(const vec<3, T> &v1, const vec<3, T> &v2)
    {
        vec<3, T> v;
        v.x = std::max(v1.x, v2.x);
        v.y = std::max(v1.y, v2.y);
        v.z = std::max(v1.z, v2.z);
        return v;
    }
    T x, y, z;
};

/////////////////////////////////////////////////////////////////////////////////

template <size_t DIM, typename T>
T operator*(const vec<DIM, T> &lhs, const vec<DIM, T> &rhs)
{
    T ret = T();
    for (size_t i = DIM; i--; ret += lhs[i] * rhs[i])
        ;
    return ret;
}

template <size_t DIM, typename T>
vec<DIM, T> dotMul(vec<DIM, T> lhs, const vec<DIM, T> &rhs)
{
    for (size_t i = DIM; i--; lhs[i] *= rhs[i])
        ;
    return lhs;
}

template <size_t DIM, typename T>
vec<DIM, T> operator+(vec<DIM, T> lhs, const vec<DIM, T> &rhs)
{
    for (size_t i = DIM; i--; lhs[i] += rhs[i])
        ;
    return lhs;
}

template <size_t DIM, typename T>
vec<DIM, T> operator-(vec<DIM, T> lhs, const vec<DIM, T> &rhs)
{
    for (size_t i = DIM; i--; lhs[i] -= rhs[i])
        ;
    return lhs;
}

template <size_t DIM, typename T>
vec<DIM, T> operator-(vec<DIM, T> v)
{
    for (size_t i = DIM; i--; v[i] = -v[i])
        ;
    return v;
}

template <size_t DIM, typename T>
T length_squared(const vec<DIM, T> &v)
{
    T ret = 0;
    for (size_t i = DIM; i--; ret += (v[i] * v[i]))
        ;
    return ret;
}

template <size_t DIM, typename T, typename U>
vec<DIM, T> operator*(vec<DIM, T> lhs, const U &rhs)
{
    for (size_t i = DIM; i--; lhs[i] *= rhs)
        ;
    return lhs;
}

template <size_t DIM, typename T, typename U>
vec<DIM, T> operator*(const U &lhs, vec<DIM, T> rhs)
{
    for (size_t i = DIM; i--; rhs[i] *= lhs)
        ;
    return rhs;
}

template <size_t DIM, typename T, typename U>
vec<DIM, T> operator/(vec<DIM, T> lhs, const U &rhs)
{
    for (size_t i = DIM; i--; lhs[i] /= rhs)
        ;
    return lhs;
}

template <size_t LEN, size_t DIM, typename T>
vec<LEN, T> embed(const vec<DIM, T> &v, T fill = 1)
{
    vec<LEN, T> ret;
    for (size_t i = LEN; i--; ret[i] = (i < DIM ? v[i] : fill))
        ;
    return ret;
}

template <size_t LEN, size_t DIM, typename T>
vec<LEN, T> proj(const vec<DIM, T> &v)
{
    vec<LEN, T> ret;
    for (size_t i = LEN; i--; ret[i] = v[i])
        ;
    return ret;
}

template <typename T>
vec<3, T> cross(vec<3, T> v1, vec<3, T> v2)
{
    return vec<3, T>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

template <typename T>
T dot(vec<3, T> v1, vec<3, T> v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <size_t DIM, typename T>
std::ostream &operator<<(std::ostream &out, vec<DIM, T> &v)
{
    for (unsigned int i = 0; i < DIM; i++)
    {
        out << v[i] << " ";
    }
    return out;
}

/////////////////////////////////////////////////////////////////////////////////

template <size_t DIM, typename T>
struct dt
{
    static T det(const mat<DIM, DIM, T> &src)
    {
        T ret = 0;
        for (size_t i = DIM; i--; ret += src[0][i] * src.cofactor(0, i))
            ;
        return ret;
    }
};

template <typename T>
struct dt<1, T>
{
    static T det(const mat<1, 1, T> &src)
    {
        return src[0][0];
    }
};

/////////////////////////////////////////////////////////////////////////////////

template <size_t DimRows, size_t DimCols, typename T>
class mat
{
    vec<DimCols, T> rows[DimRows];

public:
    mat() {}

    vec<DimCols, T> &operator[](const size_t idx)
    {
        assert(idx < DimRows);
        return rows[idx];
    }

    const vec<DimCols, T> &operator[](const size_t idx) const
    {
        assert(idx < DimRows);
        return rows[idx];
    }

    vec<DimRows, T> col(const size_t idx) const
    {
        assert(idx < DimCols);
        vec<DimRows, T> ret;
        for (size_t i = DimRows; i--; ret[i] = rows[i][idx])
            ;
        return ret;
    }

    void set_col(size_t idx, vec<DimRows, T> v)
    {
        assert(idx < DimCols);
        for (size_t i = DimRows; i--; rows[i][idx] = v[i])
            ;
    }

    static mat<DimRows, DimCols, T> identity()
    {
        mat<DimRows, DimCols, T> ret;
        for (size_t i = DimRows; i--;)
            for (size_t j = DimCols; j--; ret[i][j] = (i == j))
                ;
        return ret;
    }

    T det() const
    {
        return dt<DimCols, T>::det(*this);
    }

    mat<DimRows - 1, DimCols - 1, T> get_minor(size_t row, size_t col) const
    {
        mat<DimRows - 1, DimCols - 1, T> ret;
        for (size_t i = DimRows - 1; i--;)
            for (size_t j = DimCols - 1; j--; ret[i][j] = rows[i < row ? i : i + 1][j < col ? j : j + 1])
                ;
        return ret;
    }

    T cofactor(size_t row, size_t col) const
    {
        return get_minor(row, col).det() * ((row + col) % 2 ? -1 : 1);
    }

    mat<DimRows, DimCols, T> adjugate() const
    {
        mat<DimRows, DimCols, T> ret;
        for (size_t i = DimRows; i--;)
            for (size_t j = DimCols; j--; ret[i][j] = cofactor(i, j))
                ;
        return ret;
    }

    mat<DimRows, DimCols, T> invert_transpose()
    {
        mat<DimRows, DimCols, T> ret = adjugate();
        T tmp = ret[0] * rows[0];
        return ret / tmp;
    }

    mat<DimCols, DimRows, T> invert()
    {
        return invert_transpose().transpose();
    }

    mat<DimCols, DimRows, T> transpose()
    {
        mat<DimCols, DimRows, T> ret;
        for (size_t i = DimRows; i--; ret[i] = this->col(i))
            ;
        return ret;
    }
};

/////////////////////////////////////////////////////////////////////////////////

template <size_t DimRows, size_t DimCols, typename T>
vec<DimRows, T> operator*(const mat<DimRows, DimCols, T> &lhs, const vec<DimCols, T> &rhs)
{
    vec<DimRows, T> ret;
    for (size_t i = DimRows; i--; ret[i] = lhs[i] * rhs)
        ;
    return ret;
}

template <size_t R1, size_t C1, size_t C2, typename T>
mat<R1, C2, T> operator*(const mat<R1, C1, T> &lhs, const mat<C1, C2, T> &rhs)
{
    mat<R1, C2, T> result;
    for (size_t i = R1; i--;)
        for (size_t j = C2; j--; result[i][j] = lhs[i] * rhs.col(j))
            ;
    return result;
}

template <size_t DimRows, size_t DimCols, typename T>
mat<DimCols, DimRows, T> operator/(mat<DimRows, DimCols, T> lhs, const T &rhs)
{
    for (size_t i = DimRows; i--; lhs[i] = lhs[i] / rhs)
        ;
    return lhs;
}

template <size_t DimRows, size_t DimCols, class T>
std::ostream &operator<<(std::ostream &out, mat<DimRows, DimCols, T> &m)
{
    for (size_t i = 0; i < DimRows; i++)
        out << m[i] << std::endl;
    return out;
}

/////////////////////////////////////////////////////////////////////////////////

typedef vec<2, float> Vec2f;
typedef vec<2, int> Vec2i;
typedef vec<2, float> Vec2lf;
typedef vec<3, float> Vec3f;
typedef vec<3, int> Vec3i;
typedef vec<3, float> Vec3lf;
typedef vec<4, float> Vec4f;
typedef mat<4, 4, float> Matrix;

/////////////////////////////////////////////////////////////////////////////////

typedef vec<3, float> color;
typedef vec<3, float> Point3lf;

/////////////////////////////////////////////////////////////////////////////////

/*
template <> template <> vec<3,int>  ::vec(const vec<3,float> &v) : x(int(v.x+.5f)),y(int(v.y+.5f)),z(int(v.z+.5f)) {}
template <> template <> vec<3,float>::vec(const vec<3,int> &v)   : x(v.x),y(v.y),z(v.z) {}
template <> template <> vec<2,int>  ::vec(const vec<2,float> &v) : x(int(v.x+.5f)),y(int(v.y+.5f)) {}
template <> template <> vec<2,float>::vec(const vec<2,int> &v)   : x(v.x),y(v.y) {}
*/

/////////////////////////////////////////////////////////////////////////////////

const float infinity = std::numeric_limits<float>::infinity();
const float pi = 3.1415926535897932385;

/////////////////////////////////////////////////////////////////////////////////

inline float degrees_to_radians(float degrees)
{
    return degrees * pi / 180.0;
}

inline float clamp(float x, float min, float max)
{
    if (x < min)
        return min;
    if (x > max)
        return max;
    return x;
}

inline Vec3lf vec_reflect(const Vec3lf &v, const Vec3lf &n)
{
    return v - 2 * dot(v, n) * n;
}

Vec3lf vec_refract(const Vec3lf &I, const Vec3lf &N, const float &ior);

void fresnel(const Vec3lf &I, const Vec3lf &N, const float &ior, float &kr);

/////////////////////////////////////////////////////////////////////////////////

// [min, max)
inline int random_int(int min, int max)
{
    return min + rand() % (max - min);
}

inline float random_float()
{
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

inline float random_float(float min, float max)
{
    // Returns a random real in [min,max).
    return min + (max - min) * random_float();
}
Vec3lf random_Vec3lf();

Vec3lf random_in_unit_sphere();

Vec3lf random_unit_vector();

Vec3lf random_in_unit_disk();

const float EPSILON = 0.001;
const float DIS_EPSILON = 0.01;

inline Vec3lf &EpsilonRay(Vec3lf &orig, const Vec3lf &dir, const Vec3lf &n)
{
    return orig = orig + n * DIS_EPSILON;
}

inline float deg2rad(const float &deg) { return deg * M_PI / 180.0; }

// interpolate
Vec3lf interpolate(float alpha, float gamma, float beta, Vec3lf &vert1, Vec3lf &vert2, Vec3lf &vert3, float weight);

Vec2lf interpolate(float alpha, float gamma, float beta, Vec2lf &vert1, Vec2lf &vert2, Vec2lf &vert3, float weight);

std::tuple<float, float, float> computeBarycentric2D(float x, float y, Vec2lf &vert1, Vec2lf &vert2, Vec2lf &vert3);

Vec3lf toWorld(const Vec3lf& a, const Vec3lf& N);

float lerp(float a, float b, float t);

Vec3f lerp(const Vec3f &a, const Vec3f &b, float t);


#endif //__GEOMETRY_H__