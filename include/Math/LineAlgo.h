/*! \file LineAlgo.h
    \brief Contains extra Line algorithms.
    \author Wojciech Jarosz
*/
#ifndef MATH_LINE_ALGO_H
#define MATH_LINE_ALGO_H

#include <Math/Line.h>

namespace Math
{

//Templated Nd Box-Ray Intersection of a ray (ro,rd) and a box, within tMin and tMax
//Returns distances to the intersections hitt0 and hitt1
template <typename Vec>
inline bool
intersects(const Vec & ro, const Vec & rd, const Box<Vec>& box,
           typename Vec::BaseType tMin = Limits<typename Vec::BaseType>::min(),
           typename Vec::BaseType tMax = Limits<typename Vec::BaseType>::max(),
           typename Vec::BaseType * hitt0 = 0,
           typename Vec::BaseType * hitt1 = 0)
{
    typedef typename Vec::BaseType T;
    for (unsigned i = 0; i < Vec::dimensions(); ++i)
    {
        // Update interval for ith bounding box slab
        T invRayDir = T(1) / rd[i];
        T tNear = (box.min[i] - ro[i]) * invRayDir;
        T tFar  = (box.max[i] - ro[i]) * invRayDir;
    
        // Update parametric interval from slab intersection ts
        if (tNear > tFar)
            std::swap(tNear, tFar);
        tMin = tNear > tMin ? tNear : tMin;
        tMax = tFar  < tMax ? tFar  : tMax;
        if (tMin > tMax)
            return false;
    }

    if (hitt0)
        *hitt0 = tMin;
    if (hitt1)
        *hitt1 = tMax;
    return true;
}

//Templated Intersection test between a ray (ro,rd) and a triangle (A,B,C)
//in the interval t0,t1. Returns distance t, uv coordinates and the geometric normal Ng
template <typename T>
inline bool
intersect(const Vec3<T> & ro, const Vec3<T> & rd,
          const Vec3<T>& A,
          const Vec3<T>& B,
          const Vec3<T>& C,
          T t0, T t1,
          T * t, Vec2<T> * uv, Vec3<T> * Ng)
{
    Vec3<T> AB = B - A;
    Vec3<T> AC = C - A;
    Vec3<T> N = cross(AB, AC);
    T det = dot(rd, N);
    
    if (det == T(0))
       return false;
        
    T detInv = T(1) / det;
    Vec3<T> toA = ro - A;

    T a = dot(rd, cross(toA, AC)) * detInv;
    if (a < T(0))
        return false;

    T b = dot(rd, cross(AB, toA)) * detInv;
    if (b < T(0) || a + b > T(1))
        return false;

    T tempT = -dot(toA, N) * detInv;
    if (tempT < t0 || tempT > t1)
        return false;

    *t = tempT;
    uv->set(a, b);
    *Ng = N;
    return true;
}


//Templated (fast) Intersection test between a ray (ro,rd) and a triangle (A,B,C)
//in the interval t0,t1. Returns only true or false.
template <typename T>
inline bool
intersects(const Vec3<T> & ro, const Vec3<T> & rd,
           const Vec3<T>& A,
           const Vec3<T>& B,
           const Vec3<T>& C,
           T t0, T t1)
{
    Vec3<T> AB = B - A;
    Vec3<T> AC = C - A;
    Vec3<T> N = cross(AB, AC);
    float det = dot(rd, N);

    if (det == T(0))
       return false;

    T detInv = 1.0f / det;
    Vec3<T> toA = ro - A;

    T a = dot(rd, cross(toA, AC)) * detInv;
    if (a < T(0))
        return false;

    T b = dot(rd, cross(AB, toA)) * detInv;
    if (b < T(0) || a + b > T(1))
        return false;

    T t = -dot(toA, N) * detInv;
    if (t < t0 || t > t1)
        return false;

    return true;
}

} // namespace Math

#endif  // MATH_LINE_ALGO_H
