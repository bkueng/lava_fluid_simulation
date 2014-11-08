/*! \file FastMath.h
    \brief Fast math routines.
    \author Wojciech Jarosz
*/
#ifndef MATH_FAST_MATH_H
#define MATH_FAST_MATH_H

#include <math.h>
#include "MathT.h"

namespace Math
{

//! Interpret the bits of a float as an int
inline int
floatAsInt(float f)
{
    union {int i; float f;} u;
    u.f = f;
    return u.i;
}

//! Interpret the bits of an int as a float
inline float
intAsFloat(int i)
{
    union {int i; float f;} u;
    u.i = i;
    return u.f;
}

inline int
iLog2(float f)
{
    return ((floatAsInt(f) & 0x7f800000) >> 23) - 0x7f;
}

//! Fast base-2 logarithm of an integer
template <typename T>
inline int
iLog2(T f)
{
    return iLog2(float(f));
}

//! Fast base-2 logarithm of a float
inline float
fastLog2(float f)
{
    int i = floatAsInt(f);
    return (((i&0x7f800000)>>23)-0x7f)+(i&0x007fffff)/(float)0x800000;
}

//! Another fast base-2 logarithm of a float
inline float
betterLog2(float a)
{
    float x = (float) floatAsInt(a);
    x *= 1.1920928955078125e-07f;
    x -= 127.0f;
    
    float y = x-floorf(x);
    y = (y-y*y)*0.346607f;
    return x+y;
}


//! Fast approximation to exponential function
inline double
fastExp(double y)
{
    union
    {
        double d;
// #ifdef LITTLE_ENDIAN
        struct { int j, i; } n;
// #elseif
//         struct { int i, j; } n;
// #endif
    } eco;

    eco.n.i = (int)((1048576/M_LN2)*(y)) + (1072693248 - 60801);
    eco.n.j = 0;
    return eco.d;
}

inline float
fastPow2(float i)
{
    float y = i - floorf(i);
    y = (y - y * y) * 0.33971f;

    float x = i + 127.0f - y;
    x *= 8388608.0f;
    
    return intAsFloat(int(x));
}

//! Fast floating-point pow function
inline float
fastPow(float a, float b)
{
    return fastPow2(b*fastLog2(a));
}

//! Another fast floating-point pow function
inline float
betterPow(float a, float b)
{
    return fastPow2(b*betterLog2(a));
}

//! Fast inverse square-root function
inline float
fastInvSqrt(float x)
{
    if (fabsf(x) == 0.0f)
        return x;

    int i = floatAsInt(x);    // get bits for floating value
    i = 0x5f375a86 - (i>>1);  // gives initial guess y0
    float r = intAsFloat(i);  // convert bits back to float
    
    // Newton iteration, repeated for more accuracy
    float halfX = 0.5f*x;
    r = r * (1.5f - halfX * r * r);
    r = r * (1.5f - halfX * r * r);
    r = r * (1.5f - halfX * r * r);
    r = r * (1.5f - halfX * r * r);
    return r;
}

//! Fast square-root function
inline float
fastSqrt(float x)
{
    return x * fastInvSqrt(x);
}

//! Test if an integer is a power of 2
inline bool
isPowerOf2(int v)
{
    return (v & (v - 1)) == 0;
}

//! Round up to the next smallest power of two
inline unsigned
roundUpPow2(unsigned v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return v+1;
}

//! Round down to the previous largest power of two
inline unsigned
roundDownPow2(unsigned v)
{
	v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return v - (v >> 1);
	
    //if (v == 0)
//        return 0;
//
//    if (v == 1)
//        return 1;
//
//    return 1 << (iLog2(v));
}


#define FAST_INT 1

#ifdef FAST_INT
inline int
round2Int(double val)
{
    // this version from PBRT doesn't work for me, have to use a union instead
//     val += 6755399441055744.0;
//     return ((long*)&val)[0];
    union {long l[2]; double d;} v;
    v.d = val + 6755399441055744.0;
    return int(v.l[0]);
}

inline int
float2Int(double val)
{
    return (val < 0.0) ? round2Int(val+(0.5-1.5e-8)) :
		                 round2Int(val-(0.5-1.5e-8));
}

inline int floor2Int(double val) {return round2Int(val-(0.5-1.5e-8));}
inline int ceil2Int(double val)  {return round2Int(val+(0.5-1.5e-8));}

#else
inline int round2Int(double val) {return int(floor(val+0.5));}
inline int float2Int(double val) {return (int)val;}
inline int floor2Int(double val) {return (int)floor(val);}
inline int ceil2Int(double val)  {return (int)ceil(val);}
#endif

} // namespace Math

#endif // MATH_MATH_H
