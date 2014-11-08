/*! \file Core.h
    \brief Contains the definition of various utility functions.
    \author Wojciech Jarosz
*/
#ifndef MATH_CORE_H
#define MATH_CORE_H

#include <algorithm>

namespace Math
{

const float FEQ_EPS  = 1e-5f;
const float FEQ_EPS2 = 1e-12f;
const float FEQ_INF  = 1e12f;


template <typename T>
inline T
abs(T a)
{
    return a < 0 ? -a : a;
}


template <typename T>
inline T sign(T a) {return (a > 0) ? T (1) : (a < 0) ? T (-1) : 0;}


template <typename T>
inline bool
isZero(T a, T t)
{
    return Math::abs(a) <= t;
}


template <typename T1, typename T2, typename T3>
inline bool
equal(T1 a, T2 b, T3 t)
{
    return Math::abs(a - b) <= t;
}


//! Linear interpolation.
/*!
    Linearly interpolates between \a a and \a b, using parameter t.
    \param a A value.
    \param b Another value.
    \param t A blending factor of \a a and \a b.
    \return Linear interpolation of \a a and \b -
            a value between a and b if \a t is between 0 and 1.
*/
template <typename T, typename S>
inline T
lerp(T a, T b, S t)
{
    return T((S(1)-t) * a + t * b);
}


//! Inverse linear interpolation.
template <typename T>
inline T
lerpFactor(T a, T b, T m)
{
    return (m - a) / (b - a);
}


//! Barycentric interpolation.
template <typename T, typename S>
inline T
bary(T a, T b, T c, S s, S t)
{
    return (S(1)-s-t)*a + s*b + t*c;
}


//! Clamps a double between two bounds.
/*!
    \param a The value to test.
    \param l The lower bound.
    \param h The upper bound.
    \return The value \a a clamped to the lower and upper bounds.
    
    This function has been specially crafted to prevent NaNs from propagating.
*/
template <typename T>
inline T
clamp(T a, T l, T h)
{
    return (a >= l) ? ((a <= h) ? a : h) : l;
}


//! Returns a modulus b.
template <typename T>
inline T
mod(T a, T b)
{
    int n = (int)(a/b);
    a -= n*b;
    if (a < 0)
        a += b;
    return a;
}


//! Returns a modulus 1, assuming a is positive.
template <typename T>
inline T
mod1(T a)
{
    return a - int(a);
}


template <typename T> inline T pow2 (T x) {return x*x;}
template <typename T> inline T pow3 (T x) {return x*x*x;}
template <typename T> inline T pow4 (T x) {T x2 = x*x; return x2*x2;}
template <typename T> inline T pow5 (T x) {T x2 = x*x; return x2*x2*x;}
template <typename T> inline T sqr  (T x) {return pow2(x);}
template <typename T> inline T cube (T x) {return pow3(x);}

} // namespace Math


namespace std
{

template <typename T>
inline const T&
min(const T& a, const T& b, const T& c)
{
    return std::min(std::min(a, b), c);
}

template <typename T>
inline const T&
min(const T& a, const T& b, const T& c, const T& d)
{
    return std::min(std::min(a, b, c), d);
}

template <typename T>
inline const T&
min(const T& a, const T& b, const T& c, const T& d, const T& e)
{
    return std::min(std::min(a, b, c, d), e);
}

template <typename T>
inline const T&
max(const T& a, const T& b, const T& c)
{
    return std::max(std::max(a, b), c);
}

template <typename T>
inline const T&
max(const T& a, const T& b, const T& c, const T& d)
{
    return std::max(std::max(a, b, c), d);
}

template <typename T>
inline const T&
max(const T& a, const T& b, const T& c, const T& d, const T& e)
{
    return std::max(std::max(a, b, c, d), e);
}

} // namespace std

#endif  // MATH_CORE_H
