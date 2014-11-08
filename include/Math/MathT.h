/*! \file MathT.h
    \brief Templated versions of standard math routines
    \author Wojciech Jarosz
        
	This code is based almost completely on the ImathMath.h file from the 
	OpenEXR project. Here is their copyright notice:
 
	\copyright
    ----------------------------------------------------------------
    
    Copyright(c) 2002, Industrial Light & Magic, a division of Lucas
    Digital Ltd. LLC
    
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:
    *       Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
    *       Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following disclaimer
    in the documentation and/or other materials provided with the
    distribution.
    *       Neither the name of Industrial Light & Magic nor the names of
    its contributors may be used to endorse or promote products derived
    from this software without specific prior written permission. 
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef MATH_MATHT_H
#define MATH_MATHT_H

#include <math.h>
#include <Math/LimitsT.h>

#ifndef M_E
# define M_E            2.7182818284590452354   /* e */
#endif
#ifndef M_LOG2E
# define M_LOG2E        1.4426950408889634074   /* log_2 e */
#endif
#ifndef M_LOG10E
# define M_LOG10E       0.43429448190325182765  /* log_10 e */
#endif
#ifndef M_LN2
# define M_LN2          0.69314718055994530942  /* log_e 2 */
#endif
#ifndef M_LN10
# define M_LN10         2.30258509299404568402  /* log_e 10 */
#endif
#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif
#ifndef M_PI_2
# define M_PI_2         1.57079632679489661923  /* pi/2 */
#endif
#ifndef M_PI_4
# define M_PI_4         0.78539816339744830962  /* pi/4 */
#endif
#ifndef M_1_PI
# define M_1_PI         0.31830988618379067154f  /* 1/pi */
#endif
#ifndef M_2_PI
# define M_2_PI         0.63661977236758134308  /* 2/pi */
#endif
#ifndef M_2_SQRTPI
# define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
#endif
#ifndef M_SQRT2
# define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
#endif
#ifndef M_SQRT1_2
# define M_SQRT1_2      0.70710678118654752440  /* 1/sqrt(2) */
#endif

#ifdef _WIN32
inline double exp2(double x) {return exp(x)/exp(2.0);}
inline double expm1(double x){return exp(x-1);}
inline double log1p(double x){return log(x+1);}
inline double log2(double x) {return log(x)/log(2.0);}
inline double cbrt(double x) {return pow(x, 1.0/3.0);}

inline float exp2f(float x) {return float(exp2(x));}
inline float expm1f(float x){return float(expm1f(x));}
inline float log1pf(float x){return float(log1pf(x));}
inline float log2f(float x) {return float(log2(x));}
inline float cbrtf(float x) {return float(cbrt(x));}
//inline float hypotf(float x, float y) {return _hypotf(x, y);}
#endif

namespace Math
{

//! Templated versions of standard math routines
/*!
    The Math class provides the ability to choose the appropriate precision
    math routine based on a template parameter.
*/
template <typename T>
struct MathT
{
    static T acos (T x)      {return ::acos(double(x));}
    static T asin (T x)      {return ::asin(double(x));}
    static T atan (T x)      {return ::atan(double(x));}
    static T atan2(T x, T y) {return ::atan2(double(x), double(y));}
    static T cos  (T x)      {return ::cos(double(x));}
    static T sin  (T x)      {return ::sin(double(x));}
    static T tan  (T x)      {return ::tan(double(x));}
    static T exp  (T x)      {return ::exp(double(x));}
    static T exp2 (T x)      {return ::exp2(double(x));}
    static T expm1(T x)      {return ::expm1(double(x));}
    static T log  (T x)      {return ::log(double(x));}
    static T log10(T x)      {return ::log10(double(x));}
    static T log1p(T x)      {return ::log1p(double(x));}
    static T log2 (T x)      {return ::log2(double(x));}
    static T modf (T x, T *iptr)
    {
        double ival;
        T rval(::modf(double(x), &ival));
        *iptr = ival;
        return rval;
    }
    static T pow  (T x, T y) {return ::pow(double(x), double(y));}
    static T cbrt (T x)      {return ::cbrt(double(x));}
    static T sqrt (T x)      {return ::sqrt(double(x));}
    static T sqrt2(T x)      {return x > 0 ? ::sqrt((double)x) : 0;}
    static T sqr(T x)        {return x*x;}
    static T fabs (T x)      {return ::fabs(double(x));}
    static T ceil (T x)      {return ::ceil(double(x));}
    static T floor(T x)      {return ::floor(double(x));}
    static T fmod (T x, T y) {return ::fmod(double(x), double(y));}
    static T hypot(T x, T y) {return ::hypot(double(x), double(y));}
};


//! Specialization for single-precision math routines.
template <>
struct MathT<float>
{
    static float acos (float x)           {return ::acosf(x);}	
    static float asin (float x)           {return ::asinf(x);}
    static float atan (float x)           {return ::atanf(x);}
    static float atan2(float x, float y)  {return ::atan2f(x, y);}
    static float cos  (float x)           {return ::cosf(x);}
    static float sin  (float x)           {return ::sinf(x);}
    static float tan  (float x)           {return ::tanf(x);}
    static float exp  (float x)           {return ::expf(x);}
    static float exp2 (float x)           {return ::exp2f(x);}
    static float expm1(float x)           {return ::expm1f(x);}
    static float log  (float x)           {return ::logf(x);}
    static float log10(float x)           {return ::log10f(x);}
    static float log1p(float x)           {return ::log1pf(x);}
    static float log2 (float x)           {return ::log2f(x);}
    static float modf (float x, float *y) {return ::modff(x, y);}
    static float pow  (float x, float y)  {return ::powf(x, y);}
    static float sqrt (float x)           {return ::sqrtf(x);}
    static float sqrt2(float x)           {return x > 0 ? ::sqrtf(x) : 0;}
    static float cbrt (float x)           {return ::cbrtf(x);}
    static float ceil (float x)           {return ::ceilf(x);}
    static float fabs (float x)           {return ::fabsf(x);}
    static float floor(float x)           {return ::floorf(x);}
    static float fmod (float x, float y)  {return ::fmodf(x, y);}
    //static float hypot(float x, float y)  {return ::hypotf(x, y);}
};


/*!
    \returns hypotenuse of real(non-complex) scalars a and b by 
             avoiding underflow/overflow
             using(a * sqrt( 1 + (b/a) * (b/a))), rather than
             sqrt(a*a + b*b).
*/
template <typename T>
T hypot(const T &a, const T &b)
{
    if(a == 0)
        return MathT<T>::fabs(b);
    else
    {
        T c = b/a;
        return MathT<T>::fabs(a) * MathT<T>::sqrt(1 + c*c);
    }
}


/*!
    The sinc function:
        sinc(x) = sin(x)/x
*/
template <typename T>
inline T
sinc(T x)
{
    if (x * x < Limits<T>::epsilon())
        return T(1);
    else
        return MathT<T>::sin(x) / x;
}


//! Compare two numbers and test if they are "approximately equal".
/*!
    \returns true iff x1 is the same as x2 with an absolute error of
             no more than e,

                abs(x1 - x2) <= e
*/
template <typename T>
inline bool
equalWithAbsError(T x1, T x2, T e)
{
    return ((x1 > x2) ? x1 - x2 : x2 - x1) <= e;
}


//! Compare two numbers and test if they are "approximately equal".
/*!
    \returns true iff x1 is the same as x2 with an relative error of
             no more than e,

                abs(x1 - x2) <= e * x1
*/
template <typename T>
inline bool
equalWithRelError(T x1, T x2, T e)
{
    return ((x1 > x2) ? x1 - x2: x2 - x1) <= e * ((x1 > 0) ? x1: -x1);
}


//struct FactHelper
//{
//	//! Size of the precomputed factorial array.
//	static const int FACT_ARRAY_SIZE = 35;
//
//	//! Precomputed factorial array.
//	static const double s_factArray[FACT_ARRAY_SIZE];
//};


//! Return n-factorial as a double (table-based).
inline double factorial(int n)
{
	static const int FACT_ARRAY_SIZE = 35;
	static const double factArray[FACT_ARRAY_SIZE] = 
	{1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3.6288e+006,
	 3.99168e+007, 4.79002e+008, 6.22702e+009, 8.71783e+010,
	 1.30767e+012, 2.09228e+013, 3.55687e+014, 6.40237e+015,
	 1.21645e+017, 2.4329e+018, 5.10909e+019, 1.124e+021,
	 2.5852e+022, 6.20448e+023, 1.55112e+025, 4.03291e+026,
	 1.08889e+028, 3.04888e+029, 8.84176e+030, 2.65253e+032,
	 8.22284e+033, 2.63131e+035, 8.68332e+036, 2.95233e+038
	};

//	assert(n >= 0);

	// take the value from a table
	if (n < FACT_ARRAY_SIZE)
		return factArray[n];
		
	// compute the value traditionally starting from the
	// maximum value in the array
	double fact = 1;
	while (n >= FACT_ARRAY_SIZE)
	  fact *= n--; /* Multiply, then decrement. */
	fact *= factArray[FACT_ARRAY_SIZE-1];

	return fact;
}

inline float  radians(float deg)  {return((float)M_PI/180.f)*deg;}
inline double radians(double deg) {return M_PI/180.0*deg;}
inline float  degrees(float rad)  {return(180.f/(float)M_PI)*rad;}
inline double degrees(double rad) {return 180.0/M_PI*rad;}

} // namespace Math

#endif // MATH_MATH_H
