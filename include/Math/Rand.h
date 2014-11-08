/*! \file Rand.h
    \brief Contains the implementation of the Mersenne Twister Pseudo-Random
           number generator.
    \author Wojciech Jarosz
*/
#ifndef RAND_RAND_H
#define RAND_RAND_H

#include <string.h>


#ifdef WIN32

#ifndef snprintf
	#define snprintf _snprintf
#endif 

typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;

#else
	#include <stdint.h>
#endif

namespace Math
{

//! Encapsulated the Mersenne Twister pseudo-random number generator.
/*!
    This random number generator uses the Mersenne Twister method.
    The copyright notice for the Mersenne Twister PRNG is provided below:
    
    A C-program for MT19937, with initialization improved 2002/1/26.
    Coded by Takuji Nishimura and Makoto Matsumoto.
    
    Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
    All rights reserved.                          
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:
    
        1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
    
        2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
    
        3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written 
        permission.
*/
class RandMT
{
public:
    RandMT(uint32_t seed = 0) {init(seed);}
    ~RandMT() {}
    void init(uint32_t seed);


    //! Uniform bool in the range [false, true].
    bool          nextb();
    
    //! Uniform byte in the range [0, 255].
    unsigned char nextc();
    
    //! Uniform float in the range [0.0f, 1.0f).
    float         nextf();
    
    //! Uniform float in the range [\a min, \a max).
    float         nextf(float min, float max);
    
    //! Uniform integer in the range [0, 2^32-1].
    uint32_t      nexti();
    
    //! Uniform integer in the range [\a min, \a max].
    int32_t       nexti(int32_t min, int32_t max);
    
    //! Uniform double in the range [0.0, 1.0).
    double        nextd();
    
    //! Uniform double in the range [\a min, \a max).
    double        nextd(double min, double max);
    
    //! STL RandomNumberGenerator Functor interface.
    int32_t          operator() (int32_t n);
    
protected:
    uint32_t next();

private:
    // Period parameters
    enum constants
    {
        N = 624,
        M = 397
    };
    
    //! Constant vector a.
    static const unsigned long MATRIX_A   = 0x9908b0dfUL;
    //! Most significant w-r bits.
    static const unsigned long UPPER_MASK = 0x80000000UL;
    //! Least significant r bits.
    static const unsigned long LOWER_MASK = 0x7fffffffUL;

    unsigned long mt[N];        //!< The array for the state vector
    int mti;                    //!< mti==N+1 means mt[N] is not initialized
};


inline bool
RandMT::nextb()
{
    return nexti() & 0x80000000UL;
}

inline unsigned char
RandMT::nextc()
{
    return (nexti() & 0xE0000000UL) >> 14;
}

inline float
RandMT::nextf()
{
    return nexti() * (1.0f/4294967296.0f);
}

inline float
RandMT::nextf(float min, float max)
{
    return min + nextf() * (max-min);
}

inline uint32_t
RandMT::nexti()
{
    return next();
}

inline int32_t
RandMT::nexti(int32_t min, int32_t max)
{
    return min + long(nextf() * (max-min) + 0.5f);
}

inline double
RandMT::nextd()
{
    uint32_t a = nexti() >> 5, b = nexti() >> 6; 
    return (a*67108864.0+b)* (1.0/9007199254740992.0);
    //        2^26               2^53    
}

inline double
RandMT::nextd (double min, double max)
{
    return min + nextd() * (max-min);
}

inline int32_t
RandMT::operator()(int32_t max)
{
    return int32_t(nextf() * (max-1) + 0.5f);
}


//! Re-seed
inline void
RandMT::init(uint32_t seed)
{
    // Initialize generator state with seed
    // See Knuth TAOCP Vol 2, 3rd Ed, p.106 for multiplier.
    // In previous versions, most significant bits (MSBs) of the seed affect
    // only MSBs of the state array.  Modified 9 Jan 2002 by Makoto Matsumoto.
    register unsigned long *s = mt;
    register unsigned long *r = mt;
    *s++ = seed & 0xffffffffUL;
    for(mti = 1; mti < N; ++mti)
    {
        *s++ = (1812433253UL * (*r ^ (*r >> 30)) + mti) & 0xffffffffUL;
        r++;
    }
}

//! Generates a random number on [0,0xffffffff]-interval.
inline uint32_t
RandMT::next()
{
    uint32_t y;
    static const uint32_t mag01[2] = {0x0UL, MATRIX_A};
    // mag01[x] = x * MATRIX_A  for x=0,1

    // generate N words at one time
    if (mti >= N)
    {
        int kk;

        for (kk = 0; kk < N-M; kk++)
        {
            y = (mt[kk]&UPPER_MASK)| (mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        
        for (; kk < N-1; kk++)
        {
            y = (mt[kk]&UPPER_MASK)| (mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+ (M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        
        y = (mt[N-1]&UPPER_MASK)| (mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

} // namespace Rand

#endif // RAND_RAND_H
