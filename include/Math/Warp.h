/*! \file Warp.h
    \brief Contains mappings between various useful domains.
    \author Wojciech Jarosz
*/
#ifndef MATH_WARP_H
#define MATH_WARP_H

#include <Math/Vec2.h>
#include <Math/Vec3.h>
#include <Math/Core.h>
#include <cmath>

namespace Math
{

//! Encapsulates various useful mappings
/*!
    This class provides functions to map from various useful domains, such as
    disks, spheres, hemispheres, etc.
    
    You may find the following useful for reference:
    
        Shirley, Pete. "Nonuniform Random Points Via Warping."
            Graphics Gems III. pp. 80--83
    and
        Dutré, Philip. "Global Illumination Compendium." 
            http://www.cs.kuleuven.ac.be/~phil/GI/
            
    among other sources.
        
*/
class Warp
{
public:
    //-----------------------------------------------------------------------
    //@{ \name Mappings to the 2D square.
    //-----------------------------------------------------------------------
    static void  uniformSquare(Vec2f* v, float s, float t);
    static float uniformSquarePdf();
    //@}


    //-----------------------------------------------------------------------
    //@{ \name Disk.
    //-----------------------------------------------------------------------
    static void  uniformDisk(Vec2f* v, float s, float t);
    static float uniformDiskPdf();

    static void  concentricUniformDisk(Vec2f* v, float s, float t);
    //@}
    
    static void  uniformCylinder(Vec3f* v, float s, float t);
    static float uniformCylinderPdf();


    //-----------------------------------------------------------------------
    //@{ \name Mappings to sections of the 3D sphere.
    //-----------------------------------------------------------------------
    static void  uniformSphere(Vec3f* v, float, float);
    static float uniformSpherePdf();

    static void  uniformSphericalCap(Vec3f* v, float, float, float);
    static float uniformSphericalCapPdf(float cosThetaMax);
    
    static void  uniformHemisphere(Vec3f* v, float s, float t);
    static float uniformHemispherePdf();
    
    static void  cosineHemisphere(Vec3f* v, float s, float t);
    static float cosineHemispherePdf(const Vec3f& v);
    
    static void  phongHemisphere(Vec3f* v, float s, float t, float n);
    static float phongHemispherePdf(const Vec3f& v, double n);

    static void  roughHemisphere(Vec3f* v, float s, float t, float alphag);
    //@}
    
    //-----------------------------------------------------------------------
    //@{ \name Mappings to sections of an arbitrary trinagle.
    //-----------------------------------------------------------------------
    /*! map uniform triangle with vertices a,b,c */
    static void uniformTriangle(Vec3f* v, float s, float t, const Vec3f& a, 
    		const Vec3f& b, const Vec3f& c);
    static float uniformTrianglePdf(const Vec3f& a, 
    		const Vec3f& b, const Vec3f& c);
};


inline void
Warp::uniformSquare(Vec2f* v, float s, float t)
{
    v->x = s;
    v->y = t;
}


inline float
Warp::uniformSquarePdf()
{
    return 1.0f;
}


inline void
Warp::uniformDisk(Vec2f* v, float s, float t)
{
	float theta = t * (M_PI * 2.f);
	float r = sqrtf(s);
	v->x = r*cosf(theta);
	v->y = r*sinf(theta);
}

inline void
Warp::concentricUniformDisk(Vec2f* v, float s, float t)
{
	//concentric uniform disk (less area distortions than uniformDisk,
	//(eg. uniformDisk mapping above has very long & narrow fields at
	//the border, which is bad for stratified sampling. areas should have a
	//square-like form))

	//taken from:  pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.
	
    // Map uniform random numbers to $[-1,1]^2$
    float sx = 2 * s - 1;
    float sy = 2 * t - 1;

	float r, theta;
    // Map square to $(r,\theta)$

    // Handle degeneracy at the origin
    if (sx == 0.0 && sy == 0.0) {
    	v->x = 0.0;
    	v->y = 0.0;
        return;
    }
    if (sx >= -sy) {
        if (sx > sy) {
            // Handle first region of disk
            r = sx;
            if (sy > 0.0) theta = sy/r;
            else          theta = 8.0f + sy/r;
        } else {
            // Handle second region of disk
            r = sy;
            theta = 2.0f - sx/r;
        }
    } else {
        if (sx <= sy) {
            // Handle third region of disk
            r = -sx;
            theta = 4.0f - sy/r;
        } else {
            // Handle fourth region of disk
            r = -sy;
            theta = 6.0f + sx/r;
        }
    }
    theta *= M_PI / 4.f;
    v->x = r * cosf(theta);
    v->y = r * sinf(theta);
}

inline float
Warp::uniformDiskPdf()
{
    return 1.f/M_PI;
}


inline void
Warp::uniformCylinder(Vec3f* v, float s, float t)
{
	v->z = s*2.f - 1.f;
	float theta = t * (M_PI*2.f);
	v->x = sinf(theta);
	v->y = cosf(theta);
}


inline float
Warp::uniformCylinderPdf()
{
	//cylinder height is 2.f
	return 1.0f/(2.f*2.f*M_PI);
}


//! Samples a unit sphere uniformily.
/*!
    Generates a direction on the unit sphere proportional to solid angle.
    
    Uses Archimedes Theorem to sample the cylinder and then projects back
    onto the sphere.
*/
inline void
Warp::uniformSphere(Vec3f* v, float s, float t)
{
	v->z = s*2.f - 1.f;
	float theta = t * (M_PI*2.f);
	float r = sqrtf(1.f - v->z*v->z);
	v->x = r*sinf(theta);
	v->y = r*cosf(theta);
}


inline float
Warp::uniformSpherePdf()
{
	return 1.0f/(4.f*M_PI);
}


inline void
Warp::uniformSphericalCap(Vec3f* v, float s, float t, float cosThetaMax)
{
	//cosThetaMax in [-1,1]
	v->z = s*(cosThetaMax+1.f) - 1.f;
	float theta = t * (M_PI*2.f);
	float r = sqrtf(1.f - v->z*v->z);
	v->x = r*sinf(theta);
	v->y = r*cosf(theta);
}


inline float
Warp::uniformSphericalCapPdf(float cosThetaMax)
{
	return 1.0f/(2.f*M_PI * (cosThetaMax+1.f));
}


//! Samples a unit hemisphere uniformily.
/*!
    Generates a direction on the unit hemisphere uniformily distributed wrto
    solid angle.
*/
inline void
Warp::uniformHemisphere(Vec3f* v, float s, float t)
{
	v->z = s;
	float theta = t * (M_PI*2.f);
	float r = sqrtf(std::max((dfloat)0., (dfloat)1. - v->z*v->z));
	v->x = r*cosf(theta);
	v->y = r*sinf(theta);
}


inline float
Warp::uniformHemispherePdf()
{
	return uniformSpherePdf() * 2.f;
}


//! Samples a cosine-weighted hemisphere.
/*!
    Generates a direction on the unit hemisphere distributed proportional
    to cosine-weighted solid angle.
    
    The technique used is to just use spherical coordinates directly.
*/
inline void
Warp::cosineHemisphere(Vec3f* v, float s, float t)
{
	//uniform disk
	float theta = t * (M_PI * 2.f);
	float r = sqrtf(s);
	v->x = r*cosf(theta);
	v->y = r*sinf(theta);
	//project to hemisphere
	v->z = sqrtf(std::max((dfloat)0., (dfloat)1. - s));
}

inline float
Warp::cosineHemispherePdf(const Vec3f& v)
{
	//cos(phi) / PI, cos(phi) = sqrt(1 - (x^2 + y^2)) / r = z / r, r=1
	return std::max(Math::FEQ_EPS2, (float)(v.z / M_PI));
}


//! Samples a phong-weighted hemisphere.
/*!
    Generates a direction on the unit hemisphere distributed proportional
    to cosine^n-weighted solid angle.
*/
inline void
Warp::phongHemisphere(Vec3f* v, float s, float t, float n)
{
	// TODO: use fastPow ??
	
	v->z = pow(1.f-s, 1.f/(n+1.f));
	float r = sqrtf(std::max((dfloat)0., (dfloat)1.-v->z*v->z));
	float phi = t * (M_PI*2.f);
	v->x = r*cos(phi);
	v->y = r*sin(phi);
}


inline float
Warp::phongHemispherePdf(const Vec3f& v, double n)
{
	return (n+1.)/(2.*M_PI) *  std::max(Math::FEQ_EPS2, (float)(pow((double)v.z, n)));
}


//! Samples a unit hemisphere according to rough dielectric sampling.
//  this is described in the paper: Microfacet Models for Refraction through Rough Surfaces, Bruce Walter et al.
inline void
Warp::roughHemisphere(Vec3f* v, float s, float t, float alphag)
{
	v->z = 1.f / sqrtf(std::max(0.f, alphag*alphag * s / (1.f-s) + 1.f));
	float theta = t * (M_PI*2.f);
	float r = sqrtf(std::max((dfloat)0., (dfloat)1. - v->z*v->z));
	v->x = r*cosf(theta);
	v->y = r*sinf(theta);
}


inline void
Warp::uniformTriangle(Vec3f* v, float s, float t, const Vec3f& a, 
		const Vec3f& b, const Vec3f& c) {
	/*
	if(s+t > 1.f) {
		s = (FEQ_EPS/10.f + 1.f) - s;
		t = (FEQ_EPS/10.f + 1.f) - t;
	}
	*v = a + (b-a)*s + (c-a)*t;
	//*/
	s = 1.f - sqrtf(1.f - s);
	t = (1.f - s) * t;
	*v = a + (b-a)*s + (c-a)*t;
}

inline float
Warp::uniformTrianglePdf(const Vec3f& a, 
    		const Vec3f& b, const Vec3f& c) {
	float triangle_area = cross(b-a, c-a).length()/2.;
	return 1.f / triangle_area;
}

} // namespace Math

#endif // MATH_WARP_H
