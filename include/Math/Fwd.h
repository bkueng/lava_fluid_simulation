/*! \file Fwd.h
    \brief Forward declarations of all classes in the Math module.
    \author Wojciech Jarosz
*/
#ifndef MATH_FWD_FWD_H
#define MATH_FWD_FWD_H

#include "global.h"

/*!
	\namespace Math
	Namespace for the math module. Includes classes for vectors, matrices, colors, etc.
 */
namespace Math
{

template <typename T> class Vec2;
typedef Vec2<int>           Vec2i;
typedef Vec2<unsigned>      Vec2u;
typedef Vec2<float>         Vec2f;
typedef Vec2<double>        Vec2d;
typedef Vec2<short>         Vec2s;
typedef Vec2<char>          Vec2c;
typedef Vec2<unsigned char> Vec2uc;

template <typename T> class Vec3;
typedef Vec3<int>           Vec3i;
typedef Vec3<unsigned>      Vec3u;
typedef Vec3<dfloat>         Vec3f;
typedef Vec3<double>        Vec3d;
typedef Vec3<short>         Vec3s;
typedef Vec3<char>          Vec3c;
typedef Vec3<unsigned char> Vec3uc;

template <typename T> class Vec4;
typedef Vec4<int>           Vec4i;
typedef Vec4<unsigned>      Vec4u;
typedef Vec4<float>         Vec4f;
typedef Vec4<double>        Vec4d;
typedef Vec4<short>         Vec4s;
typedef Vec4<char>          Vec4c;
typedef Vec4<unsigned char> Vec4uc;
	
template <typename T> class Color3;
typedef Color3<int>				Color3i;
typedef Color3<unsigned>		Color3u;
typedef Color3<float>			Color3f;
typedef Color3<double>			Color3d;
typedef Color3<unsigned short>	Color3s;
typedef Color3<unsigned char>	Color3c;
typedef Color3<float>			C3f;
typedef Color3<double>			C3d;
typedef Color3<unsigned short>	C3s;
typedef Color3<unsigned char>	C3c;

template <typename T> class Color4;
typedef Color4<int>				Color4i;
typedef Color4<unsigned>		Color4u;
typedef Color4<float>			Color4f;
typedef Color4<double>			Color4d;
typedef Color4<unsigned short>	Color4s;
typedef Color4<unsigned char>	Color4c;
typedef Color4<unsigned char>	Color4uc;
typedef Color3<float>			C3f;
typedef Color3<double>			C3d;
typedef Color3<unsigned short>	C3s;
typedef Color3<unsigned char>	C3c;
typedef unsigned int			PackedColor;
class RGBE;

template <typename T> class Mat44;
typedef Mat44<float>  Mat44f;
typedef Mat44<double> Mat44d;

template <typename T> class MovingMat44;
typedef MovingMat44<float>  MovingMat44f;
typedef MovingMat44<double> MovingMat44d;

template <typename T> class Quat;
typedef Quat<float>  Quatf;
typedef Quat<double> Quatd;

template <typename T> class ONB;
typedef ONB<float>  ONBf;
typedef ONB<double> ONBd;

template <typename Vec> class Box;
typedef Box<Vec2i> Box2i;
typedef Box<Vec2u> Box2u;
typedef Box<Vec2f> Box2f;
typedef Box<Vec2d> Box2d;
typedef Box<Vec3i> Box3i;
typedef Box<Vec3u> Box3u;
typedef Box<Vec3f> Box3f;
typedef Box<Vec3d> Box3d;
typedef Box<Vec4i> Box4i;
typedef Box<Vec4f> Box4f;
typedef Box<Vec4d> Box4d;

template <typename Vec> class Line;
typedef Line<Vec2i> Line2i;
typedef Line<Vec2u> Line2u;
typedef Line<Vec2f> Line2f;
typedef Line<Vec2d> Line2d;
typedef Line<Vec3i> Line3i;
typedef Line<Vec3u> Line3u;
typedef Line<Vec3f> Line3f;
typedef Line<Vec3d> Line3d;
typedef Line<Vec4i> Line4i;
typedef Line<Vec4u> Line4u;
typedef Line<Vec4f> Line4f;
typedef Line<Vec4d> Line4d;

} // namespace Math

#endif // MATH_FWD_FWD_H
