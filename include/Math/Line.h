/*! \file Ray.h
    \brief Contains the definition of Line class.
    \author Wojciech Jarosz
*/
#ifndef MATH_LINE_H_INCLUDED
#define MATH_LINE_H_INCLUDED

#include <Math/Fwd.h>
#include <Math/MathT.h>
#include <iostream>

namespace Math
{

template <typename Vec>
class Line
{
public:
    typedef typename Vec::BaseType BaseType;

    Vec o;            //!< The origin of the ray
    Vec d;            //!< The direction of the ray

    Line() {}
    Line(const Vec & o, const Vec & d) : o(o), d(d) {}

    Vec operator()(BaseType t) const {return o + t*d;}

    BaseType computeT(const Vec& p) const
    {
        if (MathT<BaseType>::fabs(d.x) >= MathT<BaseType>::fabs(d.y) &&
            MathT<BaseType>::fabs(d.x) >= MathT<BaseType>::fabs(d.z))
            return (p.x - o.x)/d.x;
        else if (MathT<BaseType>::fabs(d.y) >= MathT<BaseType>::fabs(d.z))
            return (p.y - o.y)/d.y;
        else
            return (p.z - o.z)/d.z;
    }
};


template <typename Vec>
inline std::ostream&
operator<<(std::ostream &o, const Line<Vec> &r)
{
    return o << "(" << r.o << ") + t(" << r.d << ")";
}

} // namespace Math

#endif // MATH_LINE_H_INCLUDED
