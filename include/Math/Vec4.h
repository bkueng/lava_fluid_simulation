/*! \file Vec4.h
    \brief Contains the definition and implementation of the Vec4 class.
    \author Wojciech Jarosz
*/
#ifndef MATH_VEC4_H
#define MATH_VEC4_H

#include <Math/Fwd.h>
#include <Math/MathT.h>
#include <iostream>

namespace Math
{

//! A general 3D homogenous vector class.
/*!
    This class handles storing and manipulating homogenous 3D vectors.
*/
template <typename T>
class Vec4
{
public:
    T x, y, z, w;
    
    //-----------------------------------------------------------------------
    //! \name Constructors and assignment
	//@{
    //-----------------------------------------------------------------------
    //! Default constructor.
    /*!
        Call this constructor if you do not wish to initialize the vector to
        anything. It is more efficient than the others since it does nothing.
    */
    Vec4() {}

    //! Parameter constructor.
    /*!
        Initialize the vector with the three parameters.
        \param a Value to set the x component to.
        \param b Value to set the y component to.
        \param c Value to set the z component to.
    */
    Vec4(T a, T b, T c, T d) : x(a), y(b), z(c), w(d) {}

    //! Parameter constructor.
    /*!
        Initialize each component to \a a.
        \param a The value to initialize to.
    */
    explicit Vec4(T a) : x(a), y(a), z(a), w(a) {}
    
    //! Parameter constructor.
    /*!
        Initialize the vector from an array
        \param v The array to initialize from.
    */
    template <typename S>
    Vec4(const S* v) : x(v[0]), y(v[1]), z(v[2]), w(v[3]) {}
    
    //! Copy constructor.
    /*!
        Initialize the vector to be a copy of \a v.
        \param v The vector to create a copy of.
    */
    template <typename S>
    Vec4(const Vec4<S>& v) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(T(v.w)) {}
    
    //! Copy constructor.
    template <typename S>
    Vec4(const Vec3<S>& v, S w) : x(T(v.x)), y(T(v.y)), z(T(v.z)), w(w) {}

    //! Assignment operator.
    /*!
        Assigns the values from \a a to this Vec4.
    */
    const Vec4 & operator=(const Vec4& a) {x = a.x; y = a.y; z = a.z; w = a.w; return *this;}

    //! Assignment operator.
    /*!
        Assigns the values from \a a to this Vec4.
    */
    template <typename S>
    const Vec4 & operator=(const Vec3<S>& a) {x = T(a.x); y = T(a.y); z = T(a.z); return *this;}
    
    //! Assignment operator.
    /*!
        Sets all components of this Vec4 to \a a.
    */
    const Vec4 & operator=(T a) {x = y = z = w = a; return *this;}
    
    //! Assignment operator from C array
    template <typename S>
    const Vec4 & operator=(const S* v) {x = T(v[0]); y = T(v[1]); z = T(v[2]); w = T(v[3]); return *this;}
    //@}


    //-----------------------------------------------------------------------
    //! \name Casting operators.
	//@{
    //-----------------------------------------------------------------------
    //! Constant casting operator.
    /*!
        Casts this constant Vec4 to an array as a const T pointer.
    */
    const T* toArray() const {return (const T*)&x;}
    
    //! Casting operator.
    /*!
        Casts this Vec4 to an array as a T pointer.
    */
    T* toArray() {return (T*)&x;}
    //@}


    //-----------------------------------------------------------------------
    //! \name Element access and manipulation.
	//@{
    //-----------------------------------------------------------------------
    //! Access operator.        
    /*!
        Returns the ith component of the vector.
        \param i The component to return.
        \warning i must be either 0, 1, 2 or 3 in order to get expected results.
    */
    T & operator[](int i) {return(&x)[i];}
    
    //! Constant access operator.
    /*!
        Returns the ith component of a constant vector.
        \param i The component to return.
        \warning i must be either 0, 1, 2 or 3 in order to get expected results.
    */
    const T & operator[](int i) const {return(&x)[i];}

    void set(T a) {x = y = z = w = a;}
    void set(T a, T b, T c) {x = a; y = b; z = c;}
    void set(T a, T b, T c, T d) {x = a; y = b; z = c; w = d;}
    void set(const Vec4 v) {x = v.x; y = v.y; z = v.z; w = v.w;}
    template <typename S>
    void set(const Vec4<S>& v) {x = T(v.x); y = T(v.y); z = T(v.z); w = T(v.w);}
    //@}


    //-----------------------------------------------------------------------
    //! \name Addition.
	//@{
    //-----------------------------------------------------------------------
    //! Component-wise vector addition operator.
    Vec4 operator+(const Vec4& v) const
    {
        return Vec4(x + v.x, y + v.y, z + v.z, w + v.w);
    }
    
    //! Component-wise vector addition-assignment operator.
    const Vec4 & operator+=(const Vec4& v)
    {
        x += v.x; y += v.y; z += v.z; w += v.w; return *this;
    }

    //! Scalar addition-assignment operator.
    const Vec4 & operator+=(T a) {x += a; y += a; z += a; w += a; return *this;}
    //@}


    //-----------------------------------------------------------------------
    //! \name Subtraction.
	//@{
    //-----------------------------------------------------------------------
    //! Component-wise vector subtraction operator.
    Vec4 operator-(const Vec4& v) const
    {
        return Vec4(x - v.x, y - v.y, z - v.z, w - v.w);
    }
    
    //! Component-wise vector subtraction-assignment operator.
    const Vec4 & operator-=(const Vec4& v)
    {
        x -= v.x; y -= v.y; z -= v.z; w -= v.w; return *this;
    }
    
    //! Component-wise scalar subtraction assignment operator.
    const Vec4 & operator-=(T a) {x -= a; y -= a; z -= a; w -= w; return *this;}
    //@}


    //-----------------------------------------------------------------------
    //! \name Multiplication.
	//@{
    //-----------------------------------------------------------------------
    //! Scalar multiplication operator.
    Vec4 operator*(T a) const {return Vec4(x * a, y * a, z * a, w * a);}
    
    //! Component-wise vector multiplication operator.
    Vec4 operator*(const Vec4& v) const
    {
        return Vec4(x * v.x, y * v.y, z * v.z, z * v.z);
    }
    
    //! Scalar multiplication-assignment operator.
    const Vec4 & operator*=(T a) {x *= a; y *= a; z *= a; w *= a; return *this;}
    
    //! Component-wise vector multiplication-assignment operator.
    const Vec4 & operator*=(const Vec4& v)
    {
        x *= v.x; y *= v.y; z *= v.z; z *= v.z; w *= v.w; return *this;
    }
    
    //! Negation operator.
    Vec4 operator-() const {return Vec4(-x, -y, -z, -w);}
    const Vec4 & negate() {x = -x; y = -y; z = -z; w = -w; return *this;}
    //@}


    //-----------------------------------------------------------------------
    //! \name Division.
	//@{
    //-----------------------------------------------------------------------
    //! Scalar division operator.
    Vec4 operator/(T a) const {return (*this) * (T(1) / a);}
    
    //! Component-wise vector division operator.
    Vec4 operator/(const Vec4 & v) const
    {
        return Vec4(x / v.x, y / v.y, z / v.z, w / v.w);
    }
    
    //! Scalar division-assignment operator.
    const Vec4 & operator/=(T a) {*this *= T(1) / a; return *this;}
    
    //! Component-wise vector division-assignment operator.
    const Vec4 & operator/=(const Vec4 & v)
    {
        x /= v.x; y /= v.y; z /= v.z; w /= v.w; return *this;
    }
    //@}


    //-----------------------------------------------------------------------
    //! \name Equality.
	//@{
    //-----------------------------------------------------------------------
    //! Vector equivalence operator.
    /*!
        Tests to see if each component of \a v is equal to each component of
        this Vec4.
    */
    bool operator==(const Vec4 & v) const
    {
        return(v.x == x && v.y == y && v.z == z && v.w == w);
    }
    
    //! Vector difference operator.
    /*!
        Tests to see if any component is different between the two Vec4s.
    */
    bool operator!=(const Vec4 & v) const
    {
        return(v.x != x || v.y != y || v.z != z || v.w != w);
    }
    
    //! Compare two vectors and test if they are "approximately equal".
    /*!
        \returns true iff the coefficients of this and v are the same with
                 an absolute error of no more than e, i.e., for all i
        
                    abs(this[i] - v[i]) <= e
    */
    bool equalWithAbsError(const Vec4 &v, T e) const
    {
        for (int i = 0; i < 4; i++)
            if (!Math::equalWithAbsError((*this)[i], v[i], e))
                return false;
        return true;
    }
    
    //! Compare two vectors and test if they are "approximately equal".
    /*!
        \returns true iff the coefficients of this and v are the same with
                 an absolute error of no more than e, i.e., for all i
        
                    abs(this[i] - v[i]) <= e * abs(this[i])
    */
    bool equalWithRelError(const Vec4 &v, T e) const
    {
        for (int i = 0; i < 4; i++)
            if (!Math::equalWithRelError((*this)[i], v[i], e))
                return false;
        return true;
    }
    //@}


    //-----------------------------------------------------------------------
    //! \name Length, normalization, etc.
	//@{
    //-----------------------------------------------------------------------
    //! Length<sup>2</sup>.
    /*!
        Returns the geometric length<sup>2</sup> of the vector.
    */
    T length2() const {return dot(*this, *this);}
    
    //! Length.
    /*!
        Returns the geometric length of the vector.
    */
    T length() const {return MathT<T>::sqrt(length2());}
    
    //! Normalizes the vector and return its length.
    /*!
        Scales each component of the vector in order to get unit
        length without changing direction.
    
        \return The length of the vector prior to normalization.
    */
    T unitize()
    {
        T l = length();
        *this /= l;
        return l;
    }
    
    //! Normalize a vector and return a reference to it.
    /*!
        Scales each component of the vector in order to get unit
        length without changing direction.
    
        \return A reference to the vector.
    */
    const Vec4 & normalize()
    {
        return(*this /= length());
    }
    
    //! Return a normalized copy of the vector.
    Vec4 normalized() const
    {
        return(*this / length());
    }
    //@}

    static unsigned dimensions() {return 4;}
    typedef T    BaseType;
};


//! Multiply a scalar by a Vec4.
template <typename S, typename T>
inline Vec4<T>
operator*(S s, const Vec4<T>& v)
{
    return Vec4<T>(v.x * s, v.y * s, v.z * s, v.w * s);
}


//! The dot product of two Vec4s.
template <typename T>
inline T
dot(const Vec4<T>& a, const Vec4<T>& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}


//! Output to a stream.
/*!
    Writes \a a to the output stream \a out, and returns the modified stream.
    \return A copy of \a out with the \a a written to it.
    \sa operator>>(std::istream&, Vec4&)
    \relates Vec4
*/
template <typename T>
inline std::ostream &
operator<<(std::ostream& out, const Vec4<T>& v)
{
    return out << v.x << " " << v.y << " " << v.z << " " << v.w ;
}


//! Input from a stream.
/*!
    Reads \a a from the input stream \a in, and returns the modified stream.
    \return A copy of \a in with the \a a already read from it.
    \sa operator<<(std::ostream&, const Vec4&)
    \relates Vec4
*/
template <typename T>
inline std::istream &
operator>>(std::istream& in, Vec4<T>& v)
{
    return in >> v.x >> v.y >> v.z >> v.w ;
}

} //namespace rBmath

#endif // MATH_VEC4_H
