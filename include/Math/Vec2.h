/*! \file Vec2.h
    \brief Contains the definition and implementation of the Vec2 class.
    \author Wojciech Jarosz
*/
#ifndef MATH_VEC2_H
#define MATH_VEC2_H

#include <Math/Fwd.h>
#include <Math/MathT.h>
#include <iostream>

namespace Math
{

//! A general 2D vector class.
/*!
    This class handles storing and manipulating 2D vectors.
*/
template <typename T>
class Vec2
{
public:
    T x, y;
    
    //-----------------------------------------------------------------------
    //! \name Constructors and assignment
	//@{
    //-----------------------------------------------------------------------
    //! Default constructor.
    /*!
        Call this constructor if you do not wish to initialize the vector to
        anything. It is more efficient than the others since it does nothing.
    */
    Vec2() {}

    //! Parameter constructor.
    /*!
        Initialize the vector with the three parameters.
        \param a Value to set the x component to.
        \param b Value to set the y component to.
    */
    Vec2(T a, T b) : x(a), y(b) {}

    //! Parameter constructor.
    /*!
        Initialize each component to \a a.
        \param a The value to initialize to.
    */
    explicit Vec2(T a) : x(a), y(a) {}
    
    //! Parameter constructor.
    /*!
        Initialize the vector from an array of doubles
        \param v The array of doubles to initialize from.
    */
    Vec2(const T* v) : x(v[0]), y(v[1]) {}
    
    //! Copy constructor.
    /*!
        Initialize the vector to be a copy of \a v.
        \param v The vector to create a copy of.
    */
    template <typename S>
    Vec2(const Vec2<S>& v) : x(T(v.x)), y(T(v.y)) {}

    //! Assignment operator.
    /*!
        Assigns the values from \a a to this Vec2.
    */
    const Vec2 & operator=(const Vec2& a) {x = a.x; y = a.y; return *this;}
    
    //! Assignment operator.
    /*!
        Sets all components of this Vec2 to \a a.
    */
    const Vec2 & operator=(T a) {x = y = a; return *this;}
    //@}


    //-----------------------------------------------------------------------
    //! \name Casting operators.
	//@{
    //-----------------------------------------------------------------------
    //! Constant casting operator.
    /*!
        Casts this constant Vec2 to an array as a const T pointer.
    */
    const T* toArray() const {return (const T*)&x;}
    
    //! Casting operator.
    /*!
        Casts this Vec2 to an array as a T pointer.
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
        \warning i must be either 0, 1, or 2 in order to get expected results.
    */
    T & operator[](int i) {return(&x)[i];}
    
    //! Constant access operator.
    /*!
        Returns the ith component of a constant vector.
        \param i The component to return.
        \warning i must be either 0, 1, or 2 in order to get expected results.
    */
    const T & operator[](int i) const {return(&x)[i];}

    void set(T a) {x = y = a;}
    void set(T a, T b) {x = a; y = b;}
    void set(const Vec2 v) {x = v.x; y = v.y;}
    template <typename S>
    void set(const Vec2<S>& v) {x = T(v.x); y = T(v.y);}
    //@}


    //-----------------------------------------------------------------------
    //! \name Addition.
	//@{
    //-----------------------------------------------------------------------
    //! Component-wise vector addition operator.
    Vec2 operator+(const Vec2& v) const {return Vec2(x + v.x, y + v.y);}
    
    //! Component-wise vector addition-assignment operator.
    const Vec2 & operator+=(const Vec2& v) {x += v.x; y += v.y; return *this;}

    //! Scalar addition-assignment operator.
    const Vec2 & operator+=(T a) {x += a; y += a; return *this;}
    //@}


    //-----------------------------------------------------------------------
    //! \name Subtraction.
	//@{
    //-----------------------------------------------------------------------
    //! Component-wise vector subtraction operator.
    Vec2 operator-(const Vec2& v) const {return Vec2(x - v.x, y - v.y);}
    
    //! Component-wise vector subtraction-assignment operator.
    const Vec2 & operator-=(const Vec2& v) {x -= v.x; y -= v.y; return *this;}
    
    //! Component-wise scalar subtraction assignment operator.
    const Vec2 & operator-=(T a) {x -= a; y -= a; return *this;}
    //@}


    //-----------------------------------------------------------------------
    //! \name Multiplication.
	//@{
    //-----------------------------------------------------------------------
    //! Scalar multiplication operator.
    Vec2 operator*(T a) const {return Vec2(x * a, y * a);}
    
    //! Component-wise vector multiplication operator.
    Vec2 operator*(const Vec2& v) const {return Vec2(x * v.x, y * v.y);}
    
    //! Scalar multiplication-assignment operator.
    const Vec2 & operator*=(T a) {x *= a; y *= a; return *this;}
    
    //! Component-wise vector multiplication-assignment operator.
    const Vec2 & operator*=(const Vec2& v) {x *= v.x; y *= v.y; return *this;}
    
    //! Negation operator.
    Vec2 operator-() const {return Vec2(-x, -y);}
    const Vec2 & negate() {x = -x; y = -y; return *this;}
    //@}


    //-----------------------------------------------------------------------
    //! \name Division.
	//@{
    //-----------------------------------------------------------------------
    //! Scalar division operator.
    Vec2 operator/(T a) const
    {
        T inv = T(1)/a;
        return Vec2(x * inv, y * inv);
    }
    
    //! Component-wise vector division operator.
    Vec2 operator/(const Vec2 & v) const {return Vec2(x / v.x, y / v.y);}
    
    //! Scalar division-assignment operator.
    const Vec2 & operator/=(T a)
    {
        T inv = T(1)/a;
        x *= inv; y *= inv;
        return *this;
    }
    
    //! Component-wise vector division-assignment operator.
    const Vec2 & operator/=(const Vec2 & v) {x /= v.x; y /= v.y; return *this;}
    //@}


    //-----------------------------------------------------------------------
    //! \name Equality.
	//@{
    //-----------------------------------------------------------------------
    //! Vector equivalence operator.
    /*!
        Tests to see if each component of \a v is equal to each component of
        this Vec2.
    */
    bool operator==(const Vec2 & v) const {return(v.x == x && v.y == y);}
    
    //! Vector difference operator.
    /*!
        Tests to see if any component is different between the two Vec2s.
    */
    bool operator!=(const Vec2 & v) const {return(v.x != x || v.y != y);}
    
    //! Compare two vectors and test if they are "approximately equal".
    /*!
        \returns true iff the coefficients of this and v are the same with
                 an absolute error of no more than e, i.e., for all i
        
                    abs(this[i] - v[i]) <= e
    */
    bool equalWithAbsError(const Vec2 &v, T e) const
    {
        for (int i = 0; i < 2; i++)
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
    bool equalWithRelError(const Vec2 &v, T e) const
    {
        for (int i = 0; i < 2; i++)
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
    const Vec2 & normalize()
    {
        return(*this /= length());
    }
    
    //! Return a normalized copy of the vector.
    Vec2 normalized() const
    {
        return(*this / length());
    }

    //! Hat operator.
    void hat()
    {
        T tmp = x;
        x = -y;
        y = tmp;
    }
    
    //! Hat operator.
    Vec2 hatted() const {return Vec2(-y, x);}
    //@}

    static unsigned dimensions() {return 2;}
    typedef T    BaseType;
};


//! Multiply a scalar by a Vec2.
/*!
    \return The a Vec2 with each component multiplied by \a s.
    \relates Vec2
*/
template <typename S, typename T>
inline Vec2<T>
operator*(S s, const Vec2<T>& v)
{
    return Vec2<T>(v.x * s, v.y * s);
}


//! The dot product of two Vec2s.
/*!
    This is a convenience function that simply performs a.dot(b).
    \return The dot product of \a a and \a b.
    \relates Vec2
*/
template <typename T, typename S>
inline T
dot(const Vec2<T>& a, const Vec2<S>& b)
{
    return a.x * b.x + a.y * b.y;
}


//! Output to a stream.
/*!
    Writes \a a to the output stream \a out, and returns the modified stream.
    \return A copy of \a out with the \a a written to it.
    \sa operator>>(std::istream&, Vec2&)
    \relates Vec2
*/
template <typename T>
inline std::ostream &
operator<<(std::ostream& out, const Vec2<T>& v)
{
    return out << v.x << " " << v.y ;
}


//! Input from a stream.
/*!
    Reads \a a from the input stream \a in, and returns the modified stream.
    \return A copy of \a in with the \a a already read from it.
    \sa operator<<(std::ostream&, const Vec2&)
    \relates Vec2
*/
template <typename T>
inline std::istream &
operator>>(std::istream& in, Vec2<T>& v)
{
    return in >> v.x >> v.y ;
}

} // namespace Math

#endif // MATH_VEC2_H
