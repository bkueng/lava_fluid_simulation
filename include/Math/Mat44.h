/*! \file Mat44.h
    \brief Contains the definition and implementation of the Mat44 class.
    \author Wojciech Jarosz
*/
#ifndef MATH_MAT44_H
#define MATH_MAT44_H

#include <Math/Vec3.h>
#include <Math/Vec4.h>
#include <Math/MathT.h>
#include <Math/FastMath.h>
#include <Math/LimitsT.h>

#include <iomanip>

namespace Math
{

//! A class for storing and manipulating 3D transformation matrices.
/*!
    Matrix vector multiplication is done as follows:
\verbatim
    [ x'  y'  z'  w'] = [ Ax  Bx  Cx  Dx ]   [ x ]
                        [ Ay  By  Cy  Dy ] * [ y ]
                        [ Az  Bz  Cz  Dz ]   [ z ]
                        [ Aw  Bw  Cw  Dw ]   [ w ]
\endverbatim
    or if the matrix is expressed in ABCD vector notation:

\verbatim
    [ x'  y'  z'  w'] = [ A  B  C  D ] * [ X ]
\endverbatim

	A rigid transformation will look like:

\verbatim
    [ Ax  Bx  Cx  Dx ]
    [ Ay  By  Cy  Dy ]
    [ Az  Bz  Cz  Dz ]
    [ 0   0   0   1  ]
\endverbatim
 
	where the A, B, C vectors form an ortho-normal basis, and the D vector is
    the translation.
    
    This class is compatible with OpenGL and stores the matrix coefficients
    to allow the use of glLoadMatrix, glMultMatrix, etc. Because of
    peculiarities in GL, the matrices are actually stored transposed, so that
    the A, B, C, D column vectors are each contiguous in memory. As a C array
    this would be:

\verbatim
    { { Ax, Ay, Az, Aw },
      { Bx, By, Bz, Bw },
      { Cx, Cy, Cz, Cw },
      { Dx, Dy, Dz, Dw } }
\endverbatim
*/
template <typename T>
class Mat44
{
private:
    Vec4<T> cols[4];

public:

    //-----------------------------------------------------------------------
    //! \name Constructors and assignment.
	//@{
    //-----------------------------------------------------------------------
    static Mat44 I()
    {
        return Mat44(1);
    }
    Mat44() {}
    explicit Mat44(T s);
    Mat44(const Vec4<T> &a,
          const Vec4<T> &b,
          const Vec4<T> &c,
          const Vec4<T> &d);
    Mat44(const Vec3<T> &aVec,
          const Vec3<T> &bVec,
          const Vec3<T> &cVec,
          const Vec3<T> &dVec);
    Mat44(T ax, T bx, T cx, T dx,
          T ay, T by, T cy, T dy,
          T az, T bz, T cz, T dz,
          T aw, T bw, T cw, T dw);
    Mat44(const T a[4][4]);

    const Mat44 & operator=(T s);
    const Mat44 & operator=(const Mat44 & m);
    void  makeIdentity() {*this = I();}
    //@}


// #ifndef _WIN32
    //-----------------------------------------------------------------------
    //! \name Casting operators.
	//@{
    //-----------------------------------------------------------------------
//     operator       T*()       {return (T*)this;}
//     operator const T*() const {return (const T*)this;}
    //! Constant casting operator.
    /*!
        Casts this constant Vec3 to an array as a const T pointer.
    */
    const T* toArray() const {return (const T*)&cols[0];}
    
    //! Casting operator.
    /*!
        Casts this Mat44 to an array as a T pointer.
    */
    T* toArray() {return (T*)&cols[0];}
    //@}
// #endif // _WIN32
    
    
    //-----------------------------------------------------------------------
    //! \name Element access.
	//@{
    //-----------------------------------------------------------------------
    T *       operator[](int i)              {return (T *)&cols[i];}
    const T * operator[](int i) const        {return (const T *)&cols[i];}
    T &       operator()(int i, int j)       {return cols[i][j];}
    const T & operator()(int i, int j) const {return cols[i][j];}
    
    
    Vec3<T> A() const {return Vec3<T>(cols[0][0], cols[0][1], cols[0][2]);}
    Vec3<T> B() const {return Vec3<T>(cols[1][0], cols[1][1], cols[1][2]);}
    Vec3<T> C() const {return Vec3<T>(cols[2][0], cols[2][1], cols[2][2]);}
    Vec3<T> D() const {return Vec3<T>(cols[3][0], cols[3][1], cols[3][2]);}
    
    void setA(T i, T j, T k) {cols[0][0] = i; cols[0][1] = j; cols[0][2] = k;}
    void setA(T i, T j, T k, T l) {setA(i,j,k); cols[0][3] = l;}
    void setA(const Vec3<T> & a) {setA(a.x, a.y, a.z);}
    void setA(const Vec4<T> & a) {cols[0] = a;}
    void setB(T i, T j, T k) {cols[1][0] = i; cols[1][1] = j; cols[1][2] = k;}
    void setB(T i, T j, T k, T l) {setB(i,j,k); cols[1][3] = l;}
    void setB(const Vec3<T> & b) {setB(b.x, b.y, b.z);}
    void setB(const Vec4<T> & b) {setB(b.x, b.y, b.z, b.w);}
    void setC(T i, T j, T k) {cols[2][0] = i; cols[2][1] = j; cols[2][2] = k;}
    void setC(T i, T j, T k, T l) {setC(i,j,k); cols[2][3] = l;}
    void setC(const Vec3<T> & c) {setC(c.x, c.y, c.z);}
    void setC(const Vec4<T> & c) {setC(c.x, c.y, c.z, c.w);}
    void setD(T i, T j, T k) {cols[3][0] = i; cols[3][1] = j; cols[3][2] = k;}
    void setD(T i, T j, T k, T l) {setD(i,j,k); cols[3][3] = l;}
    void setD(const Vec3<T> & d) {setD(d.x, d.y, d.z);}
    void setD(const Vec4<T> & d) {setD(d.x, d.y, d.z, d.w);}
    
    const T * col(int i) const {return (T*)&cols[i];}
    T * col(int i) {return (const T*)&cols[i];}
    const Vec3<T> & colAsVec3(int i) const {return *(Vec3<T>*)col(i);}
    const Vec4<T> & colAsVec4(int i) const {return cols[i];}
    const Vec3<T>   rowAsVec3(int i) const {return Vec3<T>(cols[0][i], cols[1][i], cols[2][i]);}
    const Vec4<T>   rowAsVec4(int i) const {return Vec4<T>(cols[0][i], cols[1][i], cols[2][i], cols[3][i]);}
    void setRow(int i, const Vec3<T> & r) {cols[0][i] = r[0]; cols[1][i] = r[1]; cols[2][i] = r[2];}
    void setRow(int i, const Vec4<T> & r) {cols[0][i] = r[0]; cols[1][i] = r[1]; cols[2][i] = r[2]; cols[3][i] = r[3];}
    void setCol(int i, const Vec3<T> & c) {cols[i] = c;}
    void setCol(int i, const Vec4<T> & c) {cols[i] = c;}
    void set(T ax, T bx, T cx, T dx,
             T ay, T by, T cy, T dy,
             T az, T bz, T cz, T dz,
             T aw, T bw, T cw, T dw);
    //@}
    
    
    //-----------------------------------------------------------------------
    //! \name Equality.
	//@{
    //-----------------------------------------------------------------------
    bool operator==(const Mat44 &v) const;
    bool operator!=(const Mat44 &v) const;
    bool equalWithAbsError(const Mat44<T> &v, T e) const;
    bool equalWithRelError(const Mat44<T> &v, T e) const;
    //@}
    
    
    //-----------------------------------------------------------------------
    //! \name Component-wise addition and subtraction.
	//@{
    //-----------------------------------------------------------------------
    Mat44         operator+ (const Mat44 & m) const;
    Mat44         operator+ (T s) const;
    const Mat44 & operator+=(const Mat44 & m);
    
    Mat44         operator- (const Mat44 & m) const;
    Mat44         operator- (T s) const;
    const Mat44 & operator-=(const Mat44 & m);
    //@}
    
    
    //-----------------------------------------------------------------------
    //! \name Component-wise multiplication by -1
	//@{
    //-----------------------------------------------------------------------
    Mat44         operator-() const;
    const Mat44 & negate();
    //@}
    
    
    //-----------------------------------------------------------------------
    //! \name Component-wise multiplication.
	//@{
    //-----------------------------------------------------------------------
    Mat44         operator* (T s) const;
    const Mat44 & operator*=(T s);
    //@}
    
    
    //-----------------------------------------------------------------------
    //! \name Matrix-matrix multiplication.
	//@{
    //-----------------------------------------------------------------------
    Mat44         operator* (const Mat44 & m) const;
    const Mat44 & operator*=(const Mat44 & m);
    
    const Mat44 & multiply(const Mat44 &a, const Mat44 &b);
    const Mat44 & multiply3x3(const Mat44 &a, const Mat44 &b);
    //@}
    
    
    //-----------------------------------------------------------------------
    //! \name Matrix-vector multiplication and transformation.
	//@{
    //-----------------------------------------------------------------------
    template <typename S>
    void multiplyP(Vec3<S> *dst, const Vec3<S> &src) const;
    template <typename S>
    void multiplyP4x3(Vec3<S> *dst, const Vec3<S> &src) const;
    template <typename S>
    void multiplyD(Vec3<S> *dst, const Vec3<S> &src) const;
    template <typename S>
    void multiplyN(Vec3<S> *dst, const Vec3<S> &src) const;

    template <typename S>
    void multiplyP(Vec3<S> *v) const {multiplyP(v, *v);}
    template <typename S>
    void multiplyP4x3(Vec3<S> *v) const {multiplyP4x3(v, *v);}
    template <typename S>
    void multiplyD(Vec3<S> *v) const {multiplyD(v, *v);}
    template <typename S>
    void multiplyN(Vec3<S> *v) const {multiplyN(v, *v);}
    //@}
    
    
    //-----------------------------------------------------------------------
    //! \name Component-wise division.
	//@{
    //-----------------------------------------------------------------------
    Mat44         operator/ (T s) const;
    const Mat44 & operator/=(T s);
    //@}
    
    
    //-----------------------------------------------------------------------
    //! \name Matrix transpose.
	//@{
    //-----------------------------------------------------------------------
    const Mat44 & transpose();
    Mat44         transposed() const;
    //@}
    
    
    //-----------------------------------------------------------------------
    //! \name Matrix inversion.
	//@{
    //-----------------------------------------------------------------------
	//! Invert this matrix using the Gauss-Jordan method.
    const Mat44 & gjInvert();
	//! Compute the matrix inverse using the Gauss-Jordan method, and return the result.
    Mat44         gjInverse() const;
	//! Invert this matrix using determinants.
    const Mat44 & invert();
	//! Compute the matrix inverse using determinants, and return the result.
    Mat44         inverse() const;
	//! Invert the matrix assuming it is a rigid transformation.
    const Mat44 & rigidInvert();
	//! Compute the matrix inverse assuming it is a rigid transformation, and return the result
    Mat44         rigidInverse() const;
    //@}

    
    //-----------------------------------------------------------------------
    //! \name Matrix norms and determinant.
	//@{
    //-----------------------------------------------------------------------
    T             norm1() const;
    T             normInf() const;
    T             normF() const;
    T             determinant() const;
    T             trace() const;
    //@}

    
    //-----------------------------------------------------------------------
    //! \name Matrix square-root, logarithm, and exponential.
	//@{
    //-----------------------------------------------------------------------
    Mat44         sqrt(T e = T(1e-20)) const;
    Mat44         log(int i = 6, T e = T(1e-20)) const;
    Mat44         exp() const;
    //@}

    //-----------------------------------------------------------------------
    //! \name Create Rotations - these set the current matrix to a rotation.
	//@{
    //-----------------------------------------------------------------------
    template <typename S>
    const Mat44<T> &setRotateX(S angle);

    template <typename S>
    const Mat44<T> & setRotateY(S angle);

    template <typename S>
    const Mat44<T> &setRotateZ(S angle);
    
    template <typename S>
    const Mat44 & setAxisAngle(const Vec3<S> & ax, S ang);
    
    template <typename S>
    const Mat44 & setEulerAngles(const Vec3<S> & r);
    //@}


    //-----------------------------------------------------------------------
    //! \name Rotations - these rotate the current matrix by a rotation.
	//@{
    //-----------------------------------------------------------------------
    template <typename S>
    const Mat44<T> & rotateX(S angle);

    template <typename S>
    const Mat44<T> & rotateY(S angle);

    template <typename S>
    const Mat44<T> & rotateZ(S angle);

    template <typename S>
    const Mat44<T> & rotateAxisAngle(const Vec3<T> & axis, S angle);
    
    template <typename S>
    const Mat44<T> & rotateEulerAngles(const Vec3<S> & r);
    
    const Mat44<T> & rotateTo(const Vec3<T> &from, const Vec3<T> &to);
    //@}
    
    
    //-----------------------------------------------------------------------
    //! \name Scaling.
	//@{
    //-----------------------------------------------------------------------
    template <typename S>
    const Mat44 & setScale(S s) {setScale(s, s, s); return *this;}
    template <typename S>
    const Mat44 & setScale(S sx, S sy, S sz);
    template <typename S>
    const Mat44 & setScale(const Vec3<S> & s)
        {setScale(s.x, s.y, s.z); return *this;}
    template <typename S>
    const Mat44 & scale(const Vec3<S>& s);
    //@}
    

    //-----------------------------------------------------------------------
    //! \name Translations.
	//@{
    //-----------------------------------------------------------------------
    template <typename S>
    const Mat44 &   setTranslation(const Vec3<S> &t);
    //@}
    
    
    //-----------------------------------------------------------------------
    //! \name OpenGL style functions and their inverses.
	//@{
    //-----------------------------------------------------------------------
    void frustum(T l, T r, T b, T t, T n, T f);
    void invFrustum(T l, T r, T b, T t, T n, T f);
    void perspective(T yFov, T aspect, T nDist, T fDist);
    void invPerspective(T yFov, T aspect, T nDist, T fDist);
    void ortho(const Vec3<T>& min, const Vec3<T>& max);
    void invOrtho(const Vec3<T>& min, const Vec3<T>& max);
    void viewport(int ww, int wh);
    void invViewport(int ww, int wh);
    void lookAt(const Vec3<T>& eye,
                const Vec3<T>& lookAtPt,
                const Vec3<T>& viewUp);
    void invLookAt(const Vec3<T>& eye,
                   const Vec3<T>& lookAtPt,
                   const Vec3<T>& viewUp);
    //@}

    static unsigned xSize() {return 4;}
    static unsigned ySize() {return 4;}
    typedef T BaseType;
};


template <typename T>
inline
Mat44<T>::Mat44(T s)
{
    *this = s;
}


template <typename T>
inline
Mat44<T>::Mat44(const T m[4][4])
{
#ifdef VECMAT_LOOPS
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            cols[i][j] = m[i][j];
#else
    cols[0] = m[0];
    cols[1] = m[1];
    cols[2] = m[2];
    cols[3] = m[3];
#endif
}


template <typename T>
inline
Mat44<T>::Mat44(const Vec4<T> &a,
                const Vec4<T> &b,
                const Vec4<T> &c,
                const Vec4<T> &d)
{
    cols[0] = a;
    cols[1] = b;
    cols[2] = c;
    cols[3] = d;
}


template <typename T>
inline
Mat44<T>::Mat44(const Vec3<T> &a,
                const Vec3<T> &b,
                const Vec3<T> &c,
                const Vec3<T> &d)
{
    cols[0] = a; cols[0][3] = 0;
    cols[1] = b; cols[1][3] = 0;
    cols[2] = c; cols[2][3] = 0;
    cols[3] = d; cols[3][3] = 1;
}


template <typename T>
inline
Mat44<T>::Mat44(T ax, T bx, T cx, T dx,
                T ay, T by, T cy, T dy,
                T az, T bz, T cz, T dz,
                T aw, T bw, T cw, T dw)
{
    set(ax, bx, cx, dx,
        ay, by, cy, dy,
        az, bz, cz, dz,
        aw, bw, cw, dw);
}


template <typename T>
inline void
Mat44<T>::set(T ax, T bx, T cx, T dx,
              T ay, T by, T cy, T dy,
              T az, T bz, T cz, T dz,
              T aw, T bw, T cw, T dw)
{
    cols[0].set(ax, ay, az, aw);
    cols[1].set(bx, by, bz, bw);
    cols[2].set(cx, cy, cz, cw);
    cols[3].set(dx, dy, dz, dw);
}


template <class T>
inline const Mat44<T> &
Mat44<T>::operator=(const Mat44 & m)
{
#ifdef VECMAT_LOOPS
    for (int i = 0; i < 4; ++i)
        cols[i] = m[i];
#else
    cols[0] = m[0];
    cols[1] = m[1];
    cols[2] = m[2];
    cols[3] = m[3];
#endif
    return *this;
}



template <typename T>
inline const Mat44<T> &
Mat44<T>::operator=(T s)
{
    set(s, 0, 0, 0,
        0, s, 0, 0,
        0, 0, s, 0,
        0, 0, 0, 1);
    return *this;
}


//----------------------------------
// Addition and subtraction
//----------------------------------


template <typename T>
inline const Mat44<T> &
Mat44<T>::operator+=(const Mat44<T> &v)
{
    cols[0] += v.cols[0];
    cols[1] += v.cols[1];
    cols[2] += v.cols[2];
    cols[3] += v.cols[3];
    
    return *this;
}


template <typename T>
inline Mat44<T>
Mat44<T>::operator+(const Mat44<T> &v) const
{
    return Mat44(cols[0] + v.cols[0],
                 cols[1] + v.cols[1],
                 cols[2] + v.cols[2],
                 cols[3] + v.cols[3]);
}


template <typename T>
inline const Mat44<T> &
Mat44<T>::operator-=(const Mat44<T> &v)
{
    cols[0] -= v.cols[0];
    cols[1] -= v.cols[1];
    cols[2] -= v.cols[2];
    cols[3] -= v.cols[3];
    
    return *this;
}


template <typename T>
inline Mat44<T>
Mat44<T>::operator-(const Mat44<T> &v) const
{
    return Mat44(cols[0] - v.cols[0],
                 cols[1] - v.cols[1],
                 cols[2] - v.cols[2],
                 cols[3] - v.cols[3]);
}


template <typename T>
Mat44<T>
Mat44<T>::operator-() const
{
    return Mat44(-cols[0], -cols[1], -cols[2], -cols[3]);
}


template <typename T>
const Mat44<T> &
Mat44<T>::negate()
{
    cols[0].negate();
    cols[1].negate();
    cols[2].negate();
    cols[3].negate();

    return *this;
}


//---------------------------------------------------------------
// Vector-times-matrix multiplication operators
//---------------------------------------------------------------

/*!
    Treats the input vector as a position.
*/
template <typename T, typename S>
Vec3<S>
operator*(const Mat44<T> &m, const Vec3<S> &v)
{
    Vec3<S> ret;
    m.multiplyP(&ret, v);
    return ret;
}


/*!
    Treats the input vector as a position.
*/
template <typename S, typename T>
const Vec3<S> &
operator*=(Vec3<S> &v, const Mat44<T> &m)
{
    m.multiplyP(&v);
    return v;
}


template <typename T, typename S>
Vec4<S>
operator*(const Mat44<T> &m, const Vec4<S> &v)
{
#ifdef VECMAT_LOOPS
    Vec4<S> ret;
    for (int i = 0; i < 4; ++i)
        ret[i] = v[0]*m[0][i] + v[1]*m[1][i] + v[2]*m[2][i] + v[3]*m[3][i];
            
    return ret;
#else
    return Vec4<T>(v.x*m[0][0] + v.y*m[1][0] + v.z*m[2][0] + v.w*m[3][0],
                   v.x*m[0][1] + v.y*m[1][1] + v.z*m[2][1] + v.w*m[3][1],
                   v.x*m[0][2] + v.y*m[1][2] + v.z*m[2][2] + v.w*m[3][2],
                   v.x*m[0][3] + v.y*m[1][3] + v.z*m[2][3] + v.w*m[3][3]);
#endif
}


template <typename S, typename T>
const Vec4<S> &
operator*=(Vec4<S> &v, const Mat44<T> &m)
{
    v = m * v;
    return v;
}


template <typename T>
inline const Mat44<T> &
Mat44<T>::operator*=(const Mat44<T> &v)
{
    Mat44 tmp;
    tmp.multiply(*this, v);
    *this = tmp;
    return *this;
}


template <typename T>
inline Mat44<T>
Mat44<T>::operator*(const Mat44<T> &v) const
{
    Mat44 tmp;
    tmp.multiply(*this, v);
    return tmp;
}


//! Transform a position vector by this matrix.
/*!
    \pre
    \a s is a 3 dimensional position vector, i.e. it specifies a location.

    \post
    \a d contains the result of transforming \a src by this matrix.
*/
template <typename T>
template <typename S>
inline void
Mat44<T>::multiplyP(Vec3<S> *d, const Vec3<S> &s) const
{
#ifdef VECMAT_LOOPS
    S r[4];
    for (int i = 0; i < 4; ++i)
        r[i] = s[0]*cols[0][i] + s[1]*cols[1][i] + s[2]*cols[2][i] + cols[3][i];

    for (int i = 0; i < 3; ++i)
        (*d)[i] = r[i]/r[3];
#else
    S a = s.x*cols[0][0] + s.y*cols[1][0] + s.z*cols[2][0] + cols[3][0];
    S b = s.x*cols[0][1] + s.y*cols[1][1] + s.z*cols[2][1] + cols[3][1];
    S c = s.x*cols[0][2] + s.y*cols[1][2] + s.z*cols[2][2] + cols[3][2];
    S w = s.x*cols[0][3] + s.y*cols[1][3] + s.z*cols[2][3] + cols[3][3];

    d->x = a / w;
    d->y = b / w;
    d->z = c / w;
#endif
}

//! Transform a position vector by affine 4x3 portion (no projection) of this matrix.
/*!
    \pre
    \a s is a 3 dimensional position vector, i.e. it specifies a location.

    \post
    \a d contains the result of transforming \a src by this matrix.
*/
template <typename T>
template <typename S>
inline void
Mat44<T>::multiplyP4x3(Vec3<S> *d, const Vec3<S> &s) const
{
#ifdef VECMAT_LOOPS
    S r[3];
    for (int i = 0; i < 3; ++i)
        r[i] = s[0]*cols[0][i] + s[1]*cols[1][i] + s[2]*cols[2][i] + cols[3][i];

    for (int i = 0; i < 3; ++i)
        (*d)[i] = r[i];
#else
    S a = s.x*cols[0][0] + s.y*cols[1][0] + s.z*cols[2][0] + cols[3][0];
    S b = s.x*cols[0][1] + s.y*cols[1][1] + s.z*cols[2][1] + cols[3][1];
    S c = s.x*cols[0][2] + s.y*cols[1][2] + s.z*cols[2][2] + cols[3][2];

    d->x = a;
    d->y = b;
    d->z = c;
#endif
}


//! Transform a direction vector by this matrix.
/*!
    \pre
    \a s is a 3 dimensional direction vector, i.e. its homogenous
    coordinate is 0

    \post
    \a d contains the result of transforming \a src by this matrix.
*/
template <typename T>
template <typename S>
inline void
Mat44<T>::multiplyD(Vec3<S> *d, const Vec3<S> &s) const
{
#ifdef VECMAT_LOOPS
    S r[3];
    for (int i = 0; i < 3; ++i)
        r[i] = s[0]*cols[0][i] + s[1]*cols[1][i] + s[2]*cols[2][i];

    for (int i = 0; i < 3; ++i)
        (*d)[i] = r[i];
#else
    d->set(s.x * cols[0][0] + s.y * cols[1][0] + s.z * cols[2][0],
           s.x * cols[0][1] + s.y * cols[1][1] + s.z * cols[2][1],
           s.x * cols[0][2] + s.y * cols[1][2] + s.z * cols[2][2]);
#endif
}


//! Transform a direction vector by the transpose of this matrix.
/*!
    \pre
    \a s is a 3 dimensional direction vector, i.e. its homogenous
    coordinate is 0

    \post
    \a d contains the result of transforming \a src by the transpose of this
    matrix.
    
    \remarks
    This function is useful when transforming normals. It eliminates the
    explicit computation of the transpose, saving some time
*/
template <typename T>
template <typename S>
inline void
Mat44<T>::multiplyN(Vec3<S> *d, const Vec3<S> &s) const
{
#ifdef VECMAT_LOOPS
    S r[3];
    for (int i = 0; i < 3; ++i)
        r[i] = s[0]*cols[i][0] + s[1]*cols[i][1] + s[2]*cols[i][2];
        
    for (int i = 0; i < 3; ++i)
        (*d)[i] = r[i];
#else
    d->set(s.x * cols[0][0] + s.y * cols[0][1] + s.z * cols[0][2],
           s.x * cols[1][0] + s.y * cols[1][1] + s.z * cols[1][2],
           s.x * cols[2][0] + s.y * cols[2][1] + s.z * cols[2][2]);
#endif
}


//----------------------------------
// Translations
//----------------------------------

template <typename T>
template <typename S>
inline const Mat44<T> &
Mat44<T>::setTranslation(const Vec3<S> &t)
{
    makeIdentity();
    cols[3] = t;
    return *this;
}


//----------------------------------
// Transpose and inverse
//----------------------------------

template <typename T>
inline const Mat44<T> &
Mat44<T>::transpose()
{
    std::swap(cols[1][0], cols[0][1]);
    std::swap(cols[2][0], cols[0][2]);
    std::swap(cols[2][1], cols[1][2]);
    std::swap(cols[3][0], cols[0][3]);
    std::swap(cols[3][1], cols[1][3]);
    std::swap(cols[3][2], cols[2][3]);
    return *this;
}


template <typename T>
inline Mat44<T>
Mat44<T>::transposed() const
{
    return Mat44(cols[0][0], cols[0][1], cols[0][2], cols[0][3],
                 cols[1][0], cols[1][1], cols[1][2], cols[1][3],
                 cols[2][0], cols[2][1], cols[2][2], cols[2][3],
                 cols[3][0], cols[3][1], cols[3][2], cols[3][3]);
}


template <typename T>
bool
Mat44<T>::operator==(const Mat44 & m) const
{
    return cols[0][0] == m[0][0] &&
           cols[0][1] == m[0][1] &&
           cols[0][2] == m[0][2] &&
           cols[0][3] == m[0][3] &&
           cols[1][0] == m[1][0] &&
           cols[1][1] == m[1][2] &&
           cols[1][2] == m[1][2] &&
           cols[1][3] == m[1][3] &&
           cols[2][0] == m[2][0] &&
           cols[2][1] == m[2][1] &&
           cols[2][2] == m[2][2] &&
           cols[2][3] == m[2][3] &&
           cols[3][0] == m[3][0] &&
           cols[3][1] == m[3][1] &&
           cols[3][2] == m[3][2] &&
           cols[3][3] == m[3][3];
}


template <typename T>
bool
Mat44<T>::operator!=(const Mat44 & m) const
{
    return cols[0][0] != m[0][0] ||
           cols[0][1] != m[0][1] ||
           cols[0][2] != m[0][2] ||
           cols[0][3] != m[0][3] ||
           cols[1][0] != m[1][0] ||
           cols[1][1] != m[1][2] ||
           cols[1][2] != m[1][2] ||
           cols[1][3] != m[1][3] ||
           cols[2][0] != m[2][0] ||
           cols[2][1] != m[2][1] ||
           cols[2][2] != m[2][2] ||
           cols[2][3] != m[2][3] ||
           cols[3][0] != m[3][0] ||
           cols[3][1] != m[3][1] ||
           cols[3][2] != m[3][2] ||
           cols[3][3] != m[3][3];
}




//! Compare two matrices and test if they are "approximately equal".
/*!
    \returns true iff the coefficients of this and m are the same with
             an absolute error of no more than e, i.e., for all i, j

                abs(this[i][j] - m[i][j]) <= e
*/
template <typename T>
bool
Mat44<T>::equalWithAbsError(const Mat44<T> &m, T e) const
{
    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
        if (!Math::equalWithAbsError((*this)[i][j], m[i][j], e))
        return false;

    return true;
}


//! Compare two matrices and test if they are "approximately equal".
/*!
    \returns true iff the coefficients of this and m are the same with
             an relative error of no more than e, i.e., for all i, j

                abs(this[i][j] - m[i][j]) <= e * abs (this[i][j])
*/
template <typename T>
bool
Mat44<T>::equalWithRelError(const Mat44<T> &m, T e) const
{
    for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
        if (!Math::equalWithRelError((*this)[i][j], m[i][j], e))
        return false;

    return true;
}


template <typename T>
const Mat44<T> &
Mat44<T>::operator*=(T a)
{
#ifdef VECMAT_LOOPS
    for (int i = 0; i < 4; ++i)
        cols[i] *= a;
#else
    cols[0] *= a;
    cols[1] *= a;
    cols[2] *= a;
    cols[3] *= a;
#endif

    return *this;
}


template <typename T>
Mat44<T>
Mat44<T>::operator*(T a) const
{
    return Mat44(cols[0] * a,
                 cols[1] * a,
                 cols[2] * a,
                 cols[3] * a);
}


template <typename T>
Mat44<T> operator*(T a, const Mat44<T> & m)
{
    return m * a;
}


template <typename T>
const Mat44<T> &
Mat44<T>::operator/=(T a)
{
    *this *= T(1)/a;

    return *this;
}


template <typename T>
Mat44<T>
Mat44<T>::operator/(T a) const
{
    return (*this) * (T(1)/a);
}


template <typename T>
const Mat44<T> &
Mat44<T>::multiply(const Mat44 &a, const Mat44 &b)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        {
#ifdef VECMAT_LOOPS
            cols[i][j] = 0;
            for (int k = 0; k < 4; ++k)
                cols[i][j] += a[k][j] * b[i][k]; 
#else
            cols[i][j] = a[0][j]*b[i][0] + a[1][j]*b[i][1] +
                      a[2][j]*b[i][2] + a[3][j]*b[i][3];
#endif
        }
    
    return *this;
}

/*!
	Set this matrix to be the 3x3 matrix product of \a a and \a b.
	Only touches the upper 3x3 portion of all matrices.
 */
template <typename T>
const Mat44<T> &
Mat44<T>::multiply3x3(const Mat44 &a, const Mat44 &b)
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
#ifdef VECMAT_LOOPS
            cols[i][j] = 0;
            for (int k = 0; k < 3; ++k)
                cols[i][j] += a[k][j] * b[i][k]; 
#else
            cols[i][j] = a[0][j]*b[i][0] + a[1][j]*b[i][1] + a[2][j]*b[i][2];
#endif
        }
    
    return *this;
}


template <typename T>
const Mat44<T> &
Mat44<T>::gjInvert()
{
    *this = gjInverse();
    return *this;
}

/*!
	This is significantly slower than \func inverse(), but the results may
	be slightly more accurate.
 
	Taking the inverse of a singular matrix returns the identity.
 */
template <typename T>
Mat44<T>
Mat44<T>::gjInverse() const
{
    int i, j, k;
    Mat44 s(1);
    Mat44 t(*this);

    // Forward elimination
    for (i = 0; i < 3 ; i++)
    {
        int pivot = i;
    
        T pivotsize = t[i][i];
    
        if (pivotsize < 0)
            pivotsize = -pivotsize;
    
        for (j = i + 1; j < 4; j++)
        {
            T tmp = t[j][i];
    
            if (tmp < 0)
                tmp = -tmp;
    
            if (tmp > pivotsize)
            {
                pivot = j;
                pivotsize = tmp;
            }
        }
    
        if (pivotsize == 0)
            return Mat44();
    
        if (pivot != i)
        {
            for (j = 0; j < 4; j++)
            {
                T tmp;
        
                tmp = t[i][j];
                t[i][j] = t[pivot][j];
                t[pivot][j] = tmp;
        
                tmp = s[i][j];
                s[i][j] = s[pivot][j];
                s[pivot][j] = tmp;
            }
        }
    
        for (j = i + 1; j < 4; j++)
        {
            T f = t[j][i] / t[i][i];
    
            for (k = 0; k < 4; k++)
            {
                t[j][k] -= f * t[i][k];
                s[j][k] -= f * s[i][k];
            }
        }
    }

    // Backward substitution
    for (i = 3; i >= 0; --i)
    {
        T f;
    
        if ((f = t[i][i]) == 0)
            return Mat44();
    
        for (j = 0; j < 4; j++)
        {
            t[i][j] /= f;
            s[i][j] /= f;
        }
    
        for (j = 0; j < i; j++)
        {
            f = t[j][i];
    
            for (k = 0; k < 4; k++)
            {
                t[j][k] -= f * t[i][k];
                s[j][k] -= f * s[i][k];
            }
        }
    }

    return s;
}


template <typename T>
const Mat44<T> &
Mat44<T>::invert()
{
    *this = inverse();
    return *this;
}

/*!
	This is significantly faster than \func gjInverse(), but the results may
	be slightly less accurate.

	Taking the inverse of a singular matrix returns the identity.
*/
template <typename T>
Mat44<T>
Mat44<T>::inverse() const
{
    if (cols[0][3] != 0 || cols[1][3] != 0 || cols[2][3] != 0 || cols[3][3] != 1)
        return gjInverse();

    Mat44 s(cols[1][1] * cols[2][2] - cols[2][1] * cols[1][2],
            cols[2][0] * cols[1][2] - cols[1][0] * cols[2][2],
            cols[1][0] * cols[2][1] - cols[2][0] * cols[1][1],
            0,

            cols[2][1] * cols[0][2] - cols[0][1] * cols[2][2],
            cols[0][0] * cols[2][2] - cols[2][0] * cols[0][2],
            cols[2][0] * cols[0][1] - cols[0][0] * cols[2][1],
            0,

            cols[0][1] * cols[1][2] - cols[1][1] * cols[0][2],
            cols[1][0] * cols[0][2] - cols[0][0] * cols[1][2],
            cols[0][0] * cols[1][1] - cols[1][0] * cols[0][1],
            0,

            0, 0, 0, 1);

    T r = cols[0][0] * s[0][0] + cols[0][1] * s[1][0] + cols[0][2] * s[2][0];

    if (MathT<T>::fabs(r) >= 1)
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                s[i][j] /= r;
    }
    else
    {
        T mr = MathT<T>::fabs(r) / Limits<T>::smallest();

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (mr > MathT<T>::fabs(s[i][j]))
                    s[i][j] /= r;
                else
                    return Mat44();
    }

    s[3][0] = -cols[3][0] * s[0][0] - cols[3][1] * s[1][0] - cols[3][2] * s[2][0];
    s[3][1] = -cols[3][0] * s[0][1] - cols[3][1] * s[1][1] - cols[3][2] * s[2][1];
    s[3][2] = -cols[3][0] * s[0][2] - cols[3][1] * s[1][2] - cols[3][2] * s[2][2];

    return s;
}


template <typename T>
const Mat44<T> &
Mat44<T>::rigidInvert()
{
    *this = rigidInverse();
    return *this;
}

/*!
	This specialized method is extremely fast, but only works for rigid
	transformations.
 */
template <typename T>
Mat44<T>
Mat44<T>::rigidInverse() const
{
    Mat44 tmp;
    tmp.set(cols[0][0], cols[0][1], cols[0][2], -dot(D(), A()),
            cols[1][0], cols[1][1], cols[1][2], -dot(D(), B()),
            cols[2][0], cols[2][1], cols[2][2], -dot(D(), C()),
                  0,       0,       0,             1);

    return tmp;
}


template <typename T>
T
Mat44<T>::determinant() const
{
    return  cols[0][3] * cols[1][2] * cols[2][1] * cols[3][0] -
            cols[0][2] * cols[1][3] * cols[2][1] * cols[3][0] -
            cols[0][3] * cols[1][1] * cols[2][2] * cols[3][0] +
            cols[0][1] * cols[1][3] * cols[2][2] * cols[3][0] +
            cols[0][2] * cols[1][1] * cols[2][3] * cols[3][0] -
            cols[0][1] * cols[1][2] * cols[2][3] * cols[3][0] -
            cols[0][3] * cols[1][2] * cols[2][0] * cols[3][1] +
            cols[0][2] * cols[1][3] * cols[2][0] * cols[3][1] +
            cols[0][3] * cols[1][0] * cols[2][2] * cols[3][1] -
            cols[0][0] * cols[1][3] * cols[2][2] * cols[3][1] -
            cols[0][2] * cols[1][0] * cols[2][3] * cols[3][1] +
            cols[0][0] * cols[1][2] * cols[2][3] * cols[3][1] +
            cols[0][3] * cols[1][1] * cols[2][0] * cols[3][2] -
            cols[0][1] * cols[1][3] * cols[2][0] * cols[3][2] -
            cols[0][3] * cols[1][0] * cols[2][1] * cols[3][2] +
            cols[0][0] * cols[1][3] * cols[2][1] * cols[3][2] +
            cols[0][1] * cols[1][0] * cols[2][3] * cols[3][2] -
            cols[0][0] * cols[1][1] * cols[2][3] * cols[3][2] -
            cols[0][2] * cols[1][1] * cols[2][0] * cols[3][3] +
            cols[0][1] * cols[1][2] * cols[2][0] * cols[3][3] +
            cols[0][2] * cols[1][0] * cols[2][1] * cols[3][3] -
            cols[0][0] * cols[1][2] * cols[2][1] * cols[3][3] -
            cols[0][1] * cols[1][0] * cols[2][2] * cols[3][3] +
            cols[0][0] * cols[1][1] * cols[2][2] * cols[3][3];
}


template <typename T>
T
Mat44<T>::trace() const
{
    return  cols[0][0] + cols[1][1] + cols[2][2] + cols[3][3];
}


/*!
    One norm
    \return    maximum column sum.
*/
template <typename T>
T
Mat44<T>::norm1() const
{
    T f = 0;
    for (int j = 0; j < 4; j++)
    {
        T s = 0;
        for (int i = 0; i < 4; i++)
            s += MathT<T>::fabs(cols[i][j]);
        f = std::max(f, s);
    }
    return f;
}


/*!
    Infinity norm
    \return    maximum row sum.
*/
template <typename T>
T
Mat44<T>::normInf() const
{
    T f = 0;
    for (int i = 0; i < 4; i++)
    {
        T s = 0;
        for (int j = 0; j < 4; j++)
            s += MathT<T>::fabs(cols[i][j]);
        f = std::max(f, s);
    }
    return f;
}


/*!
    Frobenius norm
    \return    sqrt of sum of squares of all elements.
*/
template <typename T>
T
Mat44<T>::normF() const
{
    T f = 0;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            f += cols[i][j] * cols[i][j];
    return MathT<T>::sqrt(f);
}

template <typename T>
Mat44<T>
Mat44<T>::sqrt(T e) const
{
    const Mat44 & A = *this;
    Mat44 X = *this;
    if (A.determinant() < 0)
        // Matrix has no real roots
        return Mat44::I();

    Mat44 Y = Mat44::I();
    Mat44 iX, iY;
    int maxIter = 20;
    while ((maxIter-- > 0) && ((X*X - A).normF() > e))
    {
        iX = X.inverse();
        iY = Y.inverse();
        X = (X + iY) * T(0.5);
        Y = (Y + iX) * T(0.5);
    }
    return X;
}


template <typename T>
Mat44<T>
Mat44<T>::log(int iterations, T e) const
{
    int k = 0;
    Mat44 A = *this;
    Mat44 I = Mat44::I();

    while ((A - I).normInf() > T(0.5) && (k < 50))
    {
        A = A.sqrt(e);
        ++k;
    }

    A = I - A;
    Mat44 Z = A;
    Mat44 X = A;
    int i = 1;
    while ((Z.normF() > e) && (i < iterations))
    {
        Z *= A;
        ++i;
        X += Z * (T(1.0) / i);
    }
    X *= -T(1u << k);
    return X;
}


template <typename T>
Mat44<T>
Mat44<T>::exp() const
{
    int q = 6;
    int j = std::max(0, 1 + iLog2(normInf()));

    Mat44 A = (*this) * (T(1.0)/(1u << j));

    Mat44 D(1), N(1), X(1);
    float c = 1;
    for (int k = 1; k <= q; ++k)
    {
        c = c * (q - k + 1) / (k * (2 * q - k + 1));
        X = A * X;
        N += X * c;
        float signedC = (k % 2) ? -c : c;
        D += X * signedC;
    }
    X = D.inverse() * N;

    for (int i = 0; i < j; ++i)
        X *= X;
    return X;
}


template <typename T>
template <typename S>
const Mat44<T> &
Mat44<T>::setScale(S sx, S sy, S sz)
{
    set(sx,  0,  0,  0,
         0, sy,  0,  0,
         0,  0, sz,  0,
         0,  0,  0,  1);
    return *this;
}


template <typename T>
template <typename S>
const Mat44<T> &
Mat44<T>::scale(const Vec3<S> &s)
{
    cols[0] *= s[0];
    cols[1] *= s[1];
    cols[2] *= s[2];

    return *this;
}


//----------------------------------
// Rotations
//----------------------------------

template <typename T>
template <typename S>
const Mat44<T> &
Mat44<T>::rotateX(S angle)
{
    Mat44 mtx;
    mtx.setRotateX(angle);
    *this *= mtx;

    return *this;
}


template <typename T>
template <typename S>
const Mat44<T> &
Mat44<T>::rotateY(S angle)
{
    Mat44 mtx;
    mtx.setRotateY(angle);
    *this *= mtx;

    return *this;
}


template <typename T>
template <typename S>
const Mat44<T> &
Mat44<T>::rotateZ(S angle)
{
    Mat44 mtx;
    mtx.setRotateZ(angle);
    *this *= mtx;

    return *this;
}


template <typename T>
template <typename S>
const Mat44<T> &
Mat44<T>::rotateAxisAngle(const Vec3<T> & axis, S t)
{
    Mat44 mtx;
    mtx.setAxisAngle(axis, t);
    *this *= mtx;

    return *this;
}


template <typename T>
template <typename S>
const Mat44<T> &
Mat44<T>::rotateEulerAngles(const Vec3<S> & r)
{
    Mat44 mtx;
    mtx.setEulerAngles(r);
    *this *= mtx;

    return *this;
}


template <typename T>
const Mat44<T> &
Mat44<T>::rotateTo(const Vec3<T> &from, const Vec3<T> &to)
{
    Vec3<T> axis;
    T angle;

    if (dot(from, to) > T(0.999999)) return *this;  // Already lined up

    axis = cross(from, to);
    if (axis.length2() < T(0.0000001))
    {
        if (from.x > T(-0.9) && from.x < T(0.9))
            axis.set(1, 0, 0);
        else if (from.y > T(-0.9) && from.y < T(0.9))
            axis.set(0, 1, 0);
        else
            axis.set(0, 0, 1);
        angle = M_PI;
    }
    else
    {
        angle = MathT<T>::acos(dot(from, to));
        axis.normalize();
    }

    Mat44 mtx;
    mtx.setAxisAngle(axis, angle);
    *this *= mtx;
    
    return *this;
}


template <typename T>
template <typename S>
const Mat44<T> &
Mat44<T>::setRotateX(S angle)
{
    T cost = MathT<T>::cos(angle);
    T sint = MathT<T>::sin(angle);

    set(1,    0,     0,  0,
        0, cost, -sint,  0,
        0, sint,  cost,  0,
        0,    0,     0,  1);
    
    return *this;
}


template <typename T>
template <typename S>
const Mat44<T> &
Mat44<T>::setRotateY(S angle)
{
    T cost = MathT<T>::cos(angle);
    T sint = MathT<T>::sin(angle);

    set( cost,    0,  sint,  0,
            0,    1,     0,  0,
        -sint,    0,  cost,  0,
            0,    0,     0,  1);
    
    return *this;
}


template <typename T>
template <typename S>
const Mat44<T> &
Mat44<T>::setRotateZ(S angle)
{
    T cost = MathT<T>::cos(angle);
    T sint = MathT<T>::sin(angle);

    set(cost, -sint,  0,  0,
        sint,  cost,  0,  0,
           0,     0,  1,  0,
           0,     0,  0,  1);
    
    return *this;
}


template <typename T>
template <typename S>
const Mat44<T> &
Mat44<T>::setEulerAngles(const Vec3<S> & r)
{
    S sx = MathT<T>::sin(r.x), cx = MathT<T>::cos(r.x);
    S sy = MathT<T>::sin(r.y), cy = MathT<T>::cos(r.y);
    S sz = MathT<T>::sin(r.z), cz = MathT<T>::cos(r.z);

    set(cy*cz, sx*sy*cz-cx*sz, cx*sy*cz+sx*sz, 0,
        cy*sz, sx*sy*sz+cx*cz, cx*sy*sz-sx*cz, 0,
          -sy,          sx*cy,          cx*cy, 0,
            0,              0,              0, 1);

    return *this;
}


template <typename T>
template <typename S>
const Mat44<T> &
Mat44<T>::setAxisAngle(const Vec3<S> & v, S a)
{
    S ct  = MathT<T>::cos(a);
    S st  = MathT<T>::sin(a);
    S omc = 1.0f - ct;

    set(omc*v.x*v.x + ct    , omc*v.x*v.y - st*v.z, omc*v.x*v.z + st*v.y, 0,
        omc*v.x*v.y + st*v.z, omc*v.y*v.y + ct    , omc*v.y*v.z - st*v.x, 0,
        omc*v.x*v.z - st*v.y, omc*v.y*v.z + st*v.x, omc*v.z*v.z + ct    , 0,
        0                   , 0                   , 0                   , 1);

    return *this;
}


/*!
	Sets the matrix to the frustum matrix (see also \func invFrustum()).

	\param l, r Specify the coordinates for the left and right vertical clipping planes.
	\param b, t Specify the coordinates for the bottom and top horizontal clipping planes.
	\param n, f Specify the distances to the near and far depth clipping planes. Both distances must be positive.

	Same as glFrustum() : Perspective transformation matrix defined by a
	truncated pyramid viewing frustum that starts at the origin(eye)
	going in the -Z axis direction(viewer looks down -Z)
	with the four pyramid sides passing through the sides of a window
	defined through x=l, x=r, y=b, y=t on the viewplane at z=-n.
	The top and bottom of the pyramid are truncated by the near and far
	planes at z=-n and z=-f.
 
	A 4-vector(x,y,z,w) inside this frustum transformed by this matrix will
	have x,y,z values in the range [-w,w]. Homogeneous clipping is applied to
	restrict(x,y,z) to this range.
 
	Later, a perspective divide by w will result in an NDC
	coordinate 3-vector of the form(x/w,y/w,z/w,w/w)=(x',y',z',1) where
	x', y', and z' all are in the range [-1,1]. Perspectively divided z'
	will be in [-1,1] with -1 being the near plane and +1 being the far.
 */
template <typename T>
void
Mat44<T>::frustum(T l, T r, T b, T t, T n, T f)
{
    T inv1 = T(1) / (r-l);
    T inv2 = T(1) / (t-b);
    T inv3 = T(1) / (f-n);

    set((2*n)*inv1,          0,  (r+l)*inv1,             0,
                 0, (2*n)*inv2,  (t+b)*inv2,             0,
                 0,          0, -(f+n)*inv3, (-2*f*n)*inv3,
                 0,          0,          -1,             0);
}

	
//!	Sets the matrix to the inverse frustum matrix (see also \func frustum())
template<class T>
void
Mat44<T>::invFrustum(T l, T r, T b, T t, T n, T f)
{
    T inv1 = T(1) / (2*n);
    T inv2 = T(1) / (2*f*n);

    set((r-l)*inv1,          0,           0, (r+l)*inv1,
                 0, (t-b)*inv1,           0, (t+b)*inv1,
                 0,          0,           0,         -1,
                 0,          0, -(f-n)*inv2, (f+n)*inv2);
}

//! Sets the matrix to a perspective projection.
/*!
	Same as gluPerspective : calls frustum()
 */
template <typename T>
void
Mat44<T>::perspective(T yFov, T aspect, T nDist, T fDist)
{
    T wT = MathT<T>::tan(radians(yFov) * 0.5f) * nDist,
      wB = -wT;
    T wR = wT * aspect,
      wL = -wR;
    frustum(wL, wR, wB, wT, nDist, fDist);
}

	
//! Sets the matrix to the inverse perspective projection.
/*!
	Inverse of perspective() : calls invFrustum()
 */
template <typename T>
void
Mat44<T>::invPerspective(T yFov, T aspect, T nDist, T fDist)
{
    T wT = MathT<T>::tan(radians(yFov) * 0.5f) * nDist,
      wB = -wT;
    T wR = wT * aspect,
      wL = -wR;
    invFrustum(wL, wR, wB, wT, nDist, fDist);
}


/*!
	Set this matrix to a viewport transformation.

	Same effect as calling gluViewport(0, 0, ww, hh):
	Given Window width and height in pixels transforms the x,y,z NDC values in
	[-1,1] to (x',y') pixel values and normalized z' in [0,1]
	(near and far respectively).
 
	\param ww, hh Specify the width and height of the viewport in pixels.
 */
template <typename T>
void
Mat44<T>::viewport(int ww, int wh)
{
    T ww2 = T(ww) * T(0.5),
      wh2 = T(wh) * T(0.5);

    set(ww2,   0,   0, ww2,
          0, wh2,   0, wh2,
          0,   0, 0.5, 0.5,
          0,   0,   0,   1);
}

//! Inverse of viewport()
template <typename T>
void
Mat44<T>::invViewport(int ww, int wh)
{
    T ww2 = T(2) / T(ww),
      wh2 = T(2) / T(wh);

    set(ww2,   0,   0,  -1,
          0, wh2,   0,  -1,
          0,   0,   2,  -1,
          0,   0,   0,   1);
}


//! Set this matrix to a look-at matrix.
/*!
	Creates a viewing matrix derived from an eye point, a reference point
	indicating the center of the scene, and an up vector. Same as gluLookat().
 
	\param eye Specifies the position of the eye point.
	\param center Specifies the position of the reference point.
	\param up Specifies the direction of the up vector.

	The matrix maps the reference point to the negative z axis and the eye point
	to the origin. When a typical projection matrix is used, the center of the
	scene therefore maps to the center of the viewport. Similarly, the direction
	described by the UP vector projected onto the viewing plane is mapped to the
	positive y axis so that it points upward in the viewport. The UP vector must
	not be parallel to the line of sight from the eye point to the reference point.
 
 */
template <typename T>
void
Mat44<T>::lookAt(const Vec3<T>& eye,
                  const Vec3<T>& center,
                  const Vec3<T>& up)
{
    Vec3<T> z  = eye - center;  z.normalize();
    Vec3<T> x  = cross(up, z); x.normalize();
    Vec3<T> y  = cross(z, x);  y.normalize();
    Vec3<T> tr = -eye;
    
    set(x.x, x.y, x.z, dot(x, tr),
        y.x, y.y, y.z, dot(y, tr),
        z.x, z.y, z.z, dot(z, tr),
          0,   0,   0,           1);
}

	
//! Inverse of lookAt()
template <typename T>
void
Mat44<T>::invLookAt(const Vec3<T>& eye,
                     const Vec3<T>& center,
                     const Vec3<T>& up)
{
    Vec3<T> z = eye - center;  z.normalize();
    Vec3<T> x = cross(up, z); x.normalize();
    Vec3<T> y = cross(z, x);  y.normalize();
    
    set(x.x, y.x, z.x, eye.x,
        x.y, y.y, z.y, eye.y,
        x.z, y.z, z.z, eye.z,
          0,   0,   0,     1);
}


//! Set this matrix to a orthographic projection.
/*!
*/
template <typename T>
void
Mat44<T>::ortho(const Vec3<T>& min, const Vec3<T>& max)
{
    Vec3<T> tmp(max - min),
            inv(T(1)/tmp.x, T(1)/tmp.y, T(1)/tmp.z);
    tmp = min + max;
    
    set(2*inv.x,       0,        0, -tmp.x*inv.x,
              0, 2*inv.y,        0, -tmp.y*inv.y,
              0,       0, -2*inv.z, -tmp.z*inv.z,
              0,       0,        0,            1);
}

	
//! Inverse of ortho()
template <typename T>
void
Mat44<T>::invOrtho(const Vec3<T>& min, const Vec3<T>& max)
{
    Vec3<T> tmp1(max - min),
            tmp2(max + min);
    tmp1 *= T(0.5);
    tmp2 *= T(0.5);
    
    set(tmp1.x,      0,       0, -tmp2.x,
             0, tmp1.y,       0, -tmp2.y,
             0,      0, -tmp1.z,  tmp2.z,
             0,      0,       0,       1);
}



//----------------------------------
// Reading and writing
//----------------------------------

template <typename T>
inline std::ostream &
operator<<(std::ostream &s, const Mat44<T> & m)
{
#ifndef HAVE_IOS_BASE
    std::ios::fmtflags oldFlags = s.flags();
    long width;

    if (s.flags() & std::ios::fixed)
    {
        s.setf(std::ios::showpoint);
        width = s.precision() + 5;
    }
    else
    {
        s.setf(std::ios::scientific);
        s.setf(std::ios::showpoint);
        width = s.precision() + 8;
    }
#else
    std::ios_base::fmtflags oldFlags = s.flags();
    int width;

    if (s.flags() & std::ios_base::fixed)
    {
        s.setf(std::ios_base::showpoint);
        width = s.precision() + 5;
    }
    else
    {
        s.setf(std::ios_base::scientific);
        s.setf(std::ios_base::showpoint);
        width = s.precision() + 8;
    }
#endif

    s << " " << std::setw(width) << m[0][0] <<
         " " << std::setw(width) << m[1][0] <<
         " " << std::setw(width) << m[2][0] <<
         " " << std::setw(width) << m[3][0] << "\n" <<

         " " << std::setw(width) << m[0][1] <<
         " " << std::setw(width) << m[1][1] <<
         " " << std::setw(width) << m[2][1] <<
         " " << std::setw(width) << m[3][1] << "\n" <<

         " " << std::setw(width) << m[0][2] <<
         " " << std::setw(width) << m[1][2] <<
         " " << std::setw(width) << m[2][2] <<
         " " << std::setw(width) << m[3][2] << "\n" <<

         " " << std::setw(width) << m[0][3] <<
         " " << std::setw(width) << m[1][3] <<
         " " << std::setw(width) << m[2][3] <<
         " " << std::setw(width) << m[3][3] << "\n";

    s.flags(oldFlags);
    return s;
}


//! Input from a stream.
template <typename T>
inline std::istream &
operator>>(std::istream& in, Mat44<T>& m)
{
    return in >> m[0][0] >> m[1][0] >> m[2][0] >> m[3][0] >>
                 m[0][1] >> m[1][1] >> m[2][1] >> m[3][1] >>
                 m[0][2] >> m[1][2] >> m[2][2] >> m[3][2] >>
                 m[0][3] >> m[1][3] >> m[2][3] >> m[3][3];
}

} // namespace Math

#endif    // MATH_MAT44_H
