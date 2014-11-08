/*! \file LimitsT.h
    \brief Contains the limitations of the basic C++ numerical data types
    \author Wojciech Jarosz
*/
#ifndef MATH_LIMITS_H
#define MATH_LIMITS_H

#include <float.h>
#include <limits.h>

namespace Math
{

//! Base template for numeric Limits
/*!
    This class is only provided because the standard numeric_limits
    class does not seem to provide a means of easily attaining the
    largest negative value for a type. If I find out this is incorrect,
    this class will disappear.
    
    Specializations for built-in types are provided.
*/
template <typename T>
struct Limits
{
    static T    min();
    static T    max();
    static T    smallest();
    static T    epsilon();
    static bool isIntegral();
    static bool isSigned();
};

//! Limits specialization for char
template <>
struct Limits <char>
{
    static char min()        {return CHAR_MIN;}
    static char max()        {return CHAR_MAX;}
    static char smallest()   {return 1;}
    static char epsilon()    {return 1;}
    static bool isIntegral() {return true;}
    static bool isSigned()   {return (char) ~0 < 0;}
};

//! Limits specialization for signed char
template <>
struct Limits <signed char>
{
    static signed char min()        {return SCHAR_MIN;}
    static signed char max()        {return SCHAR_MAX;}
    static signed char smallest()   {return 1;}
    static signed char epsilon()    {return 1;}
    static bool        isIntegral() {return true;}
    static bool        isSigned()   {return true;}
};

//! Limits specialization for unsigned char
template <>
struct Limits <unsigned char>
{
    static unsigned char min()        {return 0;}
    static unsigned char max()        {return UCHAR_MAX;}
    static unsigned char smallest()   {return 1;}
    static unsigned char epsilon()    {return 1;}
    static bool          isIntegral() {return true;}
    static bool          isSigned()   {return false;}
};

//! Limits specialization for short
template <>
struct Limits <short>
{
    static short min()        {return SHRT_MIN;}
    static short max()        {return SHRT_MAX;}
    static short smallest()   {return 1;}
    static short epsilon()    {return 1;}
    static bool  isIntegral() {return true;}
    static bool  isSigned()   {return true;}
};

//! Limits specialization for unsigned short
template <>
struct Limits <unsigned short>
{
    static unsigned short min()        {return 0;}
    static unsigned short max()        {return USHRT_MAX;}
    static unsigned short smallest()   {return 1;}
    static unsigned short epsilon()    {return 1;}
    static bool           isIntegral() {return true;}
    static bool           isSigned()   {return false;}
};

//! Limits specialization for int
template <>
struct Limits <int>
{
    static int  min()        {return INT_MIN;}
    static int  max()        {return INT_MAX;}
    static int  smallest()   {return 1;}
    static int  epsilon()    {return 1;}
    static bool isIntegral() {return true;}
    static bool isSigned()   {return true;}
};

//! Limits specialization for unsigned int
template <>
struct Limits <unsigned int>
{
    static unsigned int min()        {return 0;}
    static unsigned int max()        {return UINT_MAX;}
    static unsigned int smallest()   {return 1;}
    static unsigned int epsilon()    {return 1;}
    static bool         isIntegral() {return true;}
    static bool         isSigned()   {return false;}
};

//! Limits specialization for long
template <>
struct Limits <long>
{
    static long min()        {return LONG_MIN;}
    static long max()        {return LONG_MAX;}
    static long smallest()   {return 1;}
    static long epsilon()    {return 1;}
    static bool isIntegral() {return true;}
    static bool isSigned()   {return true;}
};

//! Limits specialization for unsigned long
template <>
struct Limits <unsigned long>
{
    static unsigned long min()        {return 0;}
    static unsigned long max()        {return ULONG_MAX;}
    static unsigned long smallest()   {return 1;}
    static unsigned long epsilon()    {return 1;}
    static bool          isIntegral() {return true;}
    static bool          isSigned()   {return false;}
};

//! Limits specialization for float
template <>
struct Limits <float>
{
    static float min()        {return -FLT_MAX;}
    static float max()        {return FLT_MAX;}
    static float smallest()   {return FLT_MIN;}
    static float epsilon()    {return FLT_EPSILON;}
    static bool  isIntegral() {return false;}
    static bool  isSigned()   {return true;}
};

//! Limits specialization for double
template <>
struct Limits <double>
{
    static double min()        {return -DBL_MAX;}
    static double max()        {return DBL_MAX;}
    static double smallest()   {return DBL_MIN;}
    static double epsilon()    {return DBL_EPSILON;}
    static bool   isIntegral() {return false;}
    static bool   isSigned()   {return true;}
};

//! Limits specialization for long double
template <>
struct Limits <long double>
{
    static long double min()        {return -LDBL_MAX;}
    static long double max()        {return LDBL_MAX;}
    static long double smallest()   {return LDBL_MIN;}
    static long double epsilon()    {return LDBL_EPSILON;}
    static bool        isIntegral() {return false;}
    static bool        isSigned()   {return true;}
};

} // namespace Math

#endif // MATH_LIMITS_H
