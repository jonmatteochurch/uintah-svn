/*
 * The MIT License
 *
 * Copyright (c) 1997-2019 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */


/*
 *  Point.h: ?
 *
 *  Written by:
 *   Author?
 *   Department of Computer Science
 *   University of Utah
 *   Date?
 *
 */

#ifndef Geometry_Point_h
#define Geometry_Point_h 1

#include <Core/Math/MinMax.h>

#include <float.h>

#include   <string>
#include   <iosfwd>
#include   <cstring>

namespace Uintah {
    
class TypeDescription;
class RigorousTest;
class Vector;

class Point {
  double m_value[3];

public:
  inline Point() : m_value{DBL_MAX, DBL_MAX, DBL_MAX} {}
  inline Point(double x, double y, double z) : m_value{x, y, z} {}
  inline Point(const Point& p) { std::memcpy(m_value, p.m_value, sizeof m_value); }
  inline explicit Point(const Vector& v);
  Point(double, double, double, double);
  int operator==(const Point&) const;
  int operator!=(const Point&) const;
  inline Point& operator=(const Point&);
  inline Vector operator+(const Point&) const;
  inline Vector operator-(const Point&) const;
  inline Point operator+(const Vector&) const;
  inline Point operator-(const Vector&) const;
  inline Point operator*(double) const;
  inline Point& operator*=(const double);
  inline Point& operator+=(const Vector&);
  inline Point& operator-=(const Vector&);
  inline Point& operator/=(const double);
  inline Point operator/(const double) const;
  inline Point operator-() const;
  inline double& operator()(int idx);
  inline double operator()(int idx) const;

  inline void addscaled(const Point& p, const double scale);  // this += p * w;

  inline double x() const;
  inline double y() const;
  inline double z() const;

  inline void x( const double );
  inline void y( const double );
  inline void z( const double );

  inline Vector &asVector() const;
  inline Vector toVector() const;

  std::string get_string() const;

  //! support dynamic compilation
  static const std::string& get_h_file_path();
    
  friend class Vector;
  friend inline double Dot(const Point&, const Point&);
  friend inline double Dot(const Vector&, const Point&);
  friend inline double Dot(const Point&, const Vector&);
  //    friend inline double Dot(const Point&, const Vector&);
  friend inline Point Min(const Point&, const Point&);
  friend inline Point Max(const Point&, const Point&);
  friend Point Interpolate(const Point&, const Point&, double);
  friend Point AffineCombination(const Point&, double,
                                              const Point&, double,
                                              const Point&, double,
                                              const Point&, double);
  friend Point AffineCombination(const Point&, double,
                                              const Point&, double,
                                              const Point&, double);
  friend Point AffineCombination(const Point&, double,
                                              const Point&, double);



  // is one point within a small interval of another?

  int Overlap( double a, double b, double e );
  int InInterval( Point a, double epsilon );
    
  static void test_rigorous(RigorousTest* __test);

  friend std::ostream& operator<<(std::ostream& os, const Point& p);
  friend std::istream& operator>>(std::istream& os, Point& p);

}; // end class Point


// Actual declarations of these functions as 'friend' above doesn't
// (depending on the compiler) actually declare them.
Point Interpolate(const Point&, const Point&, double);
Point AffineCombination(const Point&, double, const Point&, double,
                                 const Point&, double, const Point&, double);
Point AffineCombination(const Point&, double, const Point&, double, const Point&, double);
Point AffineCombination(const Point&, double, const Point&, double);

std::ostream& operator<<(std::ostream& os, const Point& p);
std::istream& operator>>(std::istream& os, Point& p);

inline 
Point operator*(double d, const Point &p) {
  return p*d;
}
inline 
Point operator+(const Vector &v, const Point &p) {
  return p+v;
}

} // End namespace Uintah

// This cannot be above due to circular dependencies
#include <Core/Geometry/Vector.h>

namespace Uintah {

inline Point::Point(const Vector& v)
 : Point ( v.asPoint() )
{}

inline Point& Point::operator=(const Point& p)
{
    std::memcpy(m_value, p.m_value, sizeof m_value);
    return *this;
}

inline Vector Point::operator+(const Point& p) const
{
    return Vector(m_value[0]+p.m_value[0], m_value[1]+p.m_value[1], m_value[2]+p.m_value[2]);
}

inline Vector Point::operator-(const Point& p) const
{
    return Vector(m_value[0]-p.m_value[0], m_value[1]-p.m_value[1], m_value[2]-p.m_value[2]);
}

inline Point Point::operator+(const Vector& v) const
{
    return Point(m_value[0]+v.m_value[0], m_value[1]+v.m_value[1], m_value[2]+v.m_value[2]);
}

inline Point Point::operator-(const Vector& v) const
{
    return Point(m_value[0]-v.m_value[0], m_value[1]-v.m_value[1], m_value[2]-v.m_value[2]);
}

inline Point& Point::operator+=(const Vector& v)
{
    m_value[0]+=v.m_value[0];
    m_value[1]+=v.m_value[1];
    m_value[2]+=v.m_value[2];
    return *this;
}

inline Point& Point::operator-=(const Vector& v)
{
    m_value[0]-=v.m_value[0];
    m_value[1]-=v.m_value[1];
    m_value[2]-=v.m_value[2];
    return *this;
}

inline Point& Point::operator*=(const double d)
{
    m_value[0]*=d;
    m_value[1]*=d;
    m_value[2]*=d;
    return *this;
}

inline Point& Point::operator/=(const double d)
{
    m_value[0]/=d;
    m_value[1]/=d;
    m_value[2]/=d;
    return *this;
}

inline Point Point::operator-() const
{
    return Point(-m_value[0], -m_value[1], -m_value[2]);
}

inline Point Point::operator*(double d) const
{
    return Point(m_value[0]*d, m_value[1]*d, m_value[2]*d);
}

inline Point Point::operator/(const double d) const
{
    return Point(m_value[0]/d,m_value[1]/d,m_value[2]/d);
}

inline double& Point::operator()(int idx) {
        return (&m_value[0])[idx];
}

inline double Point::operator()(int idx) const {
        return (&m_value[0])[idx];
}

inline void Point::addscaled(const Point& p, const double scale) {
  // this += p * w;
  m_value[0] += p.m_value[0] * scale;
  m_value[1] += p.m_value[1] * scale;
  m_value[2] += p.m_value[2] * scale;
}

inline void Point::x(const double d)
{
    m_value[0]=d;
}

inline double Point::x() const
{
    return m_value[0];
}

inline void Point::y(const double d)
{
    m_value[1]=d;
}

inline double Point::y() const
{
    return m_value[1];
}

inline void Point::z(const double d)
{
    m_value[2]=d;
}

inline double Point::z() const
{
    return m_value[2];
}

inline Vector &Point::asVector() const
{
    static_assert(std::is_standard_layout<Vector>::value, "Vector and Point are pointer-incovertible only if both standard-layout classes");
    return (Vector &)(*this);
}

inline Vector Point::toVector() const
{
  return Vector(m_value[0],m_value[1],m_value[2]);
}

inline Point Min(const Point& p1, const Point& p2)
{

  double x=Min(p1.m_value[0], p2.m_value[0]);
  double y=Min(p1.m_value[1], p2.m_value[1]);
  double z=Min(p1.m_value[2], p2.m_value[2]);
  return Point(x,y,z);
}

inline Point Max(const Point& p1, const Point& p2)
{

  double x=Max(p1.m_value[0], p2.m_value[0]);
  double y=Max(p1.m_value[1], p2.m_value[1]);
  double z=Max(p1.m_value[2], p2.m_value[2]);
  return Point(x,y,z);
}

inline double Dot(const Point& p, const Vector& v)
{
    return p.m_value[0]*v.m_value[0]+p.m_value[1]*v.m_value[1]+p.m_value[2]*v.m_value[2];
}

inline double Dot(const Point& p1, const Point& p2)
{
  return p1.m_value[0]*p2.m_value[0] + p1.m_value[1]*p2.m_value[1] + p1.m_value[2]*p2.m_value[2];
}

const TypeDescription* get_type_description(Point*);

} // End namespace Uintah



#endif //ifndef Geometry_Point_h
