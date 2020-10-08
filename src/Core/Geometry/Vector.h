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
 *  Vector.h: ?
 *
 *  Written by:
 *   Author?
 *   Department of Computer Science
 *   University of Utah
 *   Date?
 *
 */

#ifndef Geometry_Vector_h
#define Geometry_Vector_h 1


#include <Core/Util/Assert.h>
#include <Core/Math/Expon.h>

#include   <string>
#include   <iosfwd>
#include   <cstring>

namespace Uintah {

class Point;
class TypeDescription;

class Vector {
  double m_value[3];

public:
  inline Vector() {}
  inline Vector(double x, double y, double z) : m_value{x, y, z} {}
  inline Vector(const Vector& v) { std::memcpy(m_value, v.m_value, sizeof m_value); }
  inline explicit Vector(double init) : m_value{init, init, init} {}
  inline explicit Vector(const Point& p);
  inline double length() const;
  inline double length2() const;
  friend inline double Dot(const Vector&, const Vector&);
  friend inline double Dot(const Point&, const Vector&);
  friend inline double Dot(const Vector&, const Point&);
  inline Vector& operator=(const Vector&);
  inline Vector& operator=(const double&);
  inline Vector& operator=(const int&);
  
  static Vector fromString( const std::string & source ); // Creates a Vector from a string that looksl like "[Num, Num, Num]".

#ifdef COMMENT_OUT
  /* !!!
  () index from 0
  [] index from 0
  !!! */

  //Note vector(0)=vector.x();vector(1)=vector.y();vector(2)=vector.z()
  inline double& operator()(int idx) {
    // Ugly, but works
    return m_value[idx];
  }

  //Note vector(0)=vector.x();vector(1)=vector.y();vector(2)=vector.z()
  inline double operator()(int idx) const {
    // Ugly, but works
    return m_value[idx];
  }
#endif

  //Note vector[0]=vector.x();vector[1]=vector.y();vector[2]=vector.z()
  inline double& operator[](int idx) {
    // Ugly, but works
    return m_value[idx];
  }

  //Note vector[0]=vector.x();vector[1]=vector.y();vector[2]=vector.z()
  inline double operator[](int idx) const {
    // Ugly, but works
    return m_value[idx];
  }

  // checks if one vector is exactly the same as another
  int operator==(const Vector&) const;
  int operator!=(const Vector&) const;

  inline Vector operator*(const double) const;
  inline Vector operator*(const Vector&) const;
  inline Vector& operator*=(const double);
  inline Vector& operator*=(const Vector&);
  inline Vector operator/(const double) const;
  inline Vector operator/(const Vector&) const;
  inline Vector& operator/=(const double);
  inline Vector operator+(const Vector&) const;
  inline Vector& operator+=(const Vector&);
  inline Vector operator-() const;
  inline Vector operator-(const Vector&) const;
  inline Vector operator-(const Point&) const;
  inline Vector& operator-=(const Vector&);
  inline double normalize();
  inline double safe_normalize();
  Vector normal() const;
  friend inline Vector Cross(const Vector&, const Vector&);
  friend inline Vector Abs(const Vector&);
  inline void x(double);
  inline double x() const;
  inline void y(double);
  inline double y() const;
  inline void z(double);
  inline double z() const;

  inline void u(double);
  inline double u() const;
  inline void v(double);
  inline double v() const;
  inline void w(double);
  inline double w() const;

  void rotz90(const int);
    
  std::string get_string() const;

  //! support dynamic compilation
  static const std::string& get_h_file_path();

  friend class Point;
    
  friend inline Vector Interpolate(const Vector&, const Vector&, double);
    
  void find_orthogonal(Vector&, Vector&) const;
  bool check_find_orthogonal(Vector&, Vector&) const;


  inline const Point &point() const;
  inline Point &asPoint() const;
  inline double minComponent() const {
    if(m_value[0]<m_value[1]){
      if(m_value[0]<m_value[2])
        return m_value[0];
      else
        return m_value[2];
    } else {
      if(m_value[1]<m_value[2])
        return m_value[1];
      else
        return m_value[2];
    }
  }
  inline double maxComponent() const {
    if(m_value[0]>m_value[1]){
      if(m_value[0]>m_value[2])
        return m_value[0];
      else
        return m_value[2];
    } else {
      if(m_value[1]>m_value[2])
        return m_value[1];
      else
        return m_value[2];
    }
  }

  inline void Set(double x, double y, double z)
    { 
      m_value[0] = x;
      m_value[1] = y;
      m_value[2] = z;
    }
      
  friend std::ostream& operator<<(std::ostream& os, const Vector& p);
  friend std::istream& operator>>(std::istream& os, Vector& p);

}; // end class Vector

// Actual declarations of these functions as 'friend' above doesn't
// (depending on the compiler) actually declare them.

std::ostream& operator<<(std::ostream& os, const Vector& p);
std::istream& operator>>(std::istream& os, Vector& p);
  
} // End namespace Uintah

// This cannot be above due to circular dependencies
#include <Core/Geometry/Point.h>

namespace Uintah {

inline Vector::Vector(const Point& p)
: Vector( p.asVector() )
{}

inline double Vector::length2() const
{
    return m_value[0]*m_value[0]+m_value[1]*m_value[1]+m_value[2]*m_value[2];
}

inline Vector& Vector::operator=(const Vector& v)
{
    std::memcpy(m_value, v.m_value, sizeof m_value);
    return *this;
}

// for initializing in dynamic code
// one often want template<class T> T val = 0.0;

inline Vector& Vector::operator=(const double& d)
{
  m_value[0] = d;
  m_value[1] = d;
  m_value[2] = d;
  return *this;
}

inline Vector& Vector::operator=(const int& d)
{
  m_value[0] = static_cast<int>(d);
  m_value[1] = static_cast<int>(d);
  m_value[2] = static_cast<int>(d);
  return *this;
}

inline bool operator<(Vector v1, Vector v2)
{
  return(v1.length()<v2.length());
}

inline bool operator<=(Vector v1, Vector v2)
{
  return(v1.length()<=v2.length());
}

inline bool operator>(Vector v1, Vector v2)
{
  return(v1.length()>v2.length());
}

inline bool operator>=(Vector v1, Vector v2)
{
  return(v1.length()>=v2.length());
}



inline Vector Vector::operator*(const double s) const
{
    return Vector(m_value[0]*s, m_value[1]*s, m_value[2]*s);
}

inline Vector& Vector::operator*=(const Vector& v)
{
  m_value[0] *= v.m_value[0];
  m_value[1] *= v.m_value[1];
  m_value[2] *= v.m_value[2];
  return *this;
}

// Allows for double * Vector so that everything doesn't have to be
// Vector * double
inline Vector operator*(const double s, const Vector& v) {
    return v*s;
}

inline Vector Vector::operator/(const double d) const
{
    return Vector(m_value[0]/d, m_value[1]/d, m_value[2]/d);
}

inline Vector Vector::operator/(const Vector& v2) const
{
    return Vector(m_value[0]/v2.m_value[0], m_value[1]/v2.m_value[1], m_value[2]/v2.m_value[2]);
}

inline Vector Vector::operator+(const Vector& v2) const
{
    return Vector(m_value[0]+v2.m_value[0], m_value[1]+v2.m_value[1], m_value[2]+v2.m_value[2]);
}

inline Vector Vector::operator*(const Vector& v2) const
{
    return Vector(m_value[0]*v2.m_value[0], m_value[1]*v2.m_value[1], m_value[2]*v2.m_value[2]);
}

inline Vector Vector::operator-(const Vector& v2) const
{
    return Vector(m_value[0]-v2.m_value[0], m_value[1]-v2.m_value[1], m_value[2]-v2.m_value[2]);
}

inline Vector Vector::operator-(const Point& v2) const
{
    return Vector(m_value[0]-v2.m_value[0], m_value[1]-v2.m_value[1], m_value[2]-v2.m_value[2]);
}

inline Vector& Vector::operator+=(const Vector& v2)
{
    m_value[0]+=v2.m_value[0];
    m_value[1]+=v2.m_value[1];
    m_value[2]+=v2.m_value[2];
    return *this;
}

inline Vector& Vector::operator-=(const Vector& v2)
{
    m_value[0]-=v2.m_value[0];
    m_value[1]-=v2.m_value[1];
    m_value[2]-=v2.m_value[2];
    return *this;
}

inline Vector Vector::operator-() const
{
    return Vector(-m_value[0],-m_value[1],-m_value[2]);
}

inline double Vector::length() const
{
    return Sqrt(m_value[0]*m_value[0]+m_value[1]*m_value[1]+m_value[2]*m_value[2]);
}

inline Vector Abs(const Vector& v)
{
    double x=v.m_value[0]<0?-v.m_value[0]:v.m_value[0];
    double y=v.m_value[1]<0?-v.m_value[1]:v.m_value[1];
    double z=v.m_value[2]<0?-v.m_value[2]:v.m_value[2];
    return Vector(x,y,z);
}

inline Vector Cross(const Vector& v1, const Vector& v2)
{
    return Vector(
        v1.m_value[1]*v2.m_value[2]-v1.m_value[2]*v2.m_value[1],
        v1.m_value[2]*v2.m_value[0]-v1.m_value[0]*v2.m_value[2],
        v1.m_value[0]*v2.m_value[1]-v1.m_value[1]*v2.m_value[0]);
}

inline Vector Interpolate(const Vector& v1, const Vector& v2,
                          double weight)
{
    double weight1=1.0-weight;
    return Vector(
        v2.m_value[0]*weight+v1.m_value[0]*weight1,
        v2.m_value[1]*weight+v1.m_value[1]*weight1,
        v2.m_value[2]*weight+v1.m_value[2]*weight1);
}

inline Vector& Vector::operator*=(const double d)
{
    m_value[0]*=d;
    m_value[1]*=d;
    m_value[2]*=d;
    return *this;
}

inline Vector& Vector::operator/=(const double d)
{
    m_value[0]/=d;
    m_value[1]/=d;
    m_value[2]/=d;
    return *this;
}

inline void Vector::x(double d)
{
    m_value[0]=d;
}

inline double Vector::x() const
{
    return m_value[0];
}

inline void Vector::y(double d)
{
    m_value[1]=d;
}

inline double Vector::y() const
{
    return m_value[1];
}

inline void Vector::z(double d)
{
    m_value[2]=d;
}

inline double Vector::z() const
{
    return m_value[2];
}



inline void Vector::u(double d)
{
    m_value[0]=d;
}

inline double Vector::u() const
{
    return m_value[0];
}

inline void Vector::v(double d)
{
    m_value[1]=d;
}

inline double Vector::v() const
{
    return m_value[1];
}

inline void Vector::w(double d)
{
    m_value[2]=d;
}

inline double Vector::w() const
{
    return m_value[2];
}

inline double Dot(const Vector& v1, const Vector& v2)
{
    return v1.m_value[0]*v2.m_value[0]+v1.m_value[1]*v2.m_value[1]+v1.m_value[2]*v2.m_value[2];
}

inline double Dot(const Vector& v, const Point& p)
{
    return v.m_value[0]*p.m_value[0]+v.m_value[1]*p.m_value[1]+v.m_value[2]*p.m_value[2];
}

inline
double Vector::normalize()
{
    double l2=m_value[0]*m_value[0]+m_value[1]*m_value[1]+m_value[2]*m_value[2];
    double l=Sqrt(l2);
    ASSERT(l>0.0);
    m_value[0]/=l;
    m_value[1]/=l;
    m_value[2]/=l;
    return l;
}

inline
double Vector::safe_normalize()
{
  double l=Sqrt(m_value[0]*m_value[0] + m_value[1]*m_value[1] + m_value[2]*m_value[2] + 1.0e-12);
  m_value[0]/=l;
  m_value[1]/=l;
  m_value[2]/=l;
  return l;
}


inline const Point &Vector::point() const {
    return (const Point &)(*this);
}

inline Point &Vector::asPoint() const {
    static_assert(std::is_standard_layout<Point>::value, "Vector and Point are pointer-incovertible only if both standard-layout classes");
    return (Point &)(*this);
}


inline Vector Min(const Vector &v1, const Vector &v2)
{
  return Vector(Min(v1.x(), v2.x()),
                Min(v1.y(), v2.y()),
                Min(v1.z(), v2.z()));
}

inline Vector Max(const Vector &v1, const Vector &v2)
{
  return Vector(Max(v1.x(), v2.x()),
                Max(v1.y(), v2.y()),
                Max(v1.z(), v2.z()));
}

const TypeDescription* get_type_description(Vector*);

} // End namespace Uintah


#endif
