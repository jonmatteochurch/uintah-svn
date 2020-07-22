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
 *  Point.cc: ?
 *
 *  Written by:
 *   Author ?
 *   Department of Computer Science
 *   University of Utah
 *   Date ?
 *
 */

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Util/Assert.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/MiscMath.h>
#include <iostream>
#include <cstdio>

using namespace std;

namespace Uintah {


Point Interpolate(const Point& p1, const Point& p2, double w)
{
    return Point(
        Interpolate(p1.m_value[0], p2.m_value[0], w),
        Interpolate(p1.m_value[1], p2.m_value[1], w),
        Interpolate(p1.m_value[2], p2.m_value[2], w));
}

string Point::get_string() const
{
    char buf[100];
    sprintf(buf, "[%g, %g, %g]", m_value[0], m_value[1], m_value[2]);
    return buf;
}

int Point::operator==(const Point& p) const
{
    return p.m_value[0] == m_value[0] && p.m_value[1] == m_value[1] && p.m_value[2] == m_value[2];
}

int Point::operator!=(const Point& p) const
{
    return p.m_value[0] != m_value[0] || p.m_value[1] != m_value[1] || p.m_value[2] != m_value[2];
}

Point::Point(double x, double y, double z, double w)
: Point()
{
    if(w==0){
        cerr << "degenerate point!" << endl;
        m_value[0]=m_value[1]=m_value[2]=0;
    } else {
        m_value[0]=x/w;
        m_value[1]=y/w;
        m_value[2]=z/w;
    }
}

Point AffineCombination(const Point& p1, double w1,
                        const Point& p2, double w2)
{
    return Point(p1.m_value[0]*w1+p2.m_value[0]*w2,
                 p1.m_value[1]*w1+p2.m_value[1]*w2,
                 p1.m_value[2]*w1+p2.m_value[2]*w2);
}

Point AffineCombination(const Point& p1, double w1,
                        const Point& p2, double w2,
                        const Point& p3, double w3)
{
    return Point(p1.m_value[0]*w1+p2.m_value[0]*w2+p3.m_value[0]*w3,
                 p1.m_value[1]*w1+p2.m_value[1]*w2+p3.m_value[1]*w3,
                 p1.m_value[2]*w1+p2.m_value[2]*w2+p3.m_value[2]*w3);
}

Point AffineCombination(const Point& p1, double w1,
                        const Point& p2, double w2,
                        const Point& p3, double w3,
                        const Point& p4, double w4)
{
    return Point(p1.m_value[0]*w1+p2.m_value[0]*w2+p3.m_value[0]*w3+p4.m_value[0]*w4,
                 p1.m_value[1]*w1+p2.m_value[1]*w2+p3.m_value[1]*w3+p4.m_value[1]*w4,
                 p1.m_value[2]*w1+p2.m_value[2]*w2+p3.m_value[2]*w3+p4.m_value[2]*w4);
}

ostream& operator<<( ostream& os, const Point& p )
{
  os << '[' << p.x() << ' ' << p.y() << ' ' << p.z() << ']';
  return os;
}

istream& operator>>( istream& is, Point& v)
{
    double x, y, z;
    char st;
  is >> st >> x >> st >> y >> st >> z >> st;
  v=Point(x,y,z);
  return is;
}

int
Point::Overlap( double a, double b, double e )
{
  double hi, lo, h, l;
  
  hi = a + e;
  lo = a - e;
  h  = b + e;
  l  = b - e;

  if ( ( hi > l ) && ( lo < h ) )
    return 1;
  else
    return 0;
}
  
int
Point::InInterval( Point a, double epsilon )
{
  if ( Overlap( m_value[0], a.x(), epsilon ) &&
      Overlap( m_value[1], a.y(), epsilon )  &&
      Overlap( m_value[2], a.z(), epsilon ) )
    return 1;
  else
    return 0;
}




} // End namespace Uintah

