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


#ifndef Packages_Uintah_Core_Grid_Stencil7_h
#define Packages_Uintah_Core_Grid_Stencil7_h

#include <Core/Disclosure/TypeUtils.h>
#include <Core/Util/FancyAssert.h>
#include <iostream>

namespace Uintah {
  class TypeDescription;
  struct Stencil7 {
    // The order of this is designed to match the order of faces in Patch
    // Do not change it!
    //     -x +x -y +y -z +z
//  double  w, e, s, n, b, t;
    // diagonal term
//  double p;

// scjmc: if an instance is optimized out adjacency is not granted thus
//        operator[] may behave unexpectedly. We need an array to
//        enforce adjacency, using references to array elements for
//        previous members

    double m_value[7];
    double &w, &e, &s, &n, &b, &t;
    double &p;

    double& operator[](int index) {
      ASSERTRANGE(index, 0, 7);
      return m_value[index];
    }
    const double& operator[](int index) const {
      ASSERTRANGE(index, 0, 7);
      return m_value[index];
    }

    void initialize(double a){
      w = a;
      e = a;
      s = a;
      n = a; 
      b = a; 
      t = a;
      p = a;
    }

    // constructors
    Stencil7() :
      w(m_value[0]),
      e(m_value[1]),
      s(m_value[2]),
      n(m_value[3]),
      b(m_value[4]),
      t(m_value[5]),
      p(m_value[6])
    {}

    inline Stencil7(double init) : 
      Stencil7()
    {
      std::fill(m_value,m_value+7,init);
    }

    Stencil7& operator=(const Stencil7& that )
    {
      std::memcpy ( m_value, that.m_value, sizeof m_value );
      return *this;
    }
  };

  std::ostream & operator << (std::ostream &out, const Uintah::Stencil7 &a);

}

namespace Uintah {
   void swapbytes( Uintah::Stencil7& );
} // namespace Uintah

#endif
