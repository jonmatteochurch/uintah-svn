/*
 * The MIT License
 *
 * Copyright (c) 1997-2020 The University of Utah
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

/**
 * @file CCA/Components/PhaseField/Lapack/Poly.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2020/04
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_Lapack_Poly_h
#define Packages_Uintah_CCA_Components_PhaseField_Lapack_Poly_h

#include <sci_defs/lapack_defs.h>

#include<string>

#ifdef HAVE_LAPACK

namespace Uintah
{
namespace PhaseField
{
namespace Lapack
{

class Poly
{
    double * raw;
    double * rr;
    double * ri;

public:
    int d;
    int nr;

    Poly ( int d );
    Poly ( const Poly & copy );
    ~Poly();

    double & cfx ( int i );
    double cfx ( int i ) const;
    bool root_in_range ( double min, double max, double & r ) const;
    double operator() ( double x0 ) const;
    double dx ( double x0 ) const;
    double ddx ( double x0 ) const;

    void fit ( int n, const double * x, const double * y );
    void roots();

    Poly ( std::initializer_list<double> l ) :
        raw ( new double [l.size()] ),
        rr ( new double [l.size() - 1] ),
        ri ( new double [l.size() - 1] ),
        d ( l.size() - 1 ),
        nr ( -1 )
    {
        std::copy ( l.begin(), l.end(), raw );
    }
    std::string rts ( int i ) const
    {
        return std::to_string ( rr[i] ) + " + " + std::to_string ( ri[i] ) + "i";
    }
};

} // namespace Lapack
} // namespace PhaseField
} // namespace Uintah

#endif
#endif
