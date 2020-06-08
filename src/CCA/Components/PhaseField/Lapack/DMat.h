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
 * @file CCA/Components/PhaseField/Lapack/DMat.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_Lapack_DMat_h
#define Packages_Uintah_CCA_Components_PhaseField_Lapack_DMat_h

#include <sci_defs/lapack_defs.h>

#ifdef HAVE_LAPACK

#include <CCA/Components/PhaseField/Lapack/DVec.h>
#include <CCA/Components/PhaseField/Lapack/PMat.h>

namespace Uintah
{
namespace PhaseField
{
namespace Lapack
{

class DMat
{
    double * raw;

    PMat * pp;
    DVec * ptau;

    DVec * ps;
    DMat * pu;
    DMat * pvt;

public:
    int m;
    int n;
    int k;

    DMat ( int m, int n );
    DMat ( int m, int n, double v );
    DMat ( const DMat &copy );
    ~DMat();

    static DMat Vander ( int m, int n, const double * c );

    double * data () const;
    double * cdata ( int j ) const;

    double & operator() ( int i, int j );
    const double & operator() ( int i, int j ) const;

    const PMat & p() const;
    const DVec & tau() const;

    const DVec & s() const;
    const DMat & u() const;
    const DMat & vt() const;

    void qr ( bool preserve = false );
    void svd ( bool preserve = false );
};

} // namespace Lapack
} // namespace PhaseField
} // namespace Uintah

#endif
#endif
