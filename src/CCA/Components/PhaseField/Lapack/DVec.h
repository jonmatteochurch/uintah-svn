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
 * @file CCA/Components/PhaseField/Lapack/DVec.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_Lapack_DVec_h
#define Packages_Uintah_CCA_Components_PhaseField_Lapack_DVec_h

#include <sci_defs/lapack_defs.h>

#ifdef HAVE_LAPACK

namespace Uintah
{
namespace PhaseField
{
namespace Lapack
{

class DMat;
class PMat;

class DVec
{
    double * raw;

public:
    int m;
    int n;

    DVec ( int m );
    DVec ( int m, double v );
    DVec ( int m, const double * v );
//     DVec ( std::initializer_list<double> l );
    DVec ( const DVec & copy );
    ~DVec();

    double * data() const;
    double & operator() ( int i );
    const double & operator() ( int i ) const;

    void premult_trsp_q ( const DMat & qrp );
    void premult_inv_r ( const DMat & qrp );
    void premult_inv_perm ( const PMat & p );
};

} // namespace Lapack
} // namespace PhaseField
} // namespace Uintah

#endif
#endif
