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
 * @file CCA/Components/PhaseField/PostProcess/ArmPostProcessorFactory.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2020/06
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_PostProcess_ArmPostProcessorFactory_h
#define Packages_Uintah_CCA_Components_PhaseField_PostProcess_ArmPostProcessorFactory_h

#include <CCA/Components/PhaseField/PostProcess/ArmPostProcessorParallelPoly.h>
#include <CCA/Components/PhaseField/PostProcess/ArmPostProcessorDiagonalPoly.h>
#include <CCA/Components/PhaseField/PostProcess/ArmPostProcessorParallelTanh.h>
#include <CCA/Components/PhaseField/PostProcess/ArmPostProcessorDiagonalTanh.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/GridP.h>

namespace Uintah
{

namespace PhaseField
{

template<VarType VAR, DimType DIM>
class ArmPostProcessorFactory
{
public:
    static ArmPostProcessor<VAR, DIM> *
    create (
        ProblemSpecP spec,
        const double & epsilon,
        bool dbg
    )
    {
        if ( !spec )
            return nullptr;

        std::string type;
        if ( !spec->getAttribute ( "type", type ) )
            SCI_THROW ( ProblemSetupException ( "Cannot find type attribute in ArmPostProcessor block within problem specification file.", __FILE__, __LINE__ ) );
        std::transform ( type.begin(), type.end(), type.begin(), ::tolower );

        if ( type == "polynomial" )
        {
            int npts0 = 6;
            spec->getWithDefault ( "npts0", npts0, npts0 );

            int npts1 = 5;
            spec->getWithDefault ( "npts1", npts1, npts1 );

            int npts2 = npts1;
            spec->getWithDefault ( "npts2", npts2, npts2 );

            int npts3 = npts0;
            spec->getWithDefault ( "npts3", npts3, npts3 );

            int deg0 = npts0 - 1;
            spec->getWithDefault ( "deg0", deg0, deg0 );

            int deg1 = 1;
            spec->getWithDefault ( "deg1", deg1, deg1 );

            int deg2 = npts2 - 1;
            spec->getWithDefault ( "deg2", deg2, deg2 );

            int deg3 = npts3 - 1;
            spec->getWithDefault ( "deg3", deg3, deg3 );

            if ( deg0 >= npts0 )
                SCI_THROW ( ProblemSetupException ( "Cannot use polynomial interpolation with deg0 >= npts0.", __FILE__, __LINE__ ) );
            if ( deg1 >= npts1 )
                SCI_THROW ( ProblemSetupException ( "Cannot use polynomial interpolation with deg1 >= npts1.", __FILE__, __LINE__ ) );
            if ( deg2 >= npts2 )
                SCI_THROW ( ProblemSetupException ( "Cannot use polynomial interpolation with deg2 >= npts2.", __FILE__, __LINE__ ) );
            if ( deg3 >= npts3 )
                SCI_THROW ( ProblemSetupException ( "Cannot use polynomial interpolation with deg3 >= npts3.", __FILE__, __LINE__ ) );

            double alpha = 1.;
            spec->getWithDefault ( "alpha", alpha, alpha );

            if ( epsilon < 0 )
                return scinew ArmPostProcessorDiagonalPoly<VAR, DIM> ( npts0, npts1, npts2, npts3, deg0, deg1, deg2, deg3, alpha, dbg );
            else
                return scinew ArmPostProcessorParallelPoly<VAR, DIM> ( npts0, npts1, npts2, npts3, deg0, deg1, deg2, deg3, alpha, dbg );
        }
        if ( type == "tanh" )
        {
            int npts0 = 6;
            int npts1 = 4;
            int npts3 = 6;
            spec->getWithDefault ( "npts0", npts0, npts0 );
            spec->getWithDefault ( "npts1", npts1, npts1 );
            spec->getWithDefault ( "npts3", npts3, npts3 );

            double ftol = 1e-8, xtol = 1e-8, gtol = 1e-8, trtol = 1e-2;
            int max_nfev = 200, max_triter = 10;
            spec->getWithDefault ( "ftol", ftol, ftol );
            spec->getWithDefault ( "xtol", xtol, xtol );
            spec->getWithDefault ( "gtol", gtol, gtol );
            spec->getWithDefault ( "trtol", trtol, trtol );
            spec->getWithDefault ( "max_nfev", max_nfev, max_nfev );
            spec->getWithDefault ( "max_triter", max_triter, max_triter );

            double alpha = 1.;
            spec->getWithDefault ( "alpha", alpha, alpha );

            if ( epsilon < 0 )
                return scinew ArmPostProcessorDiagonalTanh<VAR, DIM> ( {ftol, xtol, gtol, trtol, max_nfev, max_triter}, npts0, npts1, npts3, alpha, dbg );
            else
                return scinew ArmPostProcessorParallelTanh<VAR, DIM> ( {ftol, xtol, gtol, trtol, max_nfev, max_triter}, npts0, npts1, npts3, alpha, dbg );
        }

        SCI_THROW ( ProblemSetupException ( "Cannot Create ArmPostProcessor of type '" + type + "'", __FILE__, __LINE__ ) );
    }
}; // class ArmPostProcessorFactory

} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_DataTypes_ArmPostProcessor_h



