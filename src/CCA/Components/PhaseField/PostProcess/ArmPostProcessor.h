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
 * @file CCA/Components/PhaseField/PostProcess/ArmPostProcessor.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2020/04
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_PostProcess_ArmPostProcessor_h
#define Packages_Uintah_CCA_Components_PhaseField_PostProcess_ArmPostProcessor_h

#include <CCA/Components/PhaseField/Util/Definitions.h>
#include <Core/Exceptions/TypeMismatchException.h>

namespace Uintah
{

namespace PhaseField
{

class ArmPostProcessor
{
public:
    virtual ~ArmPostProcessor() = default;
    virtual void setLevel ( const Level * level ) = 0;
    virtual void initializeLocations () = 0;
    virtual void setLocations ( const Patch * patch, IntVector const & low, IntVector const & high, const std::list<Patch::FaceType> & faces, View < ScalarField<const double> > const & psi ) = 0;
    virtual void reduceLocations ( const ProcessorGroup * myworld ) = 0;
    virtual void printLocations ( std::ostream & out ) const = 0;
    virtual void initializeData () = 0;
    virtual void setData ( IntVector const & low, IntVector const & high, View < ScalarField<const double> > const & psi, View< ScalarField<int> > * refine_flag ) = 0;
    virtual void reduceData ( const ProcessorGroup * myworld ) = 0;
    virtual void printData ( std::ostream & out ) const = 0;
    virtual void computeTipInfo ( double & tip_position, double tip_curvatures[3] ) = 0;
};

} // namespace PhaseField
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_PostProcess_ArmPostProcessor_h



