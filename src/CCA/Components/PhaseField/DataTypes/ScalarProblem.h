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

/**
 * @file CCA/Components/PhaseField/DataTypes/HeatProblem.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_DataTypes_ScalarProblem_h
#define Packages_Uintah_CCA_Components_PhaseField_DataTypes_ScalarProblem_h

#include <CCA/Components/PhaseField/DataTypes/Problem.h>
#include <CCA/Components/PhaseField/DataTypes/ScalarField.h>

namespace Uintah
{
namespace PhaseField
{

/// Type of Problem used by Heat and Bencmark01 application
template<VarType VAR, StnType STN> using ScalarProblem = Problem < VAR, STN, ScalarField<const double> >;

} // namespace PhaseFieldSTNSTN
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_DataTypes_ScalarProblem_h



