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
 * @file CCA/Components/PhaseField/DataTypes/HeatTestProblem.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_DataTypes_HeatTestProblem_h
#define Packages_Uintah_CCA_Components_PhaseField_DataTypes_HeatTestProblem_h

#include <CCA/Components/PhaseField/DataTypes/Problem.h>
#include <CCA/Components/PhaseField/DataTypes/ScalarField.h>
#include <CCA/Components/PhaseField/DataTypes/VectorField.h>
#include <CCA/Components/PhaseField/DataTypes/ScalarProblem.h>

namespace Uintah
{
namespace PhaseField
{

/// Type of Problem used by Heat application
#ifdef PhaseField_Heat_DBG_DERIVATIVES
template<VarType VAR, StnType STN> using HeatTestProblem = Problem < VAR, STN, ScalarField<const double>, VectorField<const double, get_stn<STN>::dim >, VectorField<const double, get_stn<STN>::dim > >;
#else
template<VarType VAR, StnType STN> using HeatTestProblem = Problem < VAR, STN, ScalarField<const double> >;
#endif

} // namespace PhaseFieldSTNSTN
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_DataTypes_HeatTestProblem_h



