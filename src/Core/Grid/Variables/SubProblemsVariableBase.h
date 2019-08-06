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

#ifndef UINTAH_HOMEBREW_SubProblemsVariableBase_H
#define UINTAH_HOMEBREW_SubProblemsVariableBase_H

#include <CCA/Ports/ApplicationInterface.h>

#include <Core/Geometry/IntVector.h>
#include <Core/Grid/Variables/Variable.h>
#include <Core/Parallel/BufferInfo.h>
#include <Core/Exceptions/InternalError.h>

#include <string>

namespace Uintah
{

class Patch;
class RefCounted;

/**************************************

CLASS
 SubProblemsVariableBase

 Short description...

GENERAL INFORMATION

 SubProblemsVariableBase.h

 Jon Matteo Church
 School of Computing
 University of Leeds

KEYWORDS
 SubProblemsVariableBase

DESCRIPTION
 Long description...

WARNING

****************************************/

// inherits from Variable solely for the purpose of stuffing it in the DW
class SubProblemsVariableBase : public Variable
{
protected:
    SubProblemsVariableBase() = default;

private:
    SubProblemsVariableBase ( const SubProblemsVariableBase & ) = delete;
    SubProblemsVariableBase & operator= ( const SubProblemsVariableBase & ) = delete;

public:
    virtual ~SubProblemsVariableBase() = default;

public: // VARIABLE METHODS

    virtual const TypeDescription * virtualGetTypeDescription() const override = 0;

    virtual void emitNormal ( std::ostream & out, const IntVector & l, const IntVector & h, ProblemSpecP varnode, bool outputDoubleAsFloat ) override = 0;;

    virtual void readNormal ( std::istream & in, bool swapbytes ) override = 0;

    virtual void allocate ( const Patch * patch, const IntVector & boundary ) override;

    virtual void getSizeInfo ( std::string & elems, unsigned long & totsize, void *& ptr ) const override = 0;

    virtual size_t getDataSize() const override = 0;

    virtual bool copyOut ( void * dst ) const override = 0;

    virtual void copyPointer ( Variable & ) override = 0;

    virtual RefCounted * getRefCounted() override;

public: // SUBPROBLEMSVARIABLEBASE METHODS

    virtual void allocate ( const ApplicationInterface * ai, const VarLabel * label, int material, const Patch * patch ) = 0;

    virtual void addGhostRegion ( const IntVector & low, const IntVector & high ) = 0;

    virtual SubProblemsVariableBase * clone() const = 0;

    virtual void * getBasePointer() const = 0;

    virtual void getMPIBuffer ( BufferInfo & ) const = 0;
};

} // End namespace Uintah

#endif
