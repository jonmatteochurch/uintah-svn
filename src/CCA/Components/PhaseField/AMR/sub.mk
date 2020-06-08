#
#  The MIT License
#
#  Copyright (c) 1997-2020 The University of Utah
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to
#  deal in the Software without restriction, including without limitation the
#  rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
#  sell copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
#  IN THE SOFTWARE.
#
#
#
#
#
# Makefile fragment for this subdirectory

SRCDIR := CCA/Components/PhaseField/AMR

SRCS += \
  $(SRCDIR)/AMRFDViewScalarProblemCCP5FC0-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemCCP5FC1-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemCCP5FCSimple-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemCCP5FCLinear-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemCCP5FCBilinear-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemCCP7FC0-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemCCP7FC1-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemNCP5FC0-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemNCP5FC1-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemNCP5FCSimple-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemNCP5FCLinear-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemNCP5FCBilinear-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemNCP7FC0-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemNCP7FC1-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCP5FC0-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCP5FC1-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCP5FCSimple-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCP5FCLinear-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCP5FCBilinear-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCP7FC0-0-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCP7FC0-1-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCP7FC0-2-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCP7FC1-0-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCP7FC1-1-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCP7FC1-2-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCP5FC0-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCP5FC1-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCP5FCSimple-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCP5FCLinear-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCP5FCBilinear-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCP7FC0-0-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCP7FC0-1-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCP7FC0-2-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCP7FC1-0-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCP7FC1-1-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCP7FC1-2-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemCCP5FC0-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemCCP5FC1-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemCCP5FCSimple-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemCCP5FCLinear-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemCCP7FC0-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemCCP7FC1-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemCCP5FCBilinear-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemNCP5FC0-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemNCP5FC1-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemNCP5FCSimple-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemNCP5FCLinear-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemNCP5FCBilinear-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemNCP7FC0-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemNCP7FC1-bld.cc \

BLDDIR := $(SRCTOP)/$(SRCDIR)

BLDSRCS += \
  $(BLDDIR)/AMRFDViewScalarProblemCCP5FC0-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemCCP5FC1-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemCCP5FCSimple-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemCCP5FCLinear-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemCCP5FCBilinear-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemCCP7FC0-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemCCP7FC1-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemNCP5FC0-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemNCP5FC1-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemNCP5FCSimple-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemNCP5FCLinear-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemNCP5FCBilinear-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemNCP7FC0-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemNCP7FC1-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCP5FC0-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCP5FC1-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCP5FCSimple-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCP5FCLinear-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCP5FCBilinear-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCP7FC0-0-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCP7FC0-1-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCP7FC0-2-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCP7FC1-0-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCP7FC1-1-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCP7FC1-2-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCP5FC0-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCP5FC1-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCP5FCSimple-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCP5FCLinear-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCP5FCBilinear-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCP7FC0-0-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCP7FC0-1-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCP7FC0-2-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCP7FC1-0-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCP7FC1-1-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCP7FC1-2-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemCCP5FC0-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemCCP5FC1-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemCCP5FCSimple-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemCCP5FCLinear-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemCCP7FC0-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemCCP7FC1-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemCCP5FCBilinear-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemNCP5FC0-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemNCP5FC1-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemNCP5FCSimple-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemNCP5FCLinear-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemNCP5FCBilinear-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemNCP7FC0-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemNCP7FC1-bld.cc \

ifeq ($(HAVE_HYPRE),yes)
  SRCS += \
    $(SRCDIR)/AMRFDViewScalarProblemCCP5FC0New-bld.cc \
    $(SRCDIR)/AMRFDViewScalarProblemCCP5FC1New-bld.cc \
    $(SRCDIR)/AMRFDViewScalarProblemCCP7FC0New-bld.cc \
    $(SRCDIR)/AMRFDViewScalarProblemCCP7FC1New-bld.cc \
    $(SRCDIR)/AMRFDViewHeatTestProblemCCP5FC0New-bld.cc \
    $(SRCDIR)/AMRFDViewHeatTestProblemCCP5FC1New-bld.cc \
    $(SRCDIR)/AMRFDViewHeatTestProblemCCP7FC0New-0-bld.cc \
    $(SRCDIR)/AMRFDViewHeatTestProblemCCP7FC0New-1-bld.cc \
    $(SRCDIR)/AMRFDViewHeatTestProblemCCP7FC0New-2-bld.cc \
    $(SRCDIR)/AMRFDViewHeatTestProblemCCP7FC1New-0-bld.cc \
    $(SRCDIR)/AMRFDViewHeatTestProblemCCP7FC1New-1-bld.cc \
    $(SRCDIR)/AMRFDViewHeatTestProblemCCP7FC1New-2-bld.cc \

  BLDSRCS += \
    $(BLDDIR)/AMRFDViewScalarProblemCCP5FC0New-bld.cc \
    $(BLDDIR)/AMRFDViewScalarProblemCCP5FC1New-bld.cc \
    $(BLDDIR)/AMRFDViewScalarProblemCCP7FC0New-bld.cc \
    $(BLDDIR)/AMRFDViewScalarProblemCCP7FC1New-bld.cc \
    $(BLDDIR)/AMRFDViewHeatTestProblemCCP5FC0New-bld.cc \
    $(BLDDIR)/AMRFDViewHeatTestProblemCCP5FC1New-bld.cc \
    $(BLDDIR)/AMRFDViewHeatTestProblemCCP7FC0New-0-bld.cc \
    $(BLDDIR)/AMRFDViewHeatTestProblemCCP7FC0New-1-bld.cc \
    $(BLDDIR)/AMRFDViewHeatTestProblemCCP7FC0New-2-bld.cc \
    $(BLDDIR)/AMRFDViewHeatTestProblemCCP7FC1New-0-bld.cc \
    $(BLDDIR)/AMRFDViewHeatTestProblemCCP7FC1New-1-bld.cc \
    $(BLDDIR)/AMRFDViewHeatTestProblemCCP7FC1New-2-bld.cc \

endif # HAVE_HYPRE
