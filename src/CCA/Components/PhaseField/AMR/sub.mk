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
  $(SRCDIR)/AMRFDViewScalarProblemCCFC0-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemCCFC1-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemNCFC0-bld.cc \
  $(SRCDIR)/AMRFDViewScalarProblemNCFC1-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCFC0-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCFC1-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCFCSimple-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCFCLinear-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemCCFCBilinear-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCFC0-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCFC1-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCFCSimple-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCFCLinear-bld.cc \
  $(SRCDIR)/AMRFDViewHeatTestProblemNCFCBilinear-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemCCFC0-bld.cc \
  $(SRCDIR)/AMRFDViewPureMetalProblemCCFC1-bld.cc \
  
BLDDIR := $(SRCTOP)/$(SRCDIR)
BLDDEPS := $(SRCDIR)/AMRFDView-bld.sh

BLDSRCS += \
  $(BLDDIR)/AMRFDViewScalarProblemCCFC0-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemCCFC1-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemNCFC0-bld.cc \
  $(BLDDIR)/AMRFDViewScalarProblemNCFC1-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCFC0-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCFC1-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCFCSimple-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCFCLinear-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemCCFCBilinear-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCFC0-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCFC1-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCFCSimple-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCFCLinear-bld.cc \
  $(BLDDIR)/AMRFDViewHeatTestProblemNCFCBilinear-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemCCFC0-bld.cc \
  $(BLDDIR)/AMRFDViewPureMetalProblemCCFC1-bld.cc \

ifeq ($(HAVE_HYPRE),yes)
  SRCS += \
    $(SRCDIR)/AMRFDViewScalarProblemCCFC0New-bld.cc \
    $(SRCDIR)/AMRFDViewScalarProblemCCFC1New-bld.cc \
    $(SRCDIR)/AMRFDViewHeatTestProblemCCFC0New-bld.cc \
    $(SRCDIR)/AMRFDViewHeatTestProblemCCFC1New-bld.cc \
    $(SRCDIR)/AMRFDViewPureMetalProblemCCFC0New-bld.cc \
    $(SRCDIR)/AMRFDViewPureMetalProblemCCFC1New-bld.cc \

  BLDSRCS += \
    $(BLDDIR)/AMRFDViewScalarProblemCCFC0New-bld.cc \
    $(BLDDIR)/AMRFDViewScalarProblemCCFC1New-bld.cc \
    $(BLDDIR)/AMRFDViewHeatTestProblemCCFC0New-bld.cc \
    $(BLDDIR)/AMRFDViewHeatTestProblemCCFC1New-bld.cc \
    $(BLDDIR)/AMRFDViewPureMetalProblemCCFC0New-bld.cc \
    $(BLDDIR)/AMRFDViewPureMetalProblemCCFC1New-bld.cc \

endif # HAVE_HYPRE
