#!/bin/bash
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

PSLVs=(PCG GMRES FlexGMRES LGMRES BiCGSTAB)
SSLVs=(SysPFMG Split FAC)
PCNs=(Diagonal SysPFMG Split FAC PCG GMRES FlexGMRES LGMRES BiCGSTAB)
DIMs=(2 3)
C2Fs=(0 1)

SRC=${0%-bld.sh}-bld.cc

echo "generating $SRC"

echo '/*' > $SRC
echo ' * The MIT License' >> $SRC
echo ' *' >> $SRC
echo ' * Copyright (c) 1997-2020 The University of Utah' >> $SRC
echo ' *' >> $SRC
echo ' * Permission is hereby granted, free of charge, to any person obtaining a copy' >> $SRC
echo ' * of this software and associated documentation files (the "Software"), to' >> $SRC
echo ' * deal in the Software without restriction, including without limitation the' >> $SRC
echo ' * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or' >> $SRC
echo ' * sell copies of the Software, and to permit persons to whom the Software is' >> $SRC
echo ' * furnished to do so, subject to the following conditions:' >> $SRC
echo ' *' >> $SRC
echo ' * The above copyright notice and this permission notice shall be included in' >> $SRC
echo ' * all copies or substantial portions of the Software.' >> $SRC
echo ' *' >> $SRC
echo ' * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR' >> $SRC
echo ' * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,' >> $SRC
echo ' * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE' >> $SRC
echo ' * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER' >> $SRC
echo ' * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING' >> $SRC
echo ' * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS' >> $SRC
echo ' * IN THE SOFTWARE.' >> $SRC
echo ' */' >> $SRC
echo '' >> $SRC
echo '#include <CCA/Components/Solvers/HypreSStruct/SStructSolver.h>' >> $SRC
echo '' >> $SRC
echo 'namespace Uintah {' >> $SRC
echo 'namespace HypreSStruct {' >> $SRC
echo '' >> $SRC

for DIM in ${DIMs[@]}; do
  for C2F in ${C2Fs[@]}; do
    for SLV in ${SSLVs[@]} ${PSLVs[@]}; do
      echo "template<> const FactoryString SStructSolver<S::$SLV, $DIM, $C2F>::Name = \"${SLV,,}|$DIM|$C2F\";" >> $SRC
    done
    for SLV in ${PSLVs[@]}; do
      for PCN in ${PCNs[@]}; do
        if [ "$PCN" != "$SLV" ]; then
          echo "template<> const FactoryString SStructSolver<S::$SLV, $DIM, $C2F, P::$PCN>::Name = \"${SLV,,}|${PCN,,}|$DIM|$C2F\";" >> $SRC
        fi
      done
    done
  done
done

echo "" >> $SRC

for DIM in ${DIMs[@]}; do
  for C2F in ${C2Fs[@]}; do
    for SLV in ${SSLVs[@]} ${PSLVs[@]}; do
      echo "template class SStructSolver<S::$SLV, $DIM, $C2F>;" >> $SRC
    done
    for SLV in ${PSLVs[@]}; do
      for PCN in ${PCNs[@]}; do
        if [ "$PCN" != "$SLV" ]; then
          echo "template class SStructSolver<S::$SLV, $DIM, $C2F, P::$PCN>;" >> $SRC
        fi
      done
    done
  done
done

echo "" >> $SRC

echo '} // namespace HypreSStruct' >> $SRC
echo '} // namespace Uintah' >> $SRC
