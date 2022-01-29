#!/bin/bash
#
#  The MIT License
#
#  Copyright (c) 1997-2022 The University of Utah
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

# possible choices
VARs=(CC NC) # variable types
STNs=(P3 P5 P7) # stencils
PBss=(
  "HeatTestProblem;ScalarProblem"
  "PureMetalProblem;HeatTestProblem;ScalarProblem"
  "PureMetalProblem;HeatTestProblem;ScalarProblem"
) # problems for each stencil

# possible values
DIMs=(2 3) # dimension of stencil
CMBs=(1 3) # combination of 
DIRs=(x y z) # directions
SGNs=(minus plus) # signs
BCs=(Dirichlet Neumann) # boudary types
NFss=(
  "4"
  "4;3;1"
  "4;3;1"
) # number of fieds of problems for each stencil

declare -A Fss
Fss[PureMetalProblem]='ScalarField<const double>;VectorField<const double, $CC>'
Fss[HeatTestProblem]='ScalarField<const double>;VectorField<const double, $DD>::value>'
Fss[ScalarProblem]='ScalarField<const double>'

# generated source suffixes
PP="" # problem
VV="" # variable
SS="" # stencil

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
  -v|--var)
  VV="$2"
  shift # past argument
  shift # past value
  ;;
  -s|--stn)
  SS="$2"
  shift # past argument
  shift # past value
  ;;
  -p|--pb)
  PP="$2"
  shift # past argument
  shift # past value
  ;;
  *)    # unknown option
  POSITIONAL+=("$1") # save it in an array for later
  shift # past argument
  ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

function find_pos() {
  local value="$1"
  shift
  local array=($@)
  for ((i=0; i<${#array[@]}; i++)); do
    if [ "${value}" == "${array[i]}" ]; then
      echo $i
      return 0
    fi
  done
  return 1
}

# filter possible values
if [ -n "${VV}" ]; then
  if find_pos "${VV}" "${VARs[@]}" > /dev/null; then
    VARs=("$VV")
  else
    echo "invalid var ${VV}"
    exit -1
  fi
fi
if [ -n "${SS}" ]; then
  s=$(find_pos "${SS}" "${STNs[@]}")
  if [ "$?" -eq 0 ]; then
    STNs=(${STNs[s]});
    DIMs=(${DIMs[s]});
    CMBs=(${CMBs[s]});
    PBss=("${PBss[s]}");
    NFss=("${NFss[s]}");
  else
    echo "invalid stencil ${SS}"
    exit -1
  fi
fi
if [ -n "${PP}" ]; then
  not_found=true
  for ((s=0; s<${#STNs[@]}; s++)); do
    IFS=';' read -r -a PBs <<< "${PBss[s]}"
    IFS=';' read -r -a NFs <<< "${NFss[s]}"
    p=$(find_pos "${PP}" "${PBs[@]}")
    if [ "$?" -eq 0 ]; then
      PBss[s]="${PBs[p]}";
      NFss[s]="${NFs[p]}";
      not_found=false
    else
      PBss[s]="";
      NFss[s]="";
    fi
  done
  if $not_found; then
    echo "invalid problem ${PP}"
    exit -1
  fi
fi

# flat PBss / Fss
PBf=()
Ffs=()
for ((s=0; s<${#STNs[@]}; s++)); do
  IFS=';' read -r -a PBs <<< "${PBss[s]}"
  DD="D${DIMs[s]}"
  CC="${CMBs[s]}"
  Ffs[s]=""
  Ff=()
  for PB in ${PBs[@]}; do
    IFS=';' read -r -a Fs <<< "${Fss[$PB]}"
    if ! find_pos "${PB}" "${PBf[@]}" > /dev/null; then
      PBf+=("${PB}")
    fi
    for FF in ${Fs[@]}; do
      F=eval echo "${FF}"
      if ! find_pos "${F}" "${Ff[@]}" > /dev/null; then
        Ff+=("${F}")
      fi
    done
    for F in ${Ff[@]}; do
      Ffs[s] += "$F;"
    done
  done
done

SRC=${0%-bld.sh}${PP}${VV}${SS}${FF}-bld.cc

echo "generating $SRC"

echo '/*' > $SRC
echo ' * The MIT License' >> $SRC
echo ' *' >> $SRC
echo ' * Copyright (c) 1997-2022 The University of Utah' >> $SRC
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
for PB in ${PBf[@]}; do
  echo '#include <CCA/Components/PhaseField/DataTypes/'$PB'.h>' >> $SRC
done
echo '#include <CCA/Components/PhaseField/BoundaryConditions/BCFDView.h>' >> $SRC
echo '#include <CCA/Components/PhaseField/BoundaryConditions/BCFDViewFactory.h>' >> $SRC
echo '' >> $SRC
echo 'namespace Uintah {' >> $SRC
echo '' >> $SRC

for ((s=0; s<${#STNs[@]}; s++)); do
  IFS=';' read -r -a Ff <<< "${Ffs[$s]}"
  for F in ${Ff[@]}; do
    echo "template class Factory < PhaseField::FDView< PhaseField::$F, PhaseField::$STN >, const typename PhaseField::$F::label_type &, const VarLabel *, int, const Level *, const std::vector< PhaseField::BCInfo< PhaseField::$F > > & >;" >> $SRC
  done
  echo "" >> $SRC
done

echo '} // namespace Uintah' >> $SRC
