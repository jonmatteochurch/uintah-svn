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

# possible choices
VARs=(CC NC) # variable types
STNs=(P5 P7) # stencils
PBss=(
  "PureMetalProblem;HeatTestProblem;ScalarProblem"
  "PureMetalProblem;HeatTestProblem;ScalarProblem"
) # problems for each stencil
C2Fss=(
  "FC0;FC1;FC0New;FC1New;FCSimple;FCLinear;FCBilinear"
  "FC0;FC1;FC0New;FC1New;"
) # fine/coarse types for each stencil

# possible values
DIMs=(2 3) # dimension of stencil
DIRs=(x y z) # directions
SGNs=(minus plus) # signs
BCs=(Dirichlet Neumann) # boudary types
NFss=(
  "4;3;1"
  "4;3;1"
) # number of fieds of problems for each stencil

# generated source suffixes
PP="" # problem
VV="" # variable
SS="" # stencil
CC="" # fine/coarse
FF="" # field

# options
FIELD=-1
BC=true

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
  -c|--c2f)
  CC="$2"
  shift # past argument
  shift # past value
  ;;
  -B|--no-bc)
  BC=false
  shift # past argument
  ;;
  -f|--field)
  FIELD="$2"
  FF="-$2"
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
    PBss=("${PBss[s]}");
    NFss=("${NFss[s]}");
    C2Fss=("${C2Fss[s]}");
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
if [ -n "${CC}" ]; then
  not_found=true
  for ((s=0; s<${#STNs[@]}; s++)); do
    IFS=';' read -r -a C2Fs <<< "${C2Fss[s]}"
    if find_pos "${CC}" "${C2Fs[@]}" > /dev/null; then
      C2Fss[s]="${CC}";
      not_found=false
    else
      C2Fss[s]="";
    fi
  done
  if $not_found; then
    echo "invalid c2f ${CC}"
    exit -1
  fi
fi
if [ -n "${FF}" ]; then
  if [[ $FIELD -lt 0 ]]; then
    echo "invalid field ${FIELD}"
    exit -1
  fi
  for ((s=0; s<${#STNs[@]}; s++)); do
    IFS=';' read -r -a NFs <<< "${NFss[s]}"
    for NF in "${NFs[@]}"; do
      if [[ $FIELD -ge $NF ]]; then
        echo "invalid field ${FIELD}"
        exit -1
      fi
    done
  done
fi
if ! $BC; then
  BCs=()
fi

# flat PBss
PBf=()
for ((s=0; s<${#STNs[@]}; s++)); do
  IFS=';' read -r -a PBs <<< "${PBss[s]}"
  for PB in "${PBs[@]}"; do
    if ! find_pos "${PB}" "${PBf[@]}" > /dev/null; then
      PBf+=("${PB}")
    fi
  done
done

SRC=${0%-bld.sh}${PP}${VV}${SS}${CC}${FF}-bld.cc

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
for PB in ${PBf[@]}; do
  echo '#include <CCA/Components/PhaseField/DataTypes/'$PB'.h>' >> $SRC
done
echo '#include <CCA/Components/PhaseField/BoundaryConditions/BCFDView.h>' >> $SRC
echo '' >> $SRC
echo 'namespace Uintah {' >> $SRC
echo 'namespace PhaseField {' >> $SRC
echo '' >> $SRC

for VAR in ${VARs[@]}; do
  for ((s=0; s<${#STNs[@]}; s++)); do
    STN="${STNs[s]}"
    DIM="${DIMs[s]}"
    IFS=';' read -r -a PBs <<< "${PBss[s]}"
    IFS=';' read -r -a NFs <<< "${NFss[s]}"
    IFS=';' read -r -a C2Fs <<< "${C2Fss[s]}"
    for ((p=0; p<${#PBs[@]}; p++)); do
      PB="${PBs[p]}"
      NF="${NFs[p]}"
      PT="$PB<$VAR, $STN>"

      if [ "$PB" = "HeatTestProblem" ]; then
        echo "#ifdef PhaseField_Heat_DBG_DERIVATIVES" >> $SRC
      fi

      if [[ $FIELD -eq -1 ]]; then
        Is=`seq 0 $((NF-1))`
      else
        Is="$FIELD"
      fi
      for I in $Is; do
        for C2F in "${C2Fs[@]}"; do

          for ((d0=0; d0<$DIM; d0++)); do
            DIR0="${DIRs[d0]}"
            for SGN0 in "${SGNs[@]}"; do
              F0=$DIR0$SGN0

              # C2F on F0
              P0="Patch::$F0 | BC::FineCoarseInterface | FC::$C2F"
              echo "template<> const FactoryString BCFDView < $PT, $I, $P0 >::Name = \"$PB|$I|$VAR|$F0|$C2F|\";" >> $SRC

              for ((d1=d0+1; d1<$DIM; d1++)); do
                DIR1="${DIRs[d1]}"
                for SGN1 in "${SGNs[@]}"; do
                  F1=$DIR1$SGN1

                  # C2F on F0 and F1
                  P1="Patch::$F1 | BC::FineCoarseInterface | FC::$C2F"
                  echo "template<> const FactoryString BCFDView < $PT, $I, $P0, $P1 >::Name = \"$PB|$I|$VAR|$F0|$C2F|$F1|$C2F|\";" >> $SRC

                  for ((d2=d1+1; d2<$DIM; d2++)); do
                    DIR2="${DIRs[d2]}"
                    for SGN2 in "${SGNs[@]}"; do
                      F2=$DIR2$SGN2

                      # C2F on F0, F1 and F2
                      P2="Patch::$F2 | BC::FineCoarseInterface | FC::$C2F"
                      echo "template<> const FactoryString BCFDView < $PT, $I, $P0, $P1, $P2 >::Name = \"$PB|$I|$VAR|$F0|$C2F|$F1|$C2F|$F2|$C2F|\";" >> $SRC

                      for BC2 in "${BCs[@]}"; do

                        # C2F on F0 and F1; BC2 on F2
                        P2="Patch::$F2 | BC::$BC2"
                        echo "template<> const FactoryString BCFDView < $PT, $I, $P0, $P1, $P2 >::Name = \"$PB|$I|$VAR|$F0|$C2F|$F1|$C2F|$F2|$BC2|\";" >> $SRC

                      done # BC2

                    done
                  done # F2

                  for BC1 in "${BCs[@]}"; do

                    # C2F on F0; BC1 on F1
                    P1="Patch::$F1 | BC::$BC1"
                    echo "template<> const FactoryString BCFDView < $PT, $I, $P0, $P1 >::Name = \"$PB|$I|$VAR|$F0|$C2F|$F1|$BC1|\";" >> $SRC

                    for ((d2=d1+1; d2<$DIM; d2++)); do
                      DIR2="${DIRs[d2]}"
                      for SGN2 in "${SGNs[@]}"; do
                        F2=$DIR2$SGN2

                        # C2F on F0 and F2; BC1 on F1
                        P2="Patch::$F2 | BC::FineCoarseInterface | FC::$C2F"
                        echo "template<> const FactoryString BCFDView < $PT, $I, $P0, $P1, $P2 >::Name = \"$PB|$I|$VAR|$F0|$C2F|$F1|$BC1|$F2|$C2F|\";" >> $SRC

                        for BC2 in "${BCs[@]}"; do

                          # C2F on F0; BC1 on F1; BC2 on F2
                          P2="Patch::$F2 | BC::$BC2"
                          echo "template<> const FactoryString BCFDView < $PT, $I, $P0, $P1, $P2 >::Name = \"$PB|$I|$VAR|$F0|$C2F|$F1|$BC1|$F2|$BC2|\";" >> $SRC

                        done # BC2

                      done
                    done # F2

                  done # BC1

                done
              done # F1

              for BC0 in "${BCs[@]}"; do

                # BC0 on F0
                P0="Patch::$F0 | BC::$BC0"
                # no echo for BC only (already explicitly initialized)

                for ((d1=d0+1; d1<$DIM; d1++)); do
                  DIR1="${DIRs[d1]}"
                  for SGN1 in "${SGNs[@]}"; do
                    F1=$DIR1$SGN1

                    # BC0 on F0; C2F on F1;
                    P1="Patch::$F1 | BC::FineCoarseInterface | FC::$C2F"
                    echo "template<> const FactoryString BCFDView < $PT, $I, $P0, $P1 >::Name = \"$PB|$I|$VAR|$F0|$BC0|$F1|$C2F|\";" >> $SRC

                    for ((d2=d1+1; d2<$DIM; d2++)); do
                      DIR2="${DIRs[d2]}"
                      for SGN2 in "${SGNs[@]}"; do
                        F2=$DIR2$SGN2

                        # BC0 on F0; C2F on F1 and F2;
                        P2="Patch::$F2 | BC::FineCoarseInterface | FC::$C2F"
                        echo "template<> const FactoryString BCFDView < $PT, $I, $P0, $P1, $P2 >::Name = \"$PB|$I|$VAR|$F0|$BC0|$F1|$C2F|$F2|$C2F|\";" >> $SRC

                        for BC2 in "${BCs[@]}"; do

                          # BC0 on F0; C2F on F1; BC2 on F2;
                          P2="Patch::$F2 | BC::$BC2"
                          echo "template<> const FactoryString BCFDView < $PT, $I, $P0, $P1, $P2 >::Name = \"$PB|$I|$VAR|$F0|$BC0|$F1|$C2F|$F2|$BC2|\";" >> $SRC

                        done # BC2

                      done
                    done # F2

                    for BC1 in "${BCs[@]}"; do

                      # BC0 on F0; BC1 on F1
                      P1="Patch::$F1 | BC::$BC1"
                      # no echo for BC only (already explicitly initialized)

                      for ((d2=d1+1; d2<$DIM; d2++)); do
                        DIR2="${DIRs[d2]}"
                        for SGN2 in "${SGNs[@]}"; do
                          F2=$DIR2$SGN2

                          # BC0 on F0; BC1 on F1; C2F on F2
                          P2="Patch::$F2 | BC::FineCoarseInterface | FC::$C2F"
                          echo "template<> const FactoryString BCFDView < $PT, $I, $P0, $P1, $P2 >::Name = \"$PB|$I|$VAR|$F0|$BC0|$F1|$BC1|$F2|$C2F|\";" >> $SRC

                        done
                      done # F2

                    done # BC1

                  done
                done # F1

              done # BC0

            done
          done # F0

        done # C2F
      done # I

      if [ "$PB" = "HeatTestProblem" ]; then
        echo "#endif" >> $SRC
      fi

    done # PB
  done # STN
done # VAR

echo "" >> $SRC

for VAR in ${VARs[@]}; do
  for ((s=0; s<${#STNs[@]}; s++)); do
    STN="${STNs[s]}"
    DIM="${DIMs[s]}"
    IFS=';' read -r -a PBs <<< "${PBss[s]}"
    IFS=';' read -r -a NFs <<< "${NFss[s]}"
    IFS=';' read -r -a C2Fs <<< "${C2Fss[s]}"
    for ((p=0; p<${#PBs[@]}; p++)); do
      PB="${PBs[p]}"
      NF="${NFs[p]}"
      PT="$PB<$VAR, $STN>"

      if [ "$PB" = "HeatTestProblem" ]; then
        echo "#ifdef PhaseField_Heat_DBG_DERIVATIVES" >> $SRC
      fi

      if [[ $FIELD -eq -1 ]]; then
        Is=`seq 0 $((NF-1))`
      else
        Is="$FIELD"
      fi
      for I in $Is; do
        for C2F in "${C2Fs[@]}"; do

          for ((d0=0; d0<$DIM; d0++)); do
            DIR0="${DIRs[d0]}"
            for SGN0 in "${SGNs[@]}"; do
              F0=$DIR0$SGN0

              # C2F on F0
              P0="Patch::$F0 | BC::FineCoarseInterface | FC::$C2F"
              echo "template class BCFDView < $PT, $I, $P0 >;" >> $SRC

              for ((d1=d0+1; d1<$DIM; d1++)); do
                DIR1="${DIRs[d1]}"
                for SGN1 in "${SGNs[@]}"; do
                  F1=$DIR1$SGN1

                  # C2F on F0 and F1
                  P1="Patch::$F1 | BC::FineCoarseInterface | FC::$C2F"
                  echo "template class BCFDView < $PT, $I, $P0, $P1 >;" >> $SRC

                  for ((d2=d1+1; d2<$DIM; d2++)); do
                    DIR2="${DIRs[d2]}"
                    for SGN2 in "${SGNs[@]}"; do
                      F2=$DIR2$SGN2

                      # C2F on F0, F1 and F2
                      P2="Patch::$F2 | BC::FineCoarseInterface | FC::$C2F"
                      echo "template class BCFDView < $PT, $I, $P0, $P1, $P2 >;" >> $SRC

                      for BC2 in "${BCs[@]}"; do

                        # C2F on F0 and F1; BC2 on F2
                        P2="Patch::$F2 | BC::$BC2"
                        echo "template class BCFDView < $PT, $I, $P0, $P1, $P2 >;" >> $SRC

                      done # BC2

                    done
                  done # F2

                  for BC1 in "${BCs[@]}"; do

                    # C2F on F0; BC1 on F1
                    P1="Patch::$F1 | BC::$BC1"
                    echo "template class BCFDView < $PT, $I, $P0, $P1 >;" >> $SRC

                    for ((d2=d1+1; d2<$DIM; d2++)); do
                      DIR2="${DIRs[d2]}"
                      for SGN2 in "${SGNs[@]}"; do
                        F2=$DIR2$SGN2

                        # C2F on F0 and F2; BC1 on F1
                        P2="Patch::$F2 | BC::FineCoarseInterface | FC::$C2F"
                        echo "template class BCFDView < $PT, $I, $P0, $P1, $P2 >;" >> $SRC

                        for BC2 in "${BCs[@]}"; do

                          # C2F on F0; BC1 on F1; BC2 on F2
                          P2="Patch::$F2 | BC::$BC2"
                          echo "template class BCFDView < $PT, $I, $P0, $P1, $P2 >;" >> $SRC

                        done # BC2

                      done
                    done # F2

                  done # BC1

                done
              done # F1

              for BC0 in "${BCs[@]}"; do

                # BC0 on F0
                P0="Patch::$F0 | BC::$BC0"
                # no echo for BC only (already explicitly initialized)

                for ((d1=d0+1; d1<$DIM; d1++)); do
                  DIR1="${DIRs[d1]}"
                  for SGN1 in "${SGNs[@]}"; do
                    F1=$DIR1$SGN1

                    # BC0 on F0; C2F on F1;
                    P1="Patch::$F1 | BC::FineCoarseInterface | FC::$C2F"
                    echo "template class BCFDView < $PT, $I, $P0, $P1 >;" >> $SRC

                    for ((d2=d1+1; d2<$DIM; d2++)); do
                      DIR2="${DIRs[d2]}"
                      for SGN2 in "${SGNs[@]}"; do
                        F2=$DIR2$SGN2

                        # BC0 on F0; C2F on F1 and F2;
                        P2="Patch::$F2 | BC::FineCoarseInterface | FC::$C2F"
                        echo "template class BCFDView < $PT, $I, $P0, $P1, $P2 >;" >> $SRC

                        for BC2 in "${BCs[@]}"; do

                          # BC0 on F0; C2F on F1; BC2 on F2;
                          P2="Patch::$F2 | BC::$BC2"
                          echo "template class BCFDView < $PT, $I, $P0, $P1, $P2 >;" >> $SRC

                        done # BC2

                      done
                    done # F2

                    for BC1 in "${BCs[@]}"; do

                      # BC0 on F0; BC1 on F1
                      P1="Patch::$F1 | BC::$BC1"
                      # no echo for BC only (already explicitly initialized)

                      for ((d2=d1+1; d2<$DIM; d2++)); do
                        DIR2="${DIRs[d2]}"
                        for SGN2 in "${SGNs[@]}"; do
                          F2=$DIR2$SGN2

                          # BC0 on F0; BC1 on F1; C2F on F2
                          P2="Patch::$F2 | BC::FineCoarseInterface | FC::$C2F"
                          echo "template class BCFDView < $PT, $I, $P0, $P1, $P2 >;" >> $SRC

                        done
                      done # F2

                    done # BC1

                  done
                done # F1

              done # BC0

            done
          done # F0

        done # C2F
      done # I

      if [ "$PB" = "HeatTestProblem" ]; then
        echo "#endif" >> $SRC
      fi

    done # PB
  done # STN
done # VAR

echo "" >> $SRC

echo '} // namespace Uintah' >> $SRC
echo '} // namespace PhaseField' >> $SRC
