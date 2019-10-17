#!/bin/bash
set -e

VARs=(cc nc)
EPSs=(0.2 0.1 0.05)
MXTs=(50 220 850)
Rs=(1 2 4 8 16 32 64 128 256 512)
Pe=4
TSs=(si0)
SCHs=(semi_implicit_0)

for VAR in ${VARs[@]}; do
  for ((e=0; e<${#EPSs[@]}; e++)); do
    EPS=${EPSs[e]}
    MXT=${MXTs[e]}
    for ((n=0; n<${#Rs[@]}; n++)); do
      for ((m=0; m<${#Rs[@]}; m++)); do
        R=${Rs[n]}
        N=$(( 16 * R )) 
        M=${Rs[m]}
        K=$(bc <<< "scale=12; ( 1 / $M )");

        ONE="[1,1,0]"
        TWO="[2,2,1]"
        SXT="[16,16,1]"
        RES="[$N,$N,1]"

# explicit no amr

        TIT=$(printf "benchmark01_err_eps%1.2f_%s_fe_n%04d_m%04d" $EPS $VAR $N $M)

        sed "s|<!--title-->|<title>$TIT</title>|g;
             s|<!--var-->|<var>$VAR</var>|g;
             s|<!--delt-->|<delt>$K</delt>|g;
             s|<!--epsilon-->|<epsilon>$EPS</epsilon>|g;
             s|<!--scheme-->|<scheme>forward_euler</scheme>|g;
             s|<!--maxTime-->|<maxTime>$(bc <<< "scale=12; ( $MXT + $K )")</maxTime>|g;
             s|<!--resolution-->|<resolution>$RES</resolution>|g;
             s|<!--patches-->|<patches>$TWO</patches>|g;
             s|<!--filebase-->|<filebase>$TIT.uda</filebase>|g;
             s|<!--outputInterval-->|<outputInterval>$MXT</outputInterval>|g;
             /<!--AMR-->/d;
             /<!--Solver-->/d" benchmark01_err.template > $TIT.ups

        if [ "$VAR" == "cc" ]; then

# hypre no amr

          SLV=$(cat << SLV_END0
<!--__________________________________-->\n
    <Solver type="hypre">\n
       <Parameters variable="u">\n
          <solver>smg</solver>\n
          <preconditioner>diagonal</preconditioner>\n
          <tolerance>1.e-6</tolerance>\n
          <maxiterations>1000000</maxiterations>\n
          <precond_maxiters>1</precond_maxiters>\n
          <precond_tolerance>0</precond_tolerance>\n
          <npre>1</npre>\n
          <npost>1</npost>\n
          <skip>0</skip>\n
          <jump>0</jump>\n
          <logging>1</logging>\n
          <setupFrequency>0</setupFrequency>\n
          <updateCoefFrequency>0</updateCoefFrequency>\n
          <solveFrequency>1</solveFrequency>\n
          <relax_type>2</relax_type>\n
       </Parameters>\n
    </Solver>
SLV_END0
          )
          SLV=$(tr -d '\n' <<< "$SLV")

          for ((t=0; t<${#TSs[@]}; t++)); do
            TS=${TSs[t]}
            SCH=${SCHs[t]}

            TIT=$(printf "benchmark01_err_eps%1.2f_%s_%s_hypre_n%04d_m%04d" $EPS $VAR $TS $N $M)

            sed "s|<!--title-->|<title>$TIT</title>|g;
                s|<!--var-->|<var>$VAR</var>|g;
                s|<!--delt-->|<delt>1.</delt>|g;
                s|<!--epsilon-->|<epsilon>$EPS</epsilon>|g;
                s|<!--scheme-->|<scheme>$SCH</scheme>|g;
                s|<!--maxTime-->|<maxTime>$(bc <<< "scale=12; ( $MXT + 1 )")</maxTime>|g;
                s|<!--resolution-->|<resolution>$RES</resolution>|g;
                s|<!--patches-->|<patches>$TWO</patches>|g;
                s|<!--filebase-->|<filebase>$TIT.uda</filebase>|g;
                s|<!--outputInterval-->|<outputInterval>$MXT</outputInterval>|g;
                s|<!--Solver-->|$SLV|g;
                /<!--AMR-->/d" benchmark01_err.template > $TIT.ups

          done
        fi
      done
    done
  done
done
