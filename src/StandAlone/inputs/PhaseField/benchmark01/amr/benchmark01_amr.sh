#!/bin/bash
set -e
VARs=(cc nc)
EPSs=(0.2 0.1 0.05)
MXTs=(50 220 850)
C2Fs=(FC0 FC1 FCSimple FCLinear FCBilinear)
Rs=(64 128 256 512)
Pe=4
TSs=(si0)
SCHs=(semi_implicit_0)
for VAR in ${VARs[@]}; do
  for ((e=0; e<${#EPSs[@]}; e++)); do
    EPS=${EPSs[e]}
    MXT=${MXTs[e]}
    for ((r=0; r<${#Rs[@]}; r++)); do
      R=${Rs[r]}
      N=$(( 16 * R )) 
      M=$(( 2 * Pe * R ))
      K=$(bc <<< "scale=12; ( 1 / $M )");
      ONE="[1,1,0]"
      TWO="[2,2,1]"
      SXT="[16,16,1]"
      RES="[$N,$N,1.]"
# explicit no amr
      TIT=$(printf "benchmark01_amr_eps_%1.2f_%s_fe_n%04d_l1" $EPS $VAR $N)
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
           /<!--Solver-->/d" benchmark01_amr.template > $TIT.ups
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
          TIT=$(printf "benchmark01_amr_eps%1.2f_%s_%s_hypre_n%04d_l1" $EPS $VAR $TS $N)
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
               /<!--AMR-->/d" benchmark01_amr.template > $TIT.ups
        done
      fi
      for ((L=2; L<6-r; L++)); do
# explicit with amr
        K=$(bc <<< "scale=12; ( $K / 2 )");
        for FCI in "${C2Fs[@]}"; do
          TIT=$(printf "benchmark01_amr_eps%1.2f_%s_fe_n%04d_l%1d_%s" $EPS $VAR $N $L ${FCI,,})
          AMR=$(cat << AMR_END0
<!--__________________________________-->\n
    <AMR>\n
        <Regridder type="Tiled">\n
            <adaptive>true</adaptive>\n
            <max_levels>$L</max_levels>\n
            <cell_refinement_ratio>[$TWO]</cell_refinement_ratio>\n
            <cell_stability_dilation>$ONE</cell_stability_dilation>\n
            <min_boundary_cells>$ONE</min_boundary_cells>\n
            <min_patch_size>[$SXT]</min_patch_size>\n
        </Regridder>\n
        <FineCoarseInterfaces>\n
            <FCIType id="0" label="u" var="$FCI" />\n
        </FineCoarseInterfaces>\n
    </AMR>
AMR_END0
          )
          AMR=$(tr -d '\n' <<< "$AMR")
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
               s|<!--AMR-->|$AMR|g; 
               /<!--Solver-->/d" benchmark01_amr.template > $TIT.ups
        done
        if [ "$VAR" == "cc" ]; then
# hypre with amr
          Fs=("${C2Fs[@]:0:2}")
          SLV=$(cat << SLV_END1
<!--__________________________________-->\n
    <Solver type="hypre">\n
       <Parameters variable="u">\n
          <solver>smg</solver>\n
          <preconditioner>diagonal</preconditioner>\n
          <tolerance>1.e-6</tolerance>\n
          <maxiterations>20</maxiterations>\n
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
SLV_END1
          )
          SLV=$(tr -d '\n' <<< "$SLV")
          for FCI in "${Fs[@]}"; do
            AMR=$(cat << AMR_END1
<!--__________________________________-->\n
    <AMR>\n
        <Regridder type="Tiled">\n
            <adaptive>true</adaptive>\n
            <max_levels>$L</max_levels>\n
            <cell_refinement_ratio>[$TWO]</cell_refinement_ratio>\n
            <cell_stability_dilation>$ONE</cell_stability_dilation>\n
            <min_boundary_cells>$ONE</min_boundary_cells>\n
            <min_patch_size>[$SXT]</min_patch_size>\n
        </Regridder>\n
        <FineCoarseInterfaces>\n
            <FCIType id="0" label="u" var="${FCI}New" />\n
        </FineCoarseInterfaces>\n
    </AMR>
AMR_END1
            )
            AMR=$(tr -d '\n' <<< "$AMR")
            for ((t=0; t<${#TSs[@]}; t++)); do
              TS=${TSs[t]}
              SCH=${SCHs[t]}
              TIT=$(printf "benchmark01_amr_eps%1.2f_%s_%s_hypre_n%04d_l%1d_%snew" $EPS $VAR $TS $N $L ${FCI,,})
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
                   s|<!--AMR-->|$AMR|g;
                   s|<!--Solver-->|$SLV|g" benchmark01_amr.template > $TIT.ups
            done
          done
# hypre_sstruct with amr
          Fs=("${C2Fs[@]}")
          SLV=$(cat << SLV_END2
<!--__________________________________-->\n
    <Solver type="hypre_sstruct" ndim="2" >\n
       <Parameters variable="u">\n
          <solveFrequency>1</solveFrequency>\n
          <setupFrequency>0</setupFrequency>\n
          <updateCoefFrequency>0</updateCoefFrequency>\n
          <max_levels>$L</max_levels>\n
          <maxiterations>1000</maxiterations>\n
          <tolerance>1.e-6</tolerance>\n
          <rel_change>0</rel_change>\n
          <relax_type>2</relax_type>\n
          <npre>1</npre>\n
          <npost>1</npost>\n
          <csolver_type>2</csolver_type>\n
          <logging>1</logging>\n
       </Parameters>\n
    </Solver>
SLV_END2
          )
          SLV=$(tr -d '\n' <<< "$SLV")
          for FCI in "${Fs[@]}"; do
            AMR=$(cat << AMR_END2
<!--__________________________________-->\n
    <AMR>\n
        <Regridder type="Tiled">\n
            <adaptive>true</adaptive>\n
            <max_levels>$L</max_levels>\n
            <cell_refinement_ratio>[$TWO]</cell_refinement_ratio>\n
            <cell_stability_dilation>$ONE</cell_stability_dilation>\n
            <min_boundary_cells>$ONE</min_boundary_cells>\n
            <min_patch_size>[$SXT]</min_patch_size>\n
        </Regridder>\n
        <FineCoarseInterfaces>\n
            <FCIType id="0" label="u" var="$FCI" />\n
        </FineCoarseInterfaces>\n
    </AMR>
AMR_END2
            )
            AMR=$(tr -d '\n' <<< "$AMR")
            for ((t=0; t<${#TSs[@]}; t++)); do
              TS=${TSs[t]}
              SCH=${SCHs[t]}
              TIT=$(printf "benchmark01_amr_eps%1.2f_%s_%s_hypre_sstruct_n%04d_l%1d_%s" $EPS $VAR $TS $N $L ${FCI,,})
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
                   s|<!--AMR-->|$AMR|g;
                   s|<!--Solver-->|$SLV|g" benchmark01_amr.template > $TIT.ups
            done
          done
        fi
      done
    done
  done
done
