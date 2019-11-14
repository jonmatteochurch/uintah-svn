#!/bin/bash
set -e
VARs=(cc nc)
DIMs=(2 3)
C2Fs=(FC0 FC1 FCSimple FCLinear FCBilinear)
Rs=(1 2 4 8 16)
Pe=2
TSs=(be cn)
SCHs=(backward_euler crank_nicolson)
for VAR in ${VARs[@]}; do
  for DIM in ${DIMs[@]}; do
    for ((r=0; r<${#Rs[@]}; r++)); do
      R=${Rs[r]}
      N=$(( 16 * R )) 
      M=$(( 2 * Pe * DIM * R * R ))
      H=$(bc <<< "scale=4;  ( 1 / $R )");
      K=$(bc <<< "scale=12; ( 1 / $M )");
      if [[ $DIM -eq 3 ]]; then
        ONE="[1,1,1]"
        TWO="[2,2,2]"
        SXT="[16,16,16]"
        DST="[ 16., 16., 16.]"
        RES="[$N,$N,$N]"
        BCZ=$(cat << BCZ_END
            <Face side="z-">
                <BCType id="0" label="u" var="Neumann">
                    <value>0.</value>
                </BCType>
            </Face>
            <Face side="z+">
                <BCType id="0" label="u" var="Dirichlet">
                    <value>0.</value>
                </BCType>
            </Face>
BCZ_END
        )
      else
        ONE="[1,1,0]"
        TWO="[2,2,1]"
        SXT="[16,16,1]"
        DST="[ 16., 16.,  1.]"
        RES="[$N,$N,1]"
        BCZ=""
      fi
      BCZ=$(tr -d '\n' <<< "$BCZ")
# explicit no amr
      TIT=$(printf "heat_amr_%s_%1dd_fe_n%04d_l1" $VAR $DIM $N)
      sed "s|<!--title-->|<title>$TIT</title>|g;
           s|<!--var-->|<var>$VAR</var>|g;
           s|<!--dim-->|<dim>$DIM</dim>|g;
           s|<!--delt-->|<delt>$K</delt>|g;
           s|<!--scheme-->|<scheme>forward_euler</scheme>|g;
           s|<!--upper-->|<upper>$DST</upper>|g;
           s|<!--resolution-->|<resolution>$RES</resolution>|g;
           s|<!--patches-->|<patches>$TWO</patches>|g;
           s|<!--filebase-->|<filebase>$TIT.uda</filebase>|g;
           s|<!--BC Z Faces-->|$BCZ|g;
           /<!--AMR-->/d;
           /<!--Solver-->/d" heat_amr.template > $TIT.ups
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
          TIT=$(printf "heat_amr_%s_%1dd_%s_hypre_n%04d_l1" $VAR $DIM $TS $N)
          sed "s|<!--title-->|<title>$TIT</title>|g;
               s|<!--var-->|<var>$VAR</var>|g;
               s|<!--dim-->|<dim>$DIM</dim>|g;
               s|<!--delt-->|<delt>1.</delt>|g;
               s|<!--scheme-->|<scheme>$SCH</scheme>|g;
               s|<!--upper-->|<upper>$DST</upper>|g;
               s|<!--resolution-->|<resolution>$RES</resolution>|g;
               s|<!--patches-->|<patches>$TWO</patches>|g;
               s|<!--filebase-->|<filebase>$TIT.uda</filebase>|g;
               s|<!--BC Z Faces-->|$BCZ|g;
               s|<!--Solver-->|$SLV|g;
               /<!--AMR-->/d" heat_amr.template > $TIT.ups
        done
      fi
      for ((L=2; L<5-r; L++)); do
# explicit with amr
        K=$(bc <<< "scale=12; ( $K / 4 )");
        if [[ $DIM -eq 3 ]]; then
          Fs=("${C2Fs[@]:0:2}")
        else
          Fs=("${C2Fs[@]}")
        fi
        for FCI in "${Fs[@]}"; do
          TIT=$(printf "heat_amr_%s_%1dd_fe_n%04d_l%1d_%s" $VAR $DIM $N $L ${FCI,,})
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
               s|<!--dim-->|<dim>$DIM</dim>|g;
               s|<!--delt-->|<delt>$K</delt>|g;
               s|<!--scheme-->|<scheme>forward_euler</scheme>|g;
               s|<!--upper-->|<upper>$DST</upper>|g;
               s|<!--resolution-->|<resolution>$RES</resolution>|g;
               s|<!--patches-->|<patches>$TWO</patches>|g;
               s|<!--filebase-->|<filebase>$TIT.uda</filebase>|g;
               s|<!--BC Z Faces-->|$BCZ|g;
               s|<!--AMR-->|$AMR|g; 
               /<!--Solver-->/d" heat_amr.template > $TIT.ups
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
              TIT=$(printf "heat_amr_%s_%1dd_%s_hypre_n%04d_l%1d_%snew" $VAR $DIM $TS $N $L ${FCI,,})
              sed "s|<!--title-->|<title>$TIT</title>|g;
                   s|<!--var-->|<var>$VAR</var>|g;
                   s|<!--dim-->|<dim>$DIM</dim>|g;
                   s|<!--delt-->|<delt>1.</delt>|g;
                   s|<!--scheme-->|<scheme>$SCH</scheme>|g;
                   s|<!--upper-->|<upper>$DST</upper>|g;
                   s|<!--resolution-->|<resolution>$RES</resolution>|g;
                   s|<!--patches-->|<patches>$TWO</patches>|g;
                   s|<!--filebase-->|<filebase>$TIT.uda</filebase>|g;
                   s|<!--BC Z Faces-->|$BCZ|g;
                   s|<!--AMR-->|$AMR|g;
                   s|<!--Solver-->|$SLV|g" heat_amr.template > $TIT.ups
            done
          done
# hypre_sstruct with amr
          Fs=("${C2Fs[@]}")
          SLV=$(cat << SLV_END2
<!--__________________________________-->\n
    <Solver type="hypre_sstruct" ndim="$DIM" >\n
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
              TIT=$(printf "heat_amr_%s_%1dd_%s_hypre_sstruct_n%04d_l%1d_%s" $VAR $DIM $TS $N $L ${FCI,,})
              sed "s|<!--title-->|<title>$TIT</title>|g;
                   s|<!--var-->|<var>$VAR</var>|g;
                   s|<!--dim-->|<dim>$DIM</dim>|g;
                   s|<!--delt-->|<delt>1.</delt>|g;
                   s|<!--scheme-->|<scheme>$SCH</scheme>|g;
                   s|<!--upper-->|<upper>$DST</upper>|g;
                   s|<!--resolution-->|<resolution>$RES</resolution>|g;
                   s|<!--patches-->|<patches>$TWO</patches>|g;
                   s|<!--filebase-->|<filebase>$TIT.uda</filebase>|g;
                   s|<!--BC Z Faces-->|$BCZ|g;
                   s|<!--AMR-->|$AMR|g;
                   s|<!--Solver-->|$SLV|g" heat_amr.template > $TIT.ups
            done
          done
        fi
      done
    done
  done
done
