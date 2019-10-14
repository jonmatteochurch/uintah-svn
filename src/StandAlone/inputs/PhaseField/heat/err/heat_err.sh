#!/bin/bash
set -e

VARs=(cc nc)
DIMs=(2 3)
Rs=(1 2 4 8 16 32 64 128 256 512)
TSs=(be cn)
SCHs=(backward_euler crank_nicolson)

for VAR in ${VARs[@]}; do
  for DIM in ${DIMs[@]}; do
    for ((n=0; n<${#Rs[@]}; n++)); do
      for ((m=0; m<${#Rs[@]}; m++)); do
        R=${Rs[n]}
        N=$(( 16 * R )) 
        M=${Rs[m]}
        K=$(bc <<< "scale=12; ( 1 / $M )");

        if [[ $DIM -eq 3 ]]; then
          ONE="[1,1,1]"
          TWO="[2,2,2]"
          SXT="[16,16,16]"
          DST="[ 16., 16., 16.]"
          RES="[$N,$N,$N]"
          BCZ=$(cat << BCZ_END
            <Face side="z-">\n
                <BCType id="0" label="u" var="Neumann">\n
                    <value>0.</value>\n
                </BCType>\n
            </Face>\n
            <Face side="z+">\n
                <BCType id="0" label="u" var="Dirichlet">\n
                    <value>0.</value>\n
                </BCType>\n
            </Face>\n
BCZ_END
          );
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

        TIT=$(printf "heat_err_%s_%1dd_fe_n%03d_m%03d" $VAR $DIM $N $M)

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
              /<!--Solver-->/d" heat_err.template > $TIT.ups

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

          TIT=$(printf "heat_err_%s_%1dd_%s_n%03d_m%03d" $VAR $DIM $TS $N $M)

          sed "s|<!--title-->|<title>$TIT</title>|g;
               s|<!--var-->|<var>$VAR</var>|g;
               s|<!--dim-->|<dim>$DIM</dim>|g;
               s|<!--delt-->|<delt>$K</delt>|g;
               s|<!--scheme-->|<scheme>$SCH</scheme>|g;
               s|<!--upper-->|<upper>$DST</upper>|g;
               s|<!--resolution-->|<resolution>$RES</resolution>|g;
               s|<!--patches-->|<patches>$TWO</patches>|g;
               s|<!--filebase-->|<filebase>$TIT.uda</filebase>|g;
               s|<!--BC Z Faces-->|$BCZ|g;
               s|<!--Solver-->|$SLV|g" heat_err.template > $TIT.ups

        done
      fi
      done
    done
  done
done
