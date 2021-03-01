ARC = """
#$ -cwd -V
#$ -m be -M scjmc@leeds.ac.uk
#$ -pe ib %(np)d
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G

if [ -z "$UPATH" ]; then
    UPATH=$HOME/uintah/trunk
fi

SUS=$UPATH/arc.dbgder/build/StandAlone/sus
INP=$UPATH/src/StandAlone/inputs/PhaseField/heat/cn_struct_conv

for UPS in $INP/cn_struct_conv_%(solver)s_2d_nlvl*_%(ph)d_*.ups; do
    LOG=$(basename $UPS).log
    [[ -f $LOG && -s $LOG ]] || mpirun -np %(np)d $SUS $UPS &> $LOG
done
"""

NAME = "cn_struct_conv_%(solver)s_2d_np%(np)d"

for solver in ("pfmg","smg","cycred","gmres","flexgmres","lgmres","bicgstab"):
	for nlvl in range(1,9):
		for ph in range(0,5):
			for pk in range(0,5):
				np = 2**ph
				name = NAME % { "solver": solver, "np": np };
				arc = open(name+".arc4", "w")
				arc.write(ARC % { "np": np, "solver": solver, "ph": ph })
				arc.close()
