ARC = """
#$ -cwd -V
#$ -m be -M scjmc@leeds.ac.uk
#$ -pe ib %(np)d
#$ -l h_rt=12:00:00
#$ -l h_vmem=4.8G

if [ -z "$UPATH" ]; then
    UPATH=$HOME/uintah/trunk
fi

SUS=$UPATH/arc/build/StandAlone/sus
INP=$UPATH/src/StandAlone/inputs/PhaseField/heat/be_sstruct_conv

for LD in $LDs; do
echo $LD
	for UPS in $INP/be_sstruct_conv_%(solver)s_3d_nlvl*_%(ph)d_*.ups; do
		LOG=$(basename $UPS).log
		[[ -f $LOG && -s $LOG ]] || mpirun -np $(np)d $SUS $UPS &> $LOG
	done
done
"""

NAME = "be_sstruct_conv_%(solver)s_3d_np%(np)d"

for solver in ("fac","split","gmres","flexgmres","lgmres","bicgstab"):
	for nlvl in range(1,9):
		for ph in range(0,5):
			for pk in range(0,5):
				np = 2**ph
				name = NAME % { "solver": solver, "np": np };
				arc = open(name+".arc4", "w")
				arc.write(ARC % { "np": np, "solver": solver, "ph": ph })
				arc.close()
