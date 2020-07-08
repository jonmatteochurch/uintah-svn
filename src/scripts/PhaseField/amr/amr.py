#!/usr/bin/env python

from sys      import argv, path, exit
from os       import mkdir, environ
from os.path  import exists
from shutil   import rmtree, move

if len(argv)<2 or len(argv)>4:
  exit ("usage: %s <uintah src path> [<mpirun bin path> <xmlstarlet bin path>]" % (__file__))

src=argv[1]
mpirun="mpirun"

if len(argv)>2:
  mpirun=argv[2] + "/mpirun"

if len(argv)>3:
  environ["PATH"] += ":" + argv[3]

path.append(src + "/R_Tester")

from helpers.modUPS      import modUPS2

inp=src + "/StandAlone/inputs/PhaseField"
scr=src + "/scripts/PhaseField/amr"
sge=scr + "/sge"

if not exists(inp + "/pure_metal/amr"):
  mkdir (inp + "/pure_metal/amr")

if not exists(scr):
  mkdir (scr)

if exists(sge):
  rmtree (sge)
mkdir (sge)

for l in range(2,3): #range (1,9):
  for s in range(0,1): #range (0,4):
    for r in range(0,1): #range (0,4):
      for b in range(0,1): #range (0,4):
        for ph in range (1,6):
          for f in ("amr_parallel","amr_diagonal"):
              h = 0.4*2**(l-1)
              p = 2**(ph)
              nam = "%s_l%d_s%d_r%d_b%d_p%d" % (f,l,s,r,b,p)
              ups = modUPS2 ( "/home/jonmatteochurch/Developer/uintah/trunk/out/amr/input/", "%s.ups" % (f), [ \
                ("update", "/Uintah_specification/Meta/title:%s" % (nam) ), \
                ("update", "/Uintah_specification/Grid/Level/spacing:[%f,%f,1.]" % (h,h) ), \
                ("update", "/Uintah_specification/AMR/Regridder/max_levels:%d " % (l) ), \
                ("update", "/Uintah_specification/AMR/Regridder/cell_stability_dilation:[%d,%d,0] " % (s,s) ), \
                ("update", "/Uintah_specification/AMR/Regridder/cell_regrid_dilation:[%d,%d,0] " % (r,r) ), \
                ("update", "/Uintah_specification/AMR/Regridder/min_boundary_cells:[%d,%d,0] " % (b,b) ), \
                ("update", "/Uintah_specification/AMR/Regridder/min_patch_size:[%d,%d,1] " % (p,p) ), \
                ("update", "/Uintah_specification/DataArchiver/filebase:%s.uda " % (nam) ) \
              ])
              print("%s/%s" % ("/home/jonmatteochurch/Developer/uintah/trunk/out/amr/input/",ups))
              print("%s/amr/%s.ups" % (inp,nam))
              move("%s/%s" % ("/home/jonmatteochurch/Developer/uintah/trunk/out/amr/input/",ups), "%s/amr/%s.ups" % (inp,nam))
    
              sn = "%s/%s.arc4" % (sge,nam);
              sf = open(sn,"w+")
              sf.write("#$ -cwd -V\n")
              sf.write("#$ -m be -M scjmc@leeds.ac.uk\n")
              sf.write("#$ -pe ib 12")
              sf.write("#$ -l h_rt=5:00:00\n")
              sf.write("#$ -l h_vmem=3G")
              sf.write("\n")
              sf.write("if [ -z \"$UPATH\" ]; then\n")
              sf.write("    UPATH=$HOME/uintah/trunk\n")
              sf.write("fi\n")
              sf.write("\n")
              sf.write("SUS=$UPATH/arc/build/StandAlone/sus\n")
              sf.write("INP=%s/amr\n" % (inp))
              sf.write("\n")
              sf.write("FN=%s\n" % (nam))
              sf.write("mpirun $SUS $INP/$FN.ups > $FN.log 2> $FN.err\n")
              sf.write("\n")
