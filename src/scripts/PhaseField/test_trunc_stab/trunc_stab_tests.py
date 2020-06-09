#!/usr/bin/env python

from sys      import argv, path, exit
from os       import rename, mkdir, environ
from os.path  import exists

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
scr=src + "/scripts/PhaseField/test_trunc_stab"
sge=scr + "/sge"

if not exists(inp + "/test_trunc_stab"):
  mkdir (inp + "/test_trunc_stab")

if not exists(scr):
  mkdir (scr)

if exists(sge):
  rmtree (sge)
mkdir (sge)

ppp=8  # patches per process
cpn=40 # cores per node

# 2D CC
for ph in range (0,5):
  for pk in range (0,7):
    h = 2**(-ph)
    k = 2**(-pk)
    p = 2**(ph+1)
    i = 2*(pk)
    nam = "heat_test_cc_2d_fe_%d_%d" % (ph,pk)
    ups = modUPS2 ( inp, "heat/heat_test_cc_2d_fe.ups", [ \
      ("update", "/Uintah_specification/Meta/title:%s" % (nam) ), \
      ("update", "/Uintah_specification/PhaseField/delt:%f" % (k) ), \
      ("update", "/Uintah_specification/Grid/Level/spacing:[%f,%f,1.]" % (h,h) ), \
      ("update", "/Uintah_specification/Grid/Level/Box/patches:[%d,%d,1]" % (p,p) ), \
      ("update", "/Uintah_specification/DataArchiver/outputTimestepInterval:%d " % (i) ), \
      ("update", "/Uintah_specification/DataArchiver/filebase:%s.uda " % (nam) ), \
    ])
    rename( inp+"/"+ups, inp + "/test_trunc_stab/" + nam + ".ups" )
    np = (p*p+ppp-1)//ppp
    nodes = (np+cpn-1)//cpn
    sn = sge + "/test_trunc_stab_N%03d.arc4" % ( nodes )
    if not exists(sn):
      sf = open(sn,"w+")
      sf.write("#$ -cwd -V\n")
      sf.write("#$ -m be -M scjmc@leeds.ac.uk\n")
      sf.write("#$ -l nodes=%d\n" % nodes)
      sf.write("#$ -l h_rt=48:00:00\n")
      sf.write("\n")
      sf.write("if [ -z \"$UPATH\" ]; then\n")
      sf.write("    UPATH=$HOME/uintah/trunk\n")
      sf.write("fi\n")
      sf.write("\n")
      sf.write("SUS=$UPATH/arc4/build/StandAlone/sus\n")
      sf.write("INP=%s/test_trunc_stab\n" % (inp))
      sf.write("MPIRUN=/apps/developers/compilers/intel/19.0.4/1/default/impi/2019.4.243/intel64/bin/mpirun\n")
      sf.write("\n")
    else:
      sf = open(sn,"a+")

    sf.write("FN=%s\n" % (nam))
    sf.write("$MPIRUN -np %d $SUS $INP/$FN.ups > $FN.log 2> $FN.err\n" % (np))
    sf.write("\n")

# 2D NC
for ph in range (0,5):
  for pk in range (0,7):
    h = 2**(-ph)
    k = 2**(-pk)
    p = 2**(ph+1)
    i = 2*(pk)
    nam = "heat_test_nc_2d_fe_%d_%d" % (ph,pk)
    ups = modUPS2 ( inp, "heat/heat_test_nc_2d_fe.ups", [ \
      ("update", "/Uintah_specification/Meta/title:%s" % (nam) ), \
      ("update", "/Uintah_specification/PhaseField/delt:%f" % (k) ), \
      ("update", "/Uintah_specification/Grid/Level/spacing:[%f,%f,1.]" % (h,h) ), \
      ("update", "/Uintah_specification/Grid/Level/Box/patches:[%d,%d,1]" % (p,p) ), \
      ("update", "/Uintah_specification/DataArchiver/outputTimestepInterval:%d " % (i) ), \
      ("update", "/Uintah_specification/DataArchiver/filebase:%s.uda " % (nam) ), \
    ])
    rename( inp+"/"+ups, inp + "/test_trunc_stab/" + nam + ".ups" )
    np = (p*p+ppp-1)//ppp
    nodes = (np+cpn-1)//cpn
    sn = sge + "/test_trunc_stab_N%03d.arc4" % ( nodes )
    if not exists(sn):
      sf = open(sn,"w+")
      sf.write("#$ -cwd -V\n")
      sf.write("#$ -m be -M scjmc@leeds.ac.uk\n")
      sf.write("#$ -l nodes=%d\n" % nodes)
      sf.write("#$ -l h_rt=48:00:00\n")
      sf.write("\n")
      sf.write("if [ -z \"$UPATH\" ]; then\n")
      sf.write("    UPATH=$HOME/uintah/trunk\n")
      sf.write("fi\n")
      sf.write("\n")
      sf.write("SUS=$UPATH/arc4/build/StandAlone/sus\n")
      sf.write("INP=%s/test_trunc_stab\n" % (inp))
      sf.write("MPIRUN=/apps/developers/compilers/intel/19.0.4/1/default/impi/2019.4.243/intel64/bin/mpirun\n")
      sf.write("\n")
    else:
      sf = open(sn,"a+")

    sf.write("FN=%s\n" % (nam))
    sf.write("$MPIRUN -np %d $SUS $INP/$FN.ups > $FN.log 2> $FN.err\n" % (np))
    sf.write("\n")

# 3D CC
for ph in range (0,5):
  for pk in range (0,7):
    h = 2**(-ph)
    k = 2**(-pk)
    p = 2**(ph+1)
    i = 2*(pk)
    nam = "heat_test_cc_3d_fe_%d_%d" % (ph,pk)
    ups = modUPS2 ( inp, "heat/heat_test_cc_3d_fe.ups", [ \
      ("update", "/Uintah_specification/Meta/title:%s" % (nam) ), \
      ("update", "/Uintah_specification/PhaseField/delt:%f" % (k) ), \
      ("update", "/Uintah_specification/Grid/Level/spacing:[%f,%f,%f]" % (h,h,h) ), \
      ("update", "/Uintah_specification/Grid/Level/Box/patches:[%d,%d,%d]" % (p,p,p) ), \
      ("update", "/Uintah_specification/DataArchiver/outputTimestepInterval:%d " % (i) ), \
      ("update", "/Uintah_specification/DataArchiver/filebase:%s.uda " % (nam) ), \
    ])
    rename( inp+"/"+ups, inp + "/test_trunc_stab/" + nam + ".ups" )
    np = (p*p*p+ppp-1)//ppp
    nodes = (np+cpn-1)//cpn
    sn = sge + "/test_trunc_stab_N%03d.arc4" % ( nodes )
    if not exists(sn):
      sf = open(sn,"w+")
      sf.write("#$ -cwd -V\n")
      sf.write("#$ -m be -M scjmc@leeds.ac.uk\n")
      sf.write("#$ -l nodes=%d\n" % nodes)
      sf.write("#$ -l h_rt=48:00:00\n")
      sf.write("\n")
      sf.write("if [ -z \"$UPATH\" ]; then\n")
      sf.write("    UPATH=$HOME/uintah/trunk\n")
      sf.write("fi\n")
      sf.write("\n")
      sf.write("SUS=$UPATH/arc4/build/StandAlone/sus\n")
      sf.write("INP=%s/test_trunc_stab\n" % (inp))
      sf.write("MPIRUN=/apps/developers/compilers/intel/19.0.4/1/default/impi/2019.4.243/intel64/bin/mpirun\n")
      sf.write("\n")
    else:
      sf = open(sn,"a+")

    sf.write("FN=%s\n" % (nam))
    sf.write("$MPIRUN -np %d $SUS $INP/$FN.ups > $FN.log 2> $FN.err\n" % (np))
    sf.write("\n")

# 3D NC
for ph in range (0,5):
  for pk in range (0,7):
    h = 2**(-ph)
    k = 2**(-pk)
    p = 2**(ph+1)
    i = 2*(pk)
    nam = "heat_test_nc_3d_fe_%d_%d" % (ph,pk)
    ups = modUPS2 ( inp, "heat/heat_test_nc_3d_fe.ups", [ \
      ("update", "/Uintah_specification/Meta/title:%s" % (nam) ), \
      ("update", "/Uintah_specification/PhaseField/delt:%f" % (k) ), \
      ("update", "/Uintah_specification/Grid/Level/spacing:[%f,%f,%f]" % (h,h,h) ), \
      ("update", "/Uintah_specification/Grid/Level/Box/patches:[%d,%d,%d]" % (p,p,p) ), \
      ("update", "/Uintah_specification/DataArchiver/outputTimestepInterval:%d " % (i) ), \
      ("update", "/Uintah_specification/DataArchiver/filebase:%s.uda " % (nam) ), \
    ])
    rename( inp+"/"+ups, inp + "/test_trunc_stab/" + nam + ".ups" )
    np = (p*p*p+ppp-1)//ppp
    nodes = (np+cpn-1)//cpn
    sn = sge + "/test_trunc_stab_N%03d.arc4" % ( nodes )
    if not exists(sn):
      sf = open(sn,"w+")
      sf.write("#$ -cwd -V\n")
      sf.write("#$ -m be -M scjmc@leeds.ac.uk\n")
      sf.write("#$ -l nodes=%d\n" % nodes)
      sf.write("#$ -l h_rt=48:00:00\n")
      sf.write("\n")
      sf.write("if [ -z \"$UPATH\" ]; then\n")
      sf.write("    UPATH=$HOME/uintah/trunk\n")
      sf.write("fi\n")
      sf.write("\n")
      sf.write("SUS=$UPATH/arc4/build/StandAlone/sus\n")
      sf.write("INP=%s/test_trunc_stab\n" % (inp))
      sf.write("MPIRUN=/apps/developers/compilers/intel/19.0.4/1/default/impi/2019.4.243/intel64/bin/mpirun\n")
      sf.write("\n")
    else:
      sf = open(sn,"a+")

    sf.write("FN=%s\n" % (nam))
    sf.write("$MPIRUN -np %d $SUS $INP/$FN.ups > $FN.log 2> $FN.err\n" % (np))
    sf.write("\n")
