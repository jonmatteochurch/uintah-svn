#!/usr/bin/env python3

from os                  import environ
from sys                 import argv, exit
from helpers.runSusTests import getInputsDir, runSusTests
from helpers.modUPS      import modUPS, modUPS2

the_dir = "%s/%s" % ( getInputsDir(), "PhaseField" )

#______________________________________________________________________
#Performance and other UCF related tests

#  Test syntax: ( "folder name", "input file", # processors, "OS", ["flags1","flag2"])
#  flags:
#       gpu:                    - run test if machine is gpu enabled
#       no_uda_comparison:      - skip the uda comparisons
#       no_memoryTest:          - skip all memory checks
#       no_restart:             - skip the restart tests
#       no_dbg:                 - skip all debug compilation tests
#       no_opt:                 - skip all optimized compilation tests
#       no_cuda:                - skip test if this is a cuda enable build
#       do_performance_test:    - Run the performance test, log and plot simulation runtime.
#                                 (You cannot perform uda comparsions with this flag set)
#       doesTestRun:            - Checks if a test successfully runs
#       abs_tolerance=[double]  - absolute tolerance used in comparisons
#       rel_tolerance=[double]  - relative tolerance used in comparisons
#       exactComparison         - set absolute/relative tolerance = 0  for uda comparisons
#       postProcessRun          - start test from an existing uda in the checkpoints directory.  Compute new quantities and save them in a new uda
#       startFromCheckpoint     - start test from checkpoint. (/home/rt/CheckPoints/..../testname.uda.000)
#       sus_options="string"    - Additional command line options for sus command
#
#  Notes:
#  1) The "folder name" must be the same as input file without the extension.
#  2) Performance_tests are not run on a debug build.
#______________________________________________________________________

benchmark01_cc_ups = modUPS2 ( the_dir, "benchmark01/benchmark01_cc_eps020_n063_k3e-02.ups", [ \
  ("delete", "/Uintah_specification/DataArchiver/outputInterval" ), \
  ("update", "/Uintah_specification/Time/maxTime: 100 "), \
  ("append", "/Uintah_specification/DataArchiver/filebase:elem:outputTimestepInterval: 100 "), \
  ("append", "/Uintah_specification/Time/initTime:elem:max_Timesteps: 2001" ), \
])
benchmark01_nc_ups = modUPS2 ( the_dir, "benchmark01/benchmark01_nc_eps020_n064_k3e-02.ups", [ \
  ("delete", "/Uintah_specification/DataArchiver/outputInterval" ), \
  ("update", "/Uintah_specification/Time/maxTime: 100 "), \
  ("append", "/Uintah_specification/DataArchiver/filebase:elem:outputTimestepInterval: 100 "), \
  ("append", "/Uintah_specification/Time/initTime:elem:max_Timesteps: 2001" ), \
])
benchmark02_cc_ups = modUPS2 ( the_dir, "benchmark02/benchmark02_cc_eps010_n062_k1e-04.ups", [ \
  ("delete", "/Uintah_specification/DataArchiver/outputInterval" ), \
  ("append", "/Uintah_specification/DataArchiver/filebase:elem:outputTimestepInterval: 1000 "), \
  ("append", "/Uintah_specification/Time/initTime:elem:max_Timesteps: 20001" ), \
])
benchmark02_nc_ups = modUPS2 ( the_dir, "benchmark02/benchmark02_nc_eps010_n064_k1e-04.ups", [ \
  ("delete", "/Uintah_specification/DataArchiver/outputInterval" ), \
  ("append", "/Uintah_specification/DataArchiver/filebase:elem:outputTimestepInterval: 1000 "), \
  ("append", "/Uintah_specification/Time/initTime:elem:max_Timesteps: 20001" ), \
])
benchmark03_cc_ups = modUPS2 ( the_dir, "benchmark03/benchmark03_cc_n063_k3e-04.ups", [ \
  ("delete", "/Uintah_specification/DataArchiver/outputInterval" ), \
  ("append", "/Uintah_specification/DataArchiver/filebase:elem:outputTimestepInterval: 7500 "), \
  ("append", "/Uintah_specification/Time/initTime:elem:max_Timesteps: 150001" ), \
])
benchmark03_nc_ups = modUPS2 ( the_dir, "benchmark03/benchmark03_nc_n064_k3e-04.ups", [ \
  ("delete", "/Uintah_specification/DataArchiver/outputInterval" ), \
  ("append", "/Uintah_specification/DataArchiver/filebase:elem:outputTimestepInterval: 7500 "), \
  ("append", "/Uintah_specification/Time/initTime:elem:max_Timesteps: 150001" ), \
])
benchmark04_ups    = modUPS2( the_dir, "benchmark04/benchmark04_cc_n096.ups", [ \
  ("delete", "/Uintah_specification/DataArchiver/outputInterval" ), \
  ("append", "/Uintah_specification/DataArchiver/filebase:elem:outputTimestepInterval: 1500 "), \
  ("append", "/Uintah_specification/Time/initTime:elem:max_Timesteps: 30001" ), \
])

BENCHTEST = [
  ( "benchmark01_cc",                                   benchmark01_cc_ups,                                          4,   "All",   ["exactComparison"] ),
  ( "benchmark01_nc",                                   benchmark01_nc_ups,                                          4,   "All",   ["exactComparison"] ),
  ( "benchmark02_cc",                                   benchmark02_cc_ups,                                          4,   "All",   ["exactComparison"] ),
  ( "benchmark02_nc",                                   benchmark02_nc_ups,                                          4,   "All",   ["exactComparison"] ),
  ( "benchmark03_cc",                                   benchmark03_cc_ups,                                          2,   "All",   ["exactComparison"] ),
  ( "benchmark03_nc",                                   benchmark03_nc_ups,                                          2,   "All",   ["exactComparison"] ),
  ( "benchmark04",                                      benchmark04_ups,                                             4,   "All",   ["exactComparison"] ),
]

VAR = ["cc", "nc"]
DIM = ["2d", "3d"]
F2C = {"2d": ["fc0","fc1","fcsimple","fclinear","fcbilinear"], "3d": ["fc0","fc1"]}
FCN = ["fc0new","fc1new"]
IM = ["be", "cn"]
CC = ["cc"]

HEATTEST        = []
HEATNC3DTEST    = [] # FIXME Caught exception: ArrayIndexOutOfBoundsException
HEATMPITEST     = []
HEATNC3DMPITEST = [] # FIXME Caught exception: ArrayIndexOutOfBoundsException

freq = {"2d": 50, "3d": 25}
for dim in DIM:
  for var in VAR:
    ups = modUPS2 ( the_dir, "heat/heat_periodic_%s_%s_fe.ups" % (var,dim), [ \
      ("update", "/Uintah_specification/DataArchiver/outputTimestepInterval: %d " % (freq[dim]) ), \
      ("update", "/Uintah_specification/Time/maxTime: %d " % 100 ), \
      ("append", "/Uintah_specification/Time/initTime:elem:max_Timesteps: %d " % (20*freq[dim]+1) ), \
    ])
    if var=="nc" and dim =="3d":
      HEATNC3DTEST   .append( ( "heat_periodic_%s_%s_fe"     % (var,dim), ups, 1, "All", ["exactComparison"] ) )
      HEATNC3DMPITEST.append( ( "heat_periodic_%s_%s_fe_mpi" % (var,dim), ups, 4, "All", ["exactComparison"] ) )
    else:
      HEATTEST       .append( ( "heat_periodic_%s_%s_fe"     % (var,dim), ups, 1, "All", ["exactComparison"] ) )
      HEATMPITEST    .append( ( "heat_periodic_%s_%s_fe_mpi" % (var,dim), ups, 4, "All", ["exactComparison"] ) )

HEATBCTEST    = []
HEATBCMPITEST = []

freq = {"2d": 50, "3d": 20}
for dim in DIM:
  for var in VAR:
    ups = modUPS2 ( the_dir, "heat/heat_test_%s_%s_fe.ups" % (var,dim), [ \
      ("update", "/Uintah_specification/DataArchiver/outputTimestepInterval: %d " % (freq[dim]) ), \
      ("append", "/Uintah_specification/Time/initTime:elem:max_Timesteps: %d " % (20*freq[dim]+1) ), \
    ])
    HEATBCTEST   .append( ( "heat_test_%s_%s_fe"     % (var,dim), ups, 1, "All", ["exactComparison"] ) )
    HEATBCMPITEST.append( ( "heat_test_%s_%s_fe_mpi" % (var,dim), ups, 4, "All", ["exactComparison"] ) )

HEATAMRTEST    = [] # FIXME Caught exception: ArrayIndexOutOfBoundsException
HEATAMRMPITEST = [] # FIXME Caught exception: ArrayIndexOutOfBoundsException

freq = {"2d": 50, "3d": 20}
for dim in DIM:
  for f2c in F2C[dim]:
    for var in VAR:
      ups = modUPS2 ( the_dir, "heat/heat_periodic_%s_%s_fe_amr_%s.ups" % (var,dim,f2c), [ \
        ("update", "/Uintah_specification/DataArchiver/outputTimestepInterval: %d " % (freq[dim]) ), \
        ("append", "/Uintah_specification/Time/initTime:elem:max_Timesteps: %d " % (20*freq[dim]+1) ), \
      ])
      HEATAMRTEST   .append( ( "heat_periodic_%s_%s_fe_amr_%s"     % (var,dim,f2c), ups, 1, "All", ["exactComparison"] ) )
      HEATAMRMPITEST.append( ( "heat_periodic_%s_%s_fe_amr_%s_mpi" % (var,dim,f2c), ups, 4, "All", ["exactComparison"] ) )

HEATAMRBCTEST    = []
HEATAMRBCMPITEST = [] # FIXME Caught exception: Unknown variable: subproblems

for dim in DIM:
  for var in VAR:
    for f2c in F2C[dim]:
      ups = "heat/heat_test_%s_%s_fe_amr_%s.ups" % (var,dim,f2c)
      ups = modUPS2 ( the_dir, "heat/heat_periodic_%s_%s_fe_amr_%s.ups" % (var,dim,f2c), [ \
        ("update", "/Uintah_specification/DataArchiver/outputTimestepInterval: %d " % (freq[dim]) ), \
        ("append", "/Uintah_specification/Time/initTime:elem:max_Timesteps: %d " % (20*freq[dim]+1) ), \
      ])
      HEATAMRBCTEST   .append( ( "heat_test_%s_%s_fe_amr_%s"     % (var,dim,f2c), ups, 1, "All", ["exactComparison"] ) )
      HEATAMRBCMPITEST.append( ( "heat_test_%s_%s_fe_amr_%s_mpi" % (var,dim,f2c), ups, 4, "All", ["exactComparison"] ) )

HEATHYPREBETEST    = [] # FIXME Caught exception: Unknown variable: hypre_solver_label
HEATHYPRECNTEST    = [] # FIXME Caught exception: Unknown variable: A
HEATHYPRECNMPITEST = [] # FIXME Caught exception: Unknown variable: A
HEATHYPREBEMPITEST = [] # FIXME Caught exception: Unknown variable: hypre_solver_label

for dim in DIM:
  # for var in VAR:
    for var in CC:
      for im in IM:
        ups = "heat/heat_periodic_%s_%s_%s.ups" % (var,dim,im)
        if   im=="be":
          HEATHYPREBETEST   .append( ( "heat_periodic_%s_%s_%s"     % (var,dim,im), ups, 1, "All", ["exactComparison"] ) )
          HEATHYPREBEMPITEST.append( ( "heat_periodic_%s_%s_%s_mpi" % (var,dim,im), ups, 4, "All", ["exactComparison"] ) )
        elif im=="cn":
          HEATHYPREBETEST   .append( ( "heat_periodic_%s_%s_%s"     % (var,dim,im), ups, 1, "All", ["exactComparison"] ) )
          HEATHYPREBEMPITEST.append( ( "heat_periodic_%s_%s_%s_mpi" % (var,dim,im), ups, 4, "All", ["exactComparison"] ) )

HEATHYPREAMRTEST    = [] # FIXME Caught exception: ArrayIndexOutOfBoundsException
HEATHYPREAMRMPITEST = [] # FIXME Caught exception: ArrayIndexOutOfBoundsException

for dim in DIM:
  # for var in VAR:
    for var in CC:
      for im in IM:
        for f2c in FCN:
          ups = "heat/heat_periodic_%s_%s_%s_amr_hypre_%s.ups" % (var,dim,im,f2c)
          HEATHYPREAMRTEST   .append( ( "heat_periodic_%s_%s_%s_amr_hypre_%s"     % (var,dim,im,f2c), ups, 1, "All", ["exactComparison"] ) )
          HEATHYPREAMRMPITEST.append( ( "heat_periodic_%s_%s_%s_amr_hypre_%s_mpi" % (var,dim,im,f2c), ups, 4, "All", ["exactComparison"] ) )

HEATHYPREAMRBC2DBETEST = [] # FIXME Caught exception: Unknown variable: hypre_solver_label
HEATHYPREAMRBC2DCNTEST = [] # FIXME Caught exception: Unknown variable: A
HEATHYPREAMRBC3DTEST   = [] # FIXME Caught exception: An IntVectorOutOfBounds exception was thrown.
HEATHYPREAMRBCMPITEST  = [] # FIXME Caught exception: Unknown variable: subproblems

for dim in DIM:
  # for var in VAR:
    for var in CC:
      for im in IM:
        for f2c in FCN:
          ups = "heat/heat_test_%s_%s_%s_amr_hypre_%s.ups" % (var,dim,im,f2c)
          if   dim=="2d" and im=="be":
            HEATHYPREAMRBC2DBETEST.append( ( "heat_test_%s_%s_%s_amr_hypre_%s"     % (var,dim,im,f2c), ups, 1, "All", ["exactComparison"] ) )
          elif dim=="2d" and im=="cn":
            HEATHYPREAMRBC2DCNTEST.append( ( "heat_test_%s_%s_%s_amr_hypre_%s"     % (var,dim,im,f2c), ups, 1, "All", ["exactComparison"] ) )
          elif dim=="3d":
            HEATHYPREAMRBC3DTEST .append( ( "heat_test_%s_%s_%s_amr_hypre_%s"     % (var,dim,im,f2c), ups, 1, "All", ["exactComparison"] ) )
          HEATHYPREAMRBCMPITEST  .append( ( "heat_test_%s_%s_%s_amr_hypre_%s_mpi" % (var,dim,im,f2c), ups, 4, "All", ["exactComparison"] ) )

HEATHYPREFACTEST    = [] # TODO Error: Invalid string value for Attribute: Uintah_specification->Solver->type. 'hyprefac' not found in this list: CGSolver, hypre, hypreamr
HEATHYPREFACMPITEST = [] # TODO Error: Invalid string value for Attribute: Uintah_specification->Solver->type. 'hyprefac' not found in this list: CGSolver, hypre, hypreamr

for dim in DIM:
  # for var in VAR:
    for var in CC:
      for im in IM:
        for f2c in F2C[dim]:
          ups = "heat/heat_periodic_%s_%s_%s_amr_hyprefac_%s.ups" % (var,dim,im,f2c)
          HEATHYPREFACTEST   .append( ( "heat_periodic_%s_%s_%s_amr_hyprefac_%s"     % (var,dim,im,f2c), ups, 1, "All", ["exactComparison"] ) )
          HEATHYPREFACMPITEST.append( ( "heat_periodic_%s_%s_%s_amr_hyprefac_%s_mpi" % (var,dim,im,f2c), ups, 4, "All", ["exactComparison"] ) )

HEATHYPREFACBCTEST    = [] # TODO Error: Invalid string value for Attribute: Uintah_specification->Solver->type. 'hyprefac' not found in this list: CGSolver, hypre, hypreamr
HEATHYPREFACBCMPITEST = [] # TODO Error: Invalid string value for Attribute: Uintah_specification->Solver->type. 'hyprefac' not found in this list: CGSolver, hypre, hypreamr

for dim in DIM:
  # for var in VAR:
    for var in CC:
      for im in IM:
        for f2c in F2C[dim]:
          ups = "heat/heat_test_%s_%s_%s_amr_hyprefac_%s.ups" % (var,dim,im,f2c)
          HEATHYPREFACBCTEST   .append( ( "heat_test_%s_%s_%s_amr_hyprefac_%s"     % (var,dim,im,f2c), ups, 1, "All", ["exactComparison"] ) )
          HEATHYPREFACBCMPITEST.append( ( "heat_test_%s_%s_%s_amr_hyprefac_%s_mpi" % (var,dim,im,f2c), ups, 4, "All", ["exactComparison"] ) )

PUREMETALTEST    = [] # FIXME Caught exception: ArrayIndexOutOfBoundsException
PUREMETALMPITEST = [] # FIXME Caught exception: ArrayIndexOutOfBoundsException

freq = {"2d": 25, "3d": 1}
for dim in DIM:
  for var in VAR:
    ups = modUPS2 ( the_dir, "pure_metal/pure_metal_%s_%s.ups" % (var,dim), [ \
      ("update", "/Uintah_specification/DataArchiver/outputTimestepInterval: %d " % (freq[dim]) ), \
      ("append", "/Uintah_specification/Time/initTime:elem:max_Timesteps: %d " % (20*freq[dim]+1) ), \
    ])
    PUREMETALTEST   .append( ( "pure_metal_%s_%s"     % (var,dim), ups, 1, "All", ["exactComparison"] ) )
    PUREMETALMPITEST.append( ( "pure_metal_%s_%s_mpi" % (var,dim), ups, 4, "All", ["exactComparison"] ) )

PUREMETALAMRTEST    = [] # FIXME Caught exception: ArrayIndexOutOfBoundsException
PUREMETALAMRMPITEST = [] # FIXME Caught exception: ArrayIndexOutOfBoundsException

for dim in DIM:
  for var in VAR:
    for f2c in F2C[dim]:
      ups = "pure_metal/pure_metal_%s_%s_amr_%s.ups" % (var,dim,f2c)
      PUREMETALAMRTEST   .append( ( "pure_metal_%s_%s_amr_%s"     % (var,dim,f2c), ups, 1, "All", ["exactComparison"] ) )
      PUREMETALAMRMPITEST.append( ( "pure_metal_%s_%s_amr_%s_mpi" % (var,dim,f2c), ups, 4, "All", ["exactComparison"] ) )

#__________________________________
# The following list is parsed by the local RT script
# and allows the user to select the tests to run
#LIST: BENCHTEST HEATTEST HEATNC3DTEST HEATMPITEST HEATNC3DMPITEST HEATBCTEST HEATBCMPITEST HEATAMRTEST HEATAMRMPITEST HEATAMRBCTEST HEATAMRBCMPITEST HEATHYPREBETEST HEATHYPRECNTEST HEATHYPRECNMPITEST HEATHYPREBEMPITEST HEATHYPREAMRTEST HEATHYPREAMRMPITEST HEATHYPREAMRBC2DBETEST HEATHYPREAMRBC2DCNTEST HEATHYPREAMRBC3DTEST HEATHYPREAMRBCMPITEST HEATHYPREFACTEST HEATHYPREFACMPITEST HEATHYPREFACBCTEST HEATHYPREFACBCMPITEST PUREMETALTEST PUREMETALMPITEST PUREMETALAMRTEST PUREMETALAMRMPITEST
#___________________________________

# returns the list
def getTestList(me) :
  if me == "BENCHTEST":
    TESTS = BENCHTEST

  elif me == "HEATTEST":
    TESTS = HEATTEST
  elif me == "HEATNC3DTEST":
    TESTS = HEATNC3DTEST
  elif me == "HEATMPITEST":
    TESTS = HEATMPITEST
  elif me == "HEATNC3DMPITEST":
    TESTS = HEATNC3DMPITEST

  elif me == "HEATBCTEST":
    TESTS = HEATBCTEST
  elif me == "HEATBCMPITEST":
    TESTS = HEATBCMPITEST

  elif me == "HEATAMRTEST":
    TESTS = HEATAMRTEST
  elif me == "HEATAMRMPITEST":
    TESTS = HEATAMRMPITEST

  elif me == "HEATAMRBCTEST":
    TESTS = HEATAMRBCTEST
  elif me == "HEATAMRBCMPITEST":
    TESTS = HEATAMRBCMPITEST

  elif me == "HEATHYPREBETEST":
    TESTS = HEATHYPREBETEST
  elif me == "HEATHYPRECNTEST":
    TESTS = HEATHYPRECNTEST
  elif me == "HEATHYPREBEMPITEST":
    TESTS = HEATHYPREBEMPITEST
  elif me == "HEATHYPRECNMPITEST":
    TESTS = HEATHYPRECNMPITEST

  elif me == "HEATHYPREAMRTEST":
    TESTS = HEATHYPREAMRTEST
  elif me == "HEATHYPREAMRMPITEST":
    TESTS = HEATHYPREAMRMPITEST

  elif me == "HEATHYPREAMRBC2DBETEST":
    TESTS = HEATHYPREAMRBC2DBETEST
  elif me == "HEATHYPREAMRBC2DCNTEST":
    TESTS = HEATHYPREAMRBC2DCNTEST
  elif me == "HEATHYPREAMRBC3DTEST":
    TESTS = HEATHYPREAMRBC3DTEST
  elif me == "HEATHYPREAMRBCMPITEST":
    TESTS = HEATHYPREAMRBCMPITEST

  elif me == "HEATHYPREFACTEST":
    TESTS = HEATHYPREFACTEST
  elif me == "HEATHYPREFACMPITEST":
    TESTS = HEATHYPREFACMPITEST

  elif me == "HEATHYPREFACBCTEST":
    TESTS = HEATHYPREFACBCTEST
  elif me == "HEATHYPREFACBCMPITEST":
    TESTS = HEATHYPREFACBCMPITEST

  elif me == "PUREMETALTEST":
    TESTS = PUREMETALTEST
  elif me == "PUREMETALMPITEST":
    TESTS = PUREMETALMPITEST
    
  elif me == "PUREMETALAMRTEST":
    TESTS = PUREMETALAMRTEST
  elif me == "PUREMETALAMRMPITEST":
    TESTS = PUREMETALAMRMPITEST

  elif me == "AMRTESTS":
    TESTS = HEATAMRTEST + HEATAMRMPITEST + HEATAMRBCTEST + HEATAMRBCMPITEST + HEATHYPREAMRTEST + HEATHYPREAMRMPITEST + HEATHYPREAMRBC2DBETEST + HEATHYPREAMRBC2DCNTEST + HEATHYPREAMRBC3DTEST + HEATHYPREAMRBCMPITEST + PUREMETALAMRTEST + PUREMETALAMRMPITEST

  elif me == "DEBUGTESTS":
    TESTS = HEATNC3DTEST + HEATNC3DMPITEST + HEATAMRTEST + HEATAMRMPITEST + HEATAMRBCMPITEST + HEATHYPREBETEST + HEATHYPRECNTEST + HEATHYPRECNMPITEST + HEATHYPREBEMPITEST + HEATHYPREAMRTEST + HEATHYPREAMRMPITEST + HEATHYPREAMRBC2DBETEST + HEATHYPREAMRBC2DCNTEST + HEATHYPREAMRBC3DTEST + HEATHYPREAMRBCMPITEST + HEATHYPREFACTEST + HEATHYPREFACMPITEST + HEATHYPREFACBCTEST + HEATHYPREFACBCMPITEST + PUREMETALTEST + PUREMETALMPITEST + PUREMETALAMRTEST + PUREMETALAMRMPITEST

  elif me == "LOCALTESTS":
    TESTS = BENCHTEST + HEATMPITEST + HEATBCMPITEST + HEATAMRBCTEST
  elif me == "NIGHTLYTESTS":
    TESTS = BENCHTEST + HEATMPITEST + HEATBCMPITEST + HEATAMRBCTEST
  elif me == "BUILDBOTTESTS":
    TESTS = BENCHTEST + HEATMPITEST + HEATBCMPITEST + HEATAMRBCTEST

  else:
    print("\nERROR:PhaseField.py  getTestList:  The test list (%s) does not exist!\n\n" % me)
    exit(1)
  return TESTS
#__________________________________

if __name__ == "__main__":

  TESTS = getTestList( environ['WHICH_TESTS'] )

  result = runSusTests(argv, TESTS, "PhaseField")
  exit( result )
