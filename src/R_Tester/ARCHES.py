#!/usr/bin/env python

from os import symlink,environ, system
from sys import argv,exit,platform
from helpers.runSusTests import runSusTests, inputs_root, generatingGoldStandards
from helpers.modUPS import modUPS

the_dir = generatingGoldStandards()

if the_dir == "" :
  the_dir = "%s/ARCHES" % inputs_root()
else :
  the_dir = the_dir + "/ARCHES"

#______________________________________________________________________                            
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
#       startFromCheckpoint     - start test from checkpoint. (/home/csafe-tester/CheckPoints/..../testname.uda.000)
#       sus_options="string"    - Additional command line options for sus command
#
#  Notes: 
#  1) The "folder name" must be the same as input file without the extension.
#  2) If the processors is > 1.0 then an mpirun command will be used
#  3) Performance_tests are not run on a debug build.
#______________________________________________________________________
#  Test: purpose
#  massource_con_den: a constant density cold channel flow with a constant mass source. used to test mass conservation for pressure solver (with mass source terms)
#  massource_var_den: a two-fluid mixing cold channel flow with a constant mass source, used to test mass conservation for pressure solver (with mass source terms)
#  massource_coal_DQMOM: a coal reacting channel flow with a constant mass source term (or with coal mass source term), used to test mass conservation for pressure solver
#  methane_RCCE: 3m fire using the westbrook dryer/RCCE model

NIGHTLYTESTS = [
   ("constantMMS"                       , "mms/constantMMS.ups"                                     , 1.1 , "All"   , ["exactComparison"])   , 
   ("almgrenMMS"                        , "mms/almgrenMMS.ups"                                      , 1.1 , "All"   , ["exactComparison"])   , 
   ("oned_pulse_conv"                   , "mms/oned_pulse_conv.ups"                                 , 1.1 , "All"   , ["exactComparison"])   , 
   ("isotropic-turbulence-decay"        , "periodicTurb/isotropic-turbulence-decay.ups"             , 1.1 , "All"   , ["exactComparison"     , "no_restart"]) , 
   ("helium_1m"                         , "helium_1m.ups"                                           , 1.1 , "All"   , ["exactComparison"])   , 
   ("methane_fire"                      , "methane_fire.ups"                                        , 1.1 , "All"   , ["exactComparison"])   , 
   ("methane_fire_8patch"               , "methane_fire_8patch.ups"                                 , 8   , "All"   , ["exactComparison"])   , 
   ("methane_fire_8patch_petscrad"      , "methane_fire_8patch_petscrad.ups"                        , 8   , "All"   , ["exactComparison"     , "no_cuda" ])   ,    #11/1/6 gpu machine doesn't have petsc 
   ("rmcrt_bm1_1L"                      , "RMCRT/rmcrt_bm1_1L.ups"                                  , 1.1 , "Linux" , ["exactComparison"])   , 
   ("rmcrt_bm1_1L_bounded"               , "RMCRT/rmcrt_bm1_1L_bounded.ups"                         , 8   , "ALL"   , ["exactComparison"])   ,
   ("rmcrt_bm1_DO"                      , "RMCRT/rmcrt_bm1_DO.ups"                                  , 8   , "Linux" , ["exactComparison"])   , 
   ("rmcrt_bm1_ML"                      , "RMCRT/rmcrt_bm1_ML.ups"                                  , 1.1 , "Linux" , ["exactComparison"])   , 
   ("methane_rmcrt"                     , "RMCRT/methane_rmcrt.ups"                                 , 8   , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("dqmom_test_1"                      , "DQMOM_regression/dqmom_test_1.ups"                       , 1.1 , "Linux" , ["exactComparison"])   , 
   ("dqmom_test_2"                      , "DQMOM_regression/dqmom_test_2.ups"                       , 1.1 , "Linux" , ["exactComparison"])   , 
   ("dqmom_test_3"                      , "DQMOM_regression/dqmom_test_3.ups"                       , 1.1 , "Linux" , ["exactComparison"])   , 
   ("dqmom_test_4"                      , "DQMOM_regression/dqmom_test_4.ups"                       , 1.1 , "Linux" , ["exactComparison"])   , 
   ("dqmom_test_5"                      , "DQMOM_regression/dqmom_test_5.ups"                       , 1.1 , "Linux" , ["exactComparison"])   , 
   ("birth_test"                        , "DQMOM_regression/birth_test.ups"                         , 1.1 , "Linux" , ["exactComparison"])   , 
   ("upwind_birth_test"                 , "DQMOM_regression/upwind_birth_test.ups"                  , 1.1 , "Linux" , ["exactComparison"])   , 
   ("methane_jet"                       , "ClassicMixingTables/ups/methane_jet.ups"                 , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("acidbase_jet"                      , "ClassicMixingTables/ups/acidbase_jet.ups"                , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("acidbase_jet_2D"                   , "ClassicMixingTables/ups/acidbase_jet_2D.ups"             , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("1Dtabletest"                       , "ClassicMixingTables/ups/1DTableTest.ups"                 , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("xplus_scalar_test"                 , "ScalarTests/xplus_scalar_test.ups"                       , 6   , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("yplus_scalar_test"                 , "ScalarTests/yplus_scalar_test.ups"                       , 6   , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("zplus_scalar_test"                 , "ScalarTests/zplus_scalar_test.ups"                       , 6   , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("xminus_scalar_test"                , "ScalarTests/xminus_scalar_test.ups"                      , 6   , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("yminus_scalar_test"                , "ScalarTests/yminus_scalar_test.ups"                      , 6   , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("zminus_scalar_test"                , "ScalarTests/zminus_scalar_test.ups"                      , 6   , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("turbulent_inlet_test"              , "DigitalFilter/TurbulentInletChannel.ups"                 , 6   , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("masssource_con_den"                , "verify_masssource/source_channel_conden.ups"             , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("masssource_var_den"                , "verify_masssource/source_channel_varden.ups"             , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("heptane_pipe"                      , "heptane_pipe.ups"                                        , 1.1 , "Linux" , ["exactComparison"])   , 
   ("coal_table_pipe"                   , "coal_table_pipe.ups"                                     , 1.1 , "Linux" , ["exactComparison"])   , 
   ("scalar_var_1eqn"                   , "scalar_variance_1eqn.ups"                                , 4   , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("scalar_var_2eqn"                   , "scalar_variance_2eqn.ups"                                , 4   , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("coal_channel_FOWY"                 , "Coal/coal_channel_FOWY.ups"                              , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("pcoal_drag"                        , "Coal/pcoal_drag.ups"                                     , 1.1 , "Linux" , ["exactComparison"])   , 
   ("methane_RCCE"                      , "methane_RCCE.ups"                                        , 1.1 , "Linux" , ["exactComparison"])   , 
   ("channel_WD_CO"                     , "channel_WD_CO.ups"                                       , 1.1 , "Linux" , ["exactComparison"])   , 
   ("DOM16"                             , "DOM16.ups"                                               , 3   , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("DO_RadProps"                       , "DO_RadProps.ups"                                         , 3   , "Linux" , ["exactComparison"     , "no_cuda"   ])   ,   # 11/1/16 bug in radProps with cuda
   ("CQMOM_1x1"                         , "CQMOM_regression/CQMOM_1x1.ups"                          , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("CQMOM_scalar_transport"            , "CQMOM_regression/CQMOM_Transport.ups"                    , 6   , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("CQMOM_scalar_transport2x2x2"       , "CQMOM_regression/CQMOM_Transport_2x2x2.ups"              , 6   , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("SimpleBoxPTC"                      , "CQMOM_regression/SimpleBoxPTC.ups"                       , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("PTC_2D"                            , "CQMOM_regression/PTC_2D.ups"                             , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("PTC_3D"                            , "CQMOM_regression/PTC_3D.ups"                             , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("CQMOM_4D"                          , "CQMOM_regression/CQMOM_4D.ups"                           , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) ,
   ("CQMOM_7D"                          , "CQMOM_regression/CQMOM_7D.ups"                           , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) ,
   ("singleJet_poly"                    , "CQMOM_regression/singleJet_poly.ups"                     , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) ,
   ("angledWall"                        , "CQMOM_regression/angledWall.ups"                         , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) ,
   ("angledWall3D"                      , "CQMOM_regression/angledWall3D.ups"                       , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) ,
   ("Constant_Deposition"               , "CQMOM_regression/Constant_Deposition.ups"                , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) ,
   ("CQMOM_coal_test"                   , "CQMOM_regression/CQMOM_coal_test.ups"                    , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("channel_LagPart_inlet"             , "LagrangianParticles/channel_flow_x_lagrangian_inlet.ups" , 1.1 , "Linux" , ["exactComparison"     , "no_restart", "no_cuda"]) ,  # 11/1/16 bug with gpu support with particles  
   ("partMassFlow"                      , "Coal/partMassFlow.ups"                                   , 8   , "All"  , ["exactComparison"   ]) , 
   ("task_math"                         , "task_math.ups"                                           , 1.1 , "All"  , ["exactComparison", "no_restart"]) , 
   ("intrusion_test"                    , "intrusion_test.ups"                                      , 1.1 , "All"  , ["exactComparison"]) , 
   ("multi-patch-intrusion-test"        , "multi-patch-intrusion-test.ups"                          , 8   , "All"  , ["exactComparison"]) , 
   ("cloudBM24LS"                       , "cloudBM24LS.ups"                                         , 8 , "All"  , ["exactComparison", "no_restart"]) ,
   ("cloudBM48LS"                       , "cloudBM48LS.ups"                                         , 8 , "All"  , ["exactComparison", "no_restart"]) ,
   ("cloudBM80GLC"                      , "cloudBM80GLC.ups"                                        , 8 , "All"  , ["exactComparison", "no_restart"]) ,
   ("cloudBM80LS"                       , "cloudBM80LS.ups"                                         , 8 , "All"  , ["exactComparison", "no_restart"]) ,
   ("multibox_sweeps_coal"              , "Coal/multibox_sweeps_coal.ups"                           , 46, "All"  , ["exactComparison"]),
#__________________________________
# Kokkos tests
   ("kokkos-x-scalar"                   , "kokkos_solver_tests/kokkos-x-scalar.ups"                , 1.1 , "All" , ["exactComparison"]), 
   ("kokkos-y-scalar"                   , "kokkos_solver_tests/kokkos-y-scalar.ups"                , 1.1 , "All" , ["exactComparison"]), 
   ("kokkos-z-scalar"                   , "kokkos_solver_tests/kokkos-z-scalar.ups"                , 1.1 , "All" , ["exactComparison"]), 
#__________________________________
# THESE TESTS FAIL TO RUN TO COMPLETION ON A CUDA ENABLED BUILD   "corrupted double-linked list: 0x00000000024b8120 ***"
   ("coal_channel_hi_vel"               , "Coal/coal_channel_hi_vel.ups"                            , 1.1 , "All"  , ["exactComparison"     , "no_cuda"]) , 
   ("1GW_RT"                            , "Coal/1GW_RT.ups"                                         , 2   , "All"  , ["exactComparison"     , "no_cuda"]) , 
   ("1GW_em_tc"                         , "Coal/1GW_em_tc.ups"                                      , 2   , "All"  , ["exactComparison"     , "no_cuda"]) , 
#   ("coal_channel_Yamamoto"             , "Coal/coal_channel_Yamamoto.ups"                          , 1.1 , "Linux" , ["exactComparison"     , "no_restart"]) , 
   ("coal_RadPropsPlanck"               , "Coal/coal_RadPropsPlanck.ups"                            , 1.1 , "All"  , ["exactComparison"     , "no_cuda"]) , 
   ("OFC4"                              , "Coal/OFC4.ups"                                           , 3   , "All"  , ["exactComparison"      ,"do_performance_test",  "no_cuda"   ]) ,
   ("OFC4c"                             , "Coal/OFC4.ups"                                           , 3   , "All"  , ["exactComparison"      , "no_cuda"]) ,
   ("OFC4_smith"                        , "Coal/OFC4_smith.ups"                                     , 3   , "All"  , ["exactComparison"      , "no_cuda"]) ,
   ("OFC_smith2016"                     , "Coal/OFC_smith2016.ups"                                  , 3   , "All"  , ["exactComparison"      , "no_cuda"]) ,
   ("OFC4_hybrid"                       , "Coal/OFC4_hybrid.ups"                                    , 3   , "All"  , ["exactComparison"      , "no_cuda"]) ,  
#__________________________________


# multi-threaded NIGHTLY tests
   ("rmcrt_bm1_1L_thread"               , "RMCRT/rmcrt_bm1_1L.ups"                                  , 1.1 , "Linux"  , ["no_restart", "exactComparison", "sus_options=-nthreads 4"]),
   ("rmcrt_bm1_ML_thread"               , "RMCRT/rmcrt_bm1_ML.ups"                                  , 1.1 , "Linux"  , ["no_restart", "exactComparison", "sus_options=-nthreads 4"]),
   ("rmcrt_bm1_DO_thread"               , "RMCRT/rmcrt_bm1_DO.ups"                                  , 1.1 , "Linux"  , ["no_restart", "exactComparison", "sus_options=-nthreads 8"]),
   
   # The regrid test should be last.  It needs a checkpoint.  If you move it up the stack and run local_RT NIGHTLYTESTS then not all tests will run
   ("regridTestArches"                  , "regridTestArches"                                        , 8   , "Linux"  , ["startFromCheckpoint" , "no_restart"])
]

# Tests that are run during local regression testing
LOCALTESTS = [
   ("constantMMS"                       , "mms/constantMMS.ups"                                     , 1.1 , "All"  , ["exactComparison"]) , 
   ("almgrenMMS"                        , "mms/almgrenMMS.ups"                                      , 1.1 , "All"  , ["exactComparison"]) , 
   ("oned_pulse_conv"                   , "mms/oned_pulse_conv.ups"                                 , 1.1 , "All"  , ["exactComparison"])   , 
   ("isotropic-turbulence-decay"        , "periodicTurb/isotropic-turbulence-decay.ups"             , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("helium_1m"                         , "helium_1m.ups"                                           , 1.1 , "All"  , ["exactComparison"]) , 
   ("methane_fire"                      , "methane_fire.ups"                                        , 1.1 , "All"  , ["exactComparison"]) , 
   ("methane_fire_8patch"               , "methane_fire_8patch.ups"                                 , 8   , "All"  , ["exactComparison"]) , 
   ("methane_fire_8patch_petscrad"      , "methane_fire_8patch_petscrad.ups"                        , 8   , "All"  , ["exactComparison"]) , 
   ("rmcrt_bm1_1L"                      , "RMCRT/rmcrt_bm1_1L.ups"                                  , 1.1 , "All"  , ["exactComparison"]) , 
   ("rmcrt_bm1_1L_bounded"               , "RMCRT/rmcrt_bm1_1L_bounded.ups"                         , 8   , "ALL"  , ["exactComparison"]) ,
   ("rmcrt_bm1_DO"                      , "RMCRT/rmcrt_bm1_DO.ups"                                  , 8   , "ALL"  , ["exactComparison"]) , 
   ("rmcrt_bm1_ML"                      , "RMCRT/rmcrt_bm1_ML.ups"                                  , 1.1 , "ALL"  , ["exactComparison"]) , 
   ("methane_rmcrt"                     , "RMCRT/methane_rmcrt.ups"                                 , 8   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("dqmom_test_1"                      , "DQMOM_regression/dqmom_test_1.ups"                       , 1.1 , "All"  , ["exactComparison"]) , 
   ("dqmom_test_2"                      , "DQMOM_regression/dqmom_test_2.ups"                       , 1.1 , "All"  , ["exactComparison"]) , 
   ("dqmom_test_3"                      , "DQMOM_regression/dqmom_test_3.ups"                       , 1.1 , "All"  , ["exactComparison"]) , 
   ("dqmom_test_4"                      , "DQMOM_regression/dqmom_test_4.ups"                       , 1.1 , "All"  , ["exactComparison"]) , 
   ("dqmom_test_5"                      , "DQMOM_regression/dqmom_test_5.ups"                       , 1.1 , "All"  , ["exactComparison"]) , 
   ("birth_test"                        , "DQMOM_regression/birth_test.ups"                         , 1.1 , "All"  , ["exactComparison"])   , 
   ("upwind_birth_test"                 , "DQMOM_regression/upwind_birth_test.ups"                  , 1.1 , "All"  , ["exactComparison"])   , 
   ("methane_jet"                       , "ClassicMixingTables/ups/methane_jet.ups"                 , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("acidbase_jet"                      , "ClassicMixingTables/ups/acidbase_jet.ups"                , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("acidbase_jet_2D"                   , "ClassicMixingTables/ups/acidbase_jet_2D.ups"             , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("1DTableTest"                       , "ClassicMixingTables/ups/1DTableTest.ups"                 , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("xplus_scalar_test"                 , "ScalarTests/xplus_scalar_test.ups"                       , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("yplus_scalar_test"                 , "ScalarTests/yplus_scalar_test.ups"                       , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("zplus_scalar_test"                 , "ScalarTests/zplus_scalar_test.ups"                       , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("xminus_scalar_test"                , "ScalarTests/xminus_scalar_test.ups"                      , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("yminus_scalar_test"                , "ScalarTests/yminus_scalar_test.ups"                      , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("zminus_scalar_test"                , "ScalarTests/zminus_scalar_test.ups"                      , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("turbulent_inlet_test"              , "DigitalFilter/TurbulentInletChannel.ups"                 , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("source_channel_conden"             , "verify_masssource/source_channel_conden.ups"             , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("source_channel_varden"             , "verify_masssource/source_channel_varden.ups"             , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("heptane_pipe"                      , "heptane_pipe.ups"                                        , 1.1 , "All"  , ["exactComparison"]) , 
   ("methane_RCCE"                      , "methane_RCCE.ups"                                        , 1.1 , "All"  , ["exactComparison"]) , 
   ("channel_WD_CO"                     , "channel_WD_CO.ups"                                       , 1.1 , "All"  , ["exactComparison"]) , 
   ("coal_table_pipe"                   , "coal_table_pipe.ups"                                     , 1.1 , "All"  , ["exactComparison"]) , 
   ("pcoal_drag"                        , "Coal/pcoal_drag.ups"                                     , 1.1 , "All"  , ["exactComparison"]) , 
   ("scalar_var_1eqn"                   , "scalar_variance_1eqn.ups"                                , 4   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("scalar_var_2eqn"                   , "scalar_variance_2eqn.ups"                                , 4   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("coal_channel_FOWY"                 , "Coal/coal_channel_FOWY.ups"                              , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("DOM16"                             , "DOM16.ups"                                               , 3   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("DO_RadProps"                       , "DO_RadProps.ups"                                         , 3   , "All"  , ["exactComparison"]) , 
   ("CQMOM_1x1"                         , "CQMOM_regression/CQMOM_1x1.ups"                          , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("CQMOM_scalar_transport"            , "CQMOM_regression/CQMOM_Transport.ups"                    , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("CQMOM_scalar_transport2x2x2"       , "CQMOM_regression/CQMOM_Transport_2x2x2.ups"              , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("SimpleBoxPTC"                      , "CQMOM_regression/SimpleBoxPTC.ups"                       , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("PTC_2D"                            , "CQMOM_regression/PTC_2D.ups"                             , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("PTC_3D"                            , "CQMOM_regression/PTC_3D.ups"                             , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("CQMOM_4D"                          , "CQMOM_regression/CQMOM_4D.ups"                           , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) ,
   ("CQMOM_7D"                          , "CQMOM_regression/CQMOM_7D.ups"                           , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) ,
   ("singleJet_poly"                    , "CQMOM_regression/singleJet_poly.ups"                     , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) ,
   ("angledWall"                        , "CQMOM_regression/angledWall.ups"                         , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) ,
   ("angledWall3D"                      , "CQMOM_regression/angledWall3D.ups"                       , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) ,
   ("Constant_Deposition"               , "CQMOM_regression/Constant_Deposition.ups"                , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) ,
   ("CQMOM_coal_test"                   , "CQMOM_regression/CQMOM_coal_test.ups"                    , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) ,
   ("channel_LagPart_inlet"             , "LagrangianParticles/channel_flow_x_lagrangian_inlet.ups" , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) ,   
   ("partMassFlow"                      , "Coal/partMassFlow.ups"                                   , 8   , "All"  , ["exactComparison"   ]) , 
   ("task_math"                         , "task_math.ups"                                           , 1.1 , "All"  , ["exactComparison", "no_restart"]) , 
   ("intrusion_test"                    , "intrusion_test.ups"                                      , 1.1 , "All"  , ["exactComparison"]) , 
   ("multi-patch-intrusion-test"        , "multi-patch-intrusion-test.ups"                          , 8   , "All"  , ["exactComparison"]) , 
#__________________________________
# Kokkos tests
   ("kokkos-x-scalar"                   , "kokkos_solver_tests/kokkos-x-scalar.ups"                , 1.1 , "All" , ["exactComparison"]), 
   ("kokkos-y-scalar"                   , "kokkos_solver_tests/kokkos-y-scalar.ups"                , 1.1 , "All" , ["exactComparison"]), 
   ("kokkos-z-scalar"                   , "kokkos_solver_tests/kokkos-z-scalar.ups"                , 1.1 , "All" , ["exactComparison"]), 
#__________________________________
# THESE TESTS FAIL TO RUN TO COMPLETION ON A CUDA ENABLED BUILD, "corrupted double-linked list: 0x00000000024b8120 ***"
   ("coal_RadPropsPlanck"               , "Coal/coal_RadPropsPlanck.ups"                            , 1.1 , "All" , ["exactComparison"    ,  "no_cuda"]) ,
   ("OFC4"                              , "Coal/OFC4.ups"                                           , 3   , "All"  , ["exactComparison"   ,  "no_cuda"]) ,
   ("coal_channel_hi_vel"               , "Coal/coal_channel_hi_vel.ups"                            , 1.1 , "All"  , ["exactComparison"   ,  "no_cuda"]) , 
   ("1GW_RT"                            , "Coal/1GW_RT.ups"                                         , 2   , "All" , ["exactComparison"    ,  "no_cuda"]) , 
   ("1GW_em_tc"                         , "Coal/1GW_em_tc.ups"                                      , 2   , "All"  , ["exactComparison"     , "no_cuda"]) , 
#   ("coal_channel_Yamamoto"             , "Coal/coal_channel_Yamamoto.ups"                          , 1.1 , "All" , ["exactComparison"     , "no_restart"]) , 
   ("OFC4_smith"                        , "Coal/OFC4_smith.ups"                                     , 3   , "All"  , ["exactComparison"   ,  "no_cuda"]) , 
   ("OFC_smith2016"                     , "Coal/OFC_smith2016.ups"                                  , 3   , "All"  , ["exactComparison"      , "no_cuda"]) ,
   ("OFC4_hybrid"                       , "Coal/OFC4_hybrid.ups"                                    , 3   , "All"  , ["exactComparison"   ,  "no_cuda"]) , 
#__________________________________
   
# multi-threaded LOCAL tests
   ("rmcrt_bm1_1L_thread"                , "RMCRT/rmcrt_bm1_1L.ups"               , 1.1 , "Linux"    , ["no_restart", "exactComparison", "sus_options=-nthreads 4"]),
   ("rmcrt_bm1_ML_thread"                , "RMCRT/rmcrt_bm1_ML.ups"               , 1.1 , "Linux"    , ["no_restart", "exactComparison", "sus_options=-nthreads 4"]),
   ("rmcrt_bm1_DO_thread"                , "RMCRT/rmcrt_bm1_DO.ups"               , 1.1 , "Linux"    , ["no_restart", "exactComparison", "sus_options=-nthreads 8"]),   
]

DEBUG = [
   ("DO_RadProps"                       , "DO_RadProps.ups"                                         , 3 , "Linux" , ["exactComparison"])   , 
]

SCALARTESTS = [
   ("xplus_scalar_test"          , "ScalarTests/xplus_scalar_test.ups"           , 6   , "All"  , ["exactComparison", "no_restart"]) , 
   ("yplus_scalar_test"          , "ScalarTests/yplus_scalar_test.ups"           , 6   , "All"  , ["exactComparison", "no_restart"]) , 
   ("zplus_scalar_test"          , "ScalarTests/zplus_scalar_test.ups"           , 6   , "All"  , ["exactComparison", "no_restart"]) , 
   ("xminus_scalar_test"         , "ScalarTests/xminus_scalar_test.ups"          , 6   , "All"  , ["exactComparison", "no_restart"]) , 
   ("yminus_scalar_test"         , "ScalarTests/yminus_scalar_test.ups"          , 6   , "All"  , ["exactComparison", "no_restart"]) , 
   ("zminus_scalar_test"         , "ScalarTests/zminus_scalar_test.ups"          , 6   , "All"  , ["exactComparison", "no_restart"])
]


DQMOMTESTS = [
   ("dqmom_test_1"               , "DQMOM_regression/dqmom_test_1.ups"           , 1.1 , "All"   , ["exactComparison"]) , 
   ("dqmom_test_2"               , "DQMOM_regression/dqmom_test_2.ups"           , 1.1 , "All"   , ["exactComparison"]) , 
   ("dqmom_test_3"               , "DQMOM_regression/dqmom_test_3.ups"           , 1.1 , "All"   , ["exactComparison"]) , 
   ("dqmom_test_4"               , "DQMOM_regression/dqmom_test_4.ups"           , 1.1 , "All"   , ["exactComparison"]) , 
   ("dqmom_test_5"               , "DQMOM_regression/dqmom_test_5.ups"           , 1.1 , "All"   , ["exactComparison"]) 
]

COALTESTS = [
   ("coal_RadPropsPlanck"               , "Coal/coal_RadPropsPlanck.ups"                            , 1.1 , "All" , ["exactComparison"    ,  "no_cuda"]) ,
   ("OFC4"                              , "Coal/OFC4.ups"                                           , 3   , "All"  , ["exactComparison"   ,  "no_cuda"]) ,
   ("coal_channel_hi_vel"               , "Coal/coal_channel_hi_vel.ups"                            , 1.1 , "All"  , ["exactComparison"   ,  "no_cuda"]) , 
   ("1GW_RT"                            , "Coal/1GW_RT.ups"                                         , 2   , "All" , ["exactComparison"    ,  "no_cuda"]) , 
   ("1GW_em_tc"                         , "Coal/1GW_em_tc.ups"                                      , 2   , "All"  , ["exactComparison"     , "no_cuda"]) , 
#   ("coal_channel_Yamamoto"             , "Coal/coal_channel_Yamamoto.ups"                          , 1.1 , "All" , ["exactComparison"     , "no_restart"]) , 
   ("OFC4_smith"                        , "Coal/OFC4_smith.ups"                                     , 3   , "All"  , ["exactComparison"   ,  "no_cuda"]) , 
   ("OFC_smith2016"                     , "Coal/OFC_smith2016.ups"                                  , 3   , "All"  , ["exactComparison"      , "no_cuda"]) ,
   ("OFC4_hybrid"                       , "Coal/OFC4_hybrid.ups"                                    , 3   , "All"  , ["exactComparison"   ,  "no_cuda"]) ,
   ("OFC4_initial_hotwall"              , "Coal/OFC4_initial_hotwall.ups"                           , 3   , "All"  , ["exactComparison"   ,  "no_cuda"]) , 
]

RMCRTTESTS = [
   ("rmcrt_bm1_1L"                , "RMCRT/rmcrt_bm1_1L.ups"                      , 1.1 , "ALL"  , ["exactComparison"]) ,
   ("rmcrt_bm1_1L_bounded"         , "RMCRT/rmcrt_bm1_1L_bounded.ups"             , 8   , "ALL"  , ["exactComparison"]) ,
   ("rmcrt_bm1_ML"                , "RMCRT/rmcrt_bm1_ML.ups"                      , 1.1 , "ALL"  , ["exactComparison"]) , 
   ("rmcrt_bm1_DO"                , "RMCRT/rmcrt_bm1_DO.ups"                      , 1.1 , "ALL"  , ["exactComparison"]) ,
   ("methane_rmcrt"               , "RMCRT/methane_rmcrt.ups"                     , 8   , "ALL"  , ["exactComparison"     , "no_restart"]) ,

# multi-threaded RMCRT tests
   ("rmcrt_bm1_1L_thread"                , "RMCRT/rmcrt_bm1_1L.ups"               , 1.1 , "ALL"    , ["no_restart", "exactComparison", "sus_options=-nthreads 4"]),
   ("rmcrt_bm1_ML_thread"                , "RMCRT/rmcrt_bm1_ML.ups"               , 1.1 , "ALL"    , ["no_restart", "exactComparison", "sus_options=-nthreads 4"]),
   ("rmcrt_bm1_DO_thread"                , "RMCRT/rmcrt_bm1_DO.ups"               , 1.1 , "ALL"    , ["no_restart", "exactComparison", "sus_options=-nthreads 8"]),
]

KOKKOSTESTS = [
   ("kokkos-x-scalar"                   , "kokkos_solver_tests/kokkos-x-scalar.ups"                , 1.1 , "All" , ["exactComparison"]), 
   ("kokkos-y-scalar"                   , "kokkos_solver_tests/kokkos-y-scalar.ups"                , 1.1 , "All" , ["exactComparison"]), 
   ("kokkos-z-scalar"                   , "kokkos_solver_tests/kokkos-z-scalar.ups"                , 1.1 , "All" , ["exactComparison"]), 
]

CQMOMTESTS = [
  ("CQMOM_1x1"                   , "CQMOM_regression/CQMOM_1x1.ups"               , 1.1 , "All"  , ["exactComparison", "no_restart"]),
  ("CQMOM_scalar_transport"      , "CQMOM_regression/CQMOM_Transport.ups"         , 6   , "All"  , ["exactComparison", "no_restart"]),
  ("CQMOM_scalar_transport2x2x2" , "CQMOM_regression/CQMOM_Transport_2x2x2.ups"   , 6   , "All"  , ["exactComparison", "no_restart"]),
  ("SimpleBoxPTC"                , "CQMOM_regression/SimpleBoxPTC.ups"            , 1.1 , "All"  , ["exactComparison", "no_restart"]),
  ("PTC_2D"                      , "CQMOM_regression/PTC_2D.ups"                  , 1.1 , "All"  , ["exactComparison", "no_restart"]),
  ("PTC_3D"                      , "CQMOM_regression/PTC_3D.ups"                  , 1.1 , "All"  , ["exactComparison", "no_restart"]),
  ("CQMOM_4D"                    , "CQMOM_regression/CQMOM_4D.ups"                , 1.1 , "All"  , ["exactComparison", "no_restart"]),
  ("CQMOM_7D"                    , "CQMOM_regression/CQMOM_7D.ups"                , 1.1 , "All"  , ["exactComparison", "no_restart"]),
  ("singleJet_poly"              , "CQMOM_regression/singleJet_poly.ups"          , 1.1 , "All"  , ["exactComparison", "no_restart"]),
  ("angledWall"                  , "CQMOM_regression/angledWall.ups"              , 1.1 , "All"  , ["exactComparison", "no_restart"]),
  ("angledWall3D"                , "CQMOM_regression/angledWall3D.ups"            , 1.1 , "All"  , ["exactComparison", "no_restart"]),
  ("Constant_Deposition"         , "CQMOM_regression/Constant_Deposition.ups"     , 1.1 , "All"  , ["exactComparison", "no_restart"]),
  ("CQMOM_coal_test"             , "CQMOM_regression/CQMOM_coal_test.ups"         , 1.1 , "All"  , ["exactComparison", "no_restart"]),
]

# NO RMCRT due to the segfault on the MAC
NORMCRT = [
   ("constantMMS__NEW"                  , "mms/constantMMS__NEW.ups"                                , 1.1 , "All"  , ["exactComparison"]) , 
   ("almgrenMMS__NEW"                   , "mms/almgrenMMS__NEW.ups"                                 , 1.1 , "All"  , ["exactComparison"]) , 
   ("isotropic-turbulence-decay__NEW"   , "periodicTurb/isotropic-turbulence-decay__NEW.ups"        , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("helium_1m__NEW"                    , "helium_1m__NEW.ups"                                      , 1.1 , "All"  , ["exactComparison"]) , 
   ("methane_fire__NEW"                 , "methane_fire__NEW.ups"                                   , 1.1 , "All"  , ["exactComparison"]) , 
   ("methane_fire_8patch__NEW"          , "methane_fire_8patch__NEW.ups"                            , 8   , "All"  , ["exactComparison"]) , 
   ("methane_fire_8patch_petscrad__NEW" , "methane_fire_8patch_petscrad__NEW.ups"                   , 8   , "All"  , ["exactComparison"]) , 
   ("dqmom_test_1"                      , "DQMOM_regression/dqmom_test_1.ups"                       , 1.1 , "All"  , ["exactComparison"]) , 
   ("dqmom_test_2"                      , "DQMOM_regression/dqmom_test_2.ups"                       , 1.1 , "All"  , ["exactComparison"]) , 
   ("dqmom_test_3"                      , "DQMOM_regression/dqmom_test_3.ups"                       , 1.1 , "All"  , ["exactComparison"]) , 
   ("dqmom_test_4"                      , "DQMOM_regression/dqmom_test_4.ups"                       , 1.1 , "All"  , ["exactComparison"]) , 
   ("dqmom_test_5"                      , "DQMOM_regression/dqmom_test_5.ups"                       , 1.1 , "All"  , ["exactComparison"]) , 
   ("birth_test"                        , "DQMOM_regression/birth_test.ups"                         , 1.1 , "All"  , ["exactComparison"])   , 
   ("upwind_birth_test"                 , "DQMOM_regression/upwind_birth_test.ups"                  , 1.1 , "All"  , ["exactComparison"])   , 
   ("methane_jet"                       , "ClassicMixingTables/ups/methane_jet.ups"                 , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("acidbase_jet"                      , "ClassicMixingTables/ups/acidbase_jet.ups"                , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("acidbase_jet_2D"                   , "ClassicMixingTables/ups/acidbase_jet_2D.ups"             , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("1DTableTest"                       , "ClassicMixingTables/ups/1DTableTest.ups"                 , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("xplus_scalar_test"                 , "ScalarTests/xplus_scalar_test.ups"                       , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("yplus_scalar_test"                 , "ScalarTests/yplus_scalar_test.ups"                       , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("zplus_scalar_test"                 , "ScalarTests/zplus_scalar_test.ups"                       , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("xminus_scalar_test"                , "ScalarTests/xminus_scalar_test.ups"                      , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("yminus_scalar_test"                , "ScalarTests/yminus_scalar_test.ups"                      , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("zminus_scalar_test"                , "ScalarTests/zminus_scalar_test.ups"                      , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("turbulent_inlet_test"              , "DigitalFilter/TurbulentInletChannel.ups"                 , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("source_channel_conden"             , "verify_masssource/source_channel_conden.ups"             , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("source_channel_varden"             , "verify_masssource/source_channel_varden.ups"             , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("heptane_pipe"                      , "heptane_pipe.ups"                                        , 1.1 , "All"  , ["exactComparison"]) , 
   ("methane_RCCE"                      , "methane_RCCE.ups"                                        , 1.1 , "All " , ["exactComparison"]) , 
   ("channel_WD_CO"                     , "channel_WD_CO.ups"                                       , 1.1 , "All"  , ["exactComparison"]) , 
   ("coal_table_pipe"                   , "coal_table_pipe.ups"                                     , 1.1 , "All"  , ["exactComparison"]) , 
   ("scalar_var_1eqn"                   , "scalar_variance_1eqn.ups"                                , 4   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("scalar_var_2eqn"                   , "scalar_variance_2eqn.ups"                                , 4   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("coal_channel_FOWY"                 , "Coal/coal_channel_FOWY.ups"                              , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("coal_RadPropsPlanck"               , "Coal/coal_RadPropsPlanck.ups"                            , 1.1 , "All"  , ["exactComparison"]), 
   ("pcoal_drag"                        , "Coal/pcoal_drag.ups"                                     , 1.1 , "All"  , ["exactComparison"]) , 
   ("DOM16"                             , "DOM16.ups"                                               , 3   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("DO_RadProps"                       , "DO_RadProps.ups"                                         , 3   , "All"  , ["exactComparison"]) , 
   ("CQMOM_1x1"                         , "CQMOM_regression/CQMOM_1x1.ups"                          , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("CQMOM_scalar_transport"            , "CQMOM_regression/CQMOM_Transport.ups"                    , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("CQMOM_scalar_transport2x2x2"       , "CQMOM_regression/CQMOM_Transport_2x2x2.ups"              , 6   , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("SimpleBoxPTC"                      , "CQMOM_regression/SimpleBoxPTC.ups"                       , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("PTC_2D"                            , "CQMOM_regression/PTC_2D.ups"                             , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("PTC_3D"                            , "CQMOM_regression/PTC_3D.ups"                             , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("CQMOM_4D"                          , "CQMOM_regression/CQMOM_4D.ups"                           , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
   ("channel_LagPart_inlet"             , "LagrangianParticles/channel_flow_x_lagrangian_inlet.ups" , 1.1 , "All"  , ["exactComparison"   , "no_restart"]) , 
]

#__________________________________
# The following list is parsed by the local RT script
# and allows the user to select the tests to run
#LIST: LOCALTESTS KOKKOSTESTS RMCRTTESTS NEWTESTS SCALARTESTS DQMOMTESTS NIGHTLYTESTS CQMOMTESTS NORMCRT DEBUG COALTESTS
#__________________________________

  
# returns the list  
def getTestList(me) :
  if me == "LOCALTESTS":
    TESTS = LOCALTESTS
  elif me == "KOKKOSTESTS": 
    TESTS = KOKKOSTESTS
  elif me == "RMCRTTESTS":
    TESTS = RMCRTTESTS
  elif me == "SCALARTESTS":
    TESTS = SCALARTESTS
  elif me == "COALTESTS":
    TESTS = COALTESTS
  elif me == "DQMOMTESTS":
    TESTS = DQMOMTESTS
  elif me == "NIGHTLYTESTS":
    TESTS = NIGHTLYTESTS
  elif me == "CQMOMTESTS":
    TESTS = CQMOMTESTS
  elif me == "NORMCRT":
    TESTS = NORMCRT
  elif me == "DEBUG":
    TESTS = DEBUG
  else:
    print "\nERROR:ARCHES.py  getTestList:  The test list (%s) does not exist!\n\n" % me
    exit(1)
  return TESTS

#__________________________________

if __name__ == "__main__":

  TESTS = getTestList( environ['WHICH_TESTS'] )

  result = runSusTests(argv, TESTS, "ARCHES")
  
  # cleanup modified files
  command = "/bin/rm -rf %s/tmp > /dev/null 2>&1 " % (the_dir)
  system( command )
  
  exit( result )


