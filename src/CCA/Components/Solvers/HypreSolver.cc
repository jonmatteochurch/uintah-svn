/*
 * The MIT License
 *
 * Copyright (c) 1997-2019 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */



#include <CCA/Components/Solvers/HypreSolver.h>
#include <CCA/Components/Solvers/MatrixUtil.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/SoleVariable.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/Stencil4.h>
#include <Core/Grid/Variables/Stencil7.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/ConvergenceFailure.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/LoadBalancer.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Math/MiscMath.h>
#include <Core/Math/MinMax.h>
#include <Core/Util/DebugStream.h>
#include <Core/Util/StringUtil.h>
#include <Core/Util/Timers/Timers.hpp>
#include <iomanip>

// hypre includes
#include <_hypre_struct_mv.h>
#include <_hypre_utilities.h>
#include <HYPRE_struct_ls.h>
#include <krylov.h>

#ifndef HYPRE_TIMING
#ifndef hypre_ClearTiming
// This isn't in utilities.h for some reason...
#define hypre_ClearTiming()
#endif
#endif

using namespace std;
using namespace Uintah;

//__________________________________
//  To turn on normal output
//  setenv SCI_DEBUG "SOLVER_DOING_COUT:+"

static DebugStream cout_doing("SOLVER_DOING_COUT", false);
static DebugStream cout_assembling("SOLVER ASSEMBLING", false);

namespace Uintah {

  void swapbytes( Uintah::hypre_solver_structP& ) {
    SCI_THROW(InternalError("Swap bytes for hypre_solver_structP is not implemented", __FILE__, __LINE__));
  }

  //==============================================================================
  //
  // Class HypreStencil7
  //
  //==============================================================================
  template<class GridVarType>
  class HypreStencil7 : public RefCounted {
  public:
    HypreStencil7( const Level          * level_in
                 , const MaterialSet    * matlset_in
                 , const VarLabel       * A_in
                 ,       Task::WhichDW    which_A_dw_in
                 , const VarLabel       * x_in
                 ,       bool             modifies_X_in
                 , const VarLabel       * b_in
                 ,       Task::WhichDW    which_b_dw_in
                 , const VarLabel       * guess_in
                 ,       Task::WhichDW    which_guess_dw_in
                 , const HypreParams    * params_in
                 ,       bool             isFirstSolve_in
                 )
      : m_level(level_in)
      , m_matlset(matlset_in)
      , m_A_label(A_in)
      , m_which_A_dw(which_A_dw_in)
      , m_X_label(x_in)
      , m_modifies_X(modifies_X_in)
      , m_b_label(b_in)
      , m_which_b_dw(which_b_dw_in)
      , m_guess_label(guess_in)
      , m_which_guess_dw(which_guess_dw_in)
      , m_params(params_in)
      , m_isFirstSolve(isFirstSolve_in)
    {
      // Time Step
      m_timeStepLabel    = VarLabel::create(timeStep_name, timeStep_vartype::getTypeDescription() );
      m_firstPassThrough = true;
      m_movingAverage    = 0.0;
    }

    //---------------------------------------------------------------------------------------------

    virtual ~HypreStencil7() {
      VarLabel::destroy(m_timeStepLabel);
    }


    //---------------------------------------------------------------------------------------------
    void getPatchExtents( const Patch      * patch
                          ,     IntVector  & lo
                          ,     IntVector  & hi)
    {
      typedef typename GridVarType::double_type double_type;
      Patch::VariableBasis basis = Patch::translateTypeToBasis(double_type::getTypeDescription()->getType(), true);

      if( m_params->getSolveOnExtraCells()) {
        lo  = patch->getExtraLowIndex(  basis, IntVector(0,0,0) );
        hi  = patch->getExtraHighIndex( basis, IntVector(0,0,0) );
      } else {
        lo = patch->getLowIndex(  basis );
        hi = patch->getHighIndex( basis );
      }
    }

    //---------------------------------------------------------------------------------------------
    //   Create and populate a Hypre struct vector,
    HYPRE_StructVector
    createPopulateHypreVector(   const timeStep_vartype   timeStep
                               , const bool               recompute
                               , const bool               do_setup
                               , const ProcessorGroup   * pg
                               , const HYPRE_StructGrid & grid
                               , const PatchSubset      * patches
                               , const int                matl
                               , const VarLabel         * Q_label
                               , DataWarehouse          * Q_dw
                               , HYPRE_StructVector     * HQ)
    {
      //__________________________________
      // Create the vector
      if (timeStep == 1 || recompute) {
        cout_assembling << "hyper "<< Q_label->getName() <<" vector setup" << std::endl;
        HYPRE_StructVectorCreate( pg->getComm(), grid, HQ );
        HYPRE_StructVectorInitialize( *HQ );
      }
      else if (do_setup) {
        cout_assembling << "hyper "<< Q_label->getName() <<" vector destroy + setup" << std::endl;
        HYPRE_StructVectorDestroy( *HQ );
        HYPRE_StructVectorCreate( pg->getComm(), grid, HQ );
        HYPRE_StructVectorInitialize( *HQ );
      }

      cout_assembling << "hyper "<< Q_label->getName() <<" vector update coefficients" << std::endl;
      for(int p=0;p<patches->size();p++){
        const Patch* patch = patches->get(p);

        //__________________________________
        // Get Q
        if( Q_label ){
          ostringstream msg;
          msg<< "HypreSolver:createPopulateHypreVector ("<< Q_label->getName() <<")\n";
          printTask( patches, patch, cout_doing, msg.str() );

          typename GridVarType::const_double_type Q;
          Q_dw->get( Q, Q_label, matl, patch, Ghost::None, 0);

          // find box range
          IntVector lo;
          IntVector hi;
          getPatchExtents( patch, lo, hi );

          //__________________________________
          // Feed Q variable to Hypre
          for(int z=lo.z(); z<hi.z(); z++){
            for(int y=lo.y(); y<hi.y(); y++){

              IntVector l(lo.x(), y, z);
              IntVector h(hi.x()-1, y, z);

              const double* values = &Q[l];

              HYPRE_StructVectorSetBoxValues( *HQ,
                                             l.get_pointer(), h.get_pointer(),
                                             const_cast<double*>(values));
            }
          }
        }  // label exist?
      }  // patch loop

      if (timeStep == 1 || recompute || do_setup){
        cout_assembling << "hyper "<< Q_label->getName() <<" vector assemble" << std::endl;
        HYPRE_StructVectorAssemble( *HQ );
      }

      return *HQ;
    }

    //---------------------------------------------------------------------------------------------

    void solve( const ProcessorGroup * pg
              , const PatchSubset    * patches
              , const MaterialSubset * matls
              ,       DataWarehouse  * old_dw
              ,       DataWarehouse  * new_dw
              , const VarLabel       * hypre_solver_label
              ,       Handle<HypreStencil7<GridVarType>>
              )
    {
      //__________________________________
      //   timers
      m_tHypreAll = hypre_InitializeTiming("Total Hypre time");
      hypre_BeginTiming(m_tHypreAll);

      m_tMatVecSetup = hypre_InitializeTiming("Matrix + Vector setup");
      m_tSolveOnly   = hypre_InitializeTiming("Solve time");


       //________________________________________________________
      // get struct from data warehouse
      struct hypre_solver_struct* hypre_solver_s = 0;

      if ( new_dw->exists( hypre_solver_label ) ) {
        cout_assembling << "hyper solver exists in new DW" << std::endl;
        new_dw->get( m_hypre_solverP, hypre_solver_label );
      }
      else if ( old_dw->exists( hypre_solver_label ) ) {
        cout_assembling << "hyper solver moved from old to new DW" << std::endl;
        old_dw->get( m_hypre_solverP, hypre_solver_label );
        new_dw->put( m_hypre_solverP, hypre_solver_label );
      }
      else {
        cout_assembling << "hyper solver crated in new DW" << std::endl;

        hypre_solver_struct* hypre_struct = scinew hypre_solver_struct;

        hypre_struct->solver_p         = scinew HYPRE_StructSolver( nullptr );
        hypre_struct->precond_solver_p = scinew HYPRE_StructSolver( nullptr );
        hypre_struct->HA_p             = scinew HYPRE_StructMatrix( nullptr );
        hypre_struct->HX_p             = scinew HYPRE_StructVector( nullptr );
        hypre_struct->HB_p             = scinew HYPRE_StructVector( nullptr );
        hypre_struct->solver_type         = HypreSolver2::stringToSolverType( m_params->solvertype );
        hypre_struct->precond_solver_type = HypreSolver2::stringToSolverType( m_params->precondtype );

        hypre_struct->isRecomputeTimeStep = true;

        m_hypre_solverP.setData( hypre_struct );
        new_dw->put( m_hypre_solverP, hypre_solver_label );
      }
      hypre_solver_s = m_hypre_solverP.get().get_rep();

      bool recompute = hypre_solver_s->isRecomputeTimeStep;

      //__________________________________
      // timestep can come from the old_dw or parentOldDW
      timeStep_vartype timeStep(0);

      Task::WhichDW myOldDW = m_params->getWhichOldDW();
      DataWarehouse* pOldDW = new_dw->getOtherDataWarehouse(myOldDW);

      pOldDW->get(timeStep, m_timeStepLabel);

      //________________________________________________________
      // Solve frequency
      //
      const int solvFreq = m_params->solveFrequency;
      // note - the first timeStep in hypre is timeStep 1
      if (solvFreq == 0 || timeStep % solvFreq ) {
        new_dw->transferFrom(old_dw, m_X_label, patches, matls, true);
        return;
      }

      //________________________________________________________
      // Matrix setup frequency - this will destroy and recreate a new Hypre matrix at the specified setupFrequency
      //
      int suFreq = m_params->getSetupFrequency();
      bool do_setup = false;
      if (suFreq != 0){
        do_setup = (timeStep % suFreq == 0);
      }

      //________________________________________________________
      // update coefficient frequency - This will ONLY UPDATE the matrix coefficients without destroying/recreating the Hypre Matrix
      //
      const int updateCoefFreq = m_params->getUpdateCoefFrequency();
     bool updateCoefs = true;
      if (updateCoefFreq != 0){
        updateCoefs = (timeStep % updateCoefFreq == 0);
      }

      // if it the first pass through ignore the flags
      if(timeStep == 1 || recompute){
        updateCoefs = false;
        do_setup    = false;
      }

      //std::cout << "      HypreSolve  timestep: " << timeStep << " recompute: " << recompute << " m_firstPassThrough: " << m_firstPassThrough <<  " m_isFirstSolve: " << m_isFirstSolve <<" do_setup: " << do_setup << " updateCoefs: " << updateCoefs << std::endl;

      DataWarehouse* A_dw     = new_dw->getOtherDataWarehouse( m_which_A_dw );
      DataWarehouse* b_dw     = new_dw->getOtherDataWarehouse( m_which_b_dw );
      DataWarehouse* guess_dw = new_dw->getOtherDataWarehouse( m_which_guess_dw );

      // ASSERTEQ(sizeof(Stencil7), 7*sizeof(double));

      Timers::Simple timer;
      timer.start();

      for(int m = 0;m<matls->size();m++){
        int matl = matls->get(m);

        hypre_BeginTiming(m_tMatVecSetup);
        //__________________________________
        // Setup grid
        HYPRE_StructGrid grid;
        if (timeStep == 1 || do_setup || recompute) {
          cout_assembling << "hyper grid setup" << std::endl;
          HYPRE_StructGridCreate(pg->getComm(), 3, &grid);

          for(int p=0;p<patches->size();p++){
            const Patch* patch = patches->get(p);

            IntVector lo;
            IntVector hi;
            getPatchExtents( patch, lo, hi );
            hi -= IntVector(1,1,1);

            HYPRE_StructGridSetExtents(grid, lo.get_pointer(), hi.get_pointer());
          }

          // Periodic boundaries
          if (!patches->empty()) {
            const Level* level = getLevel(patches);
            IntVector periodic_vector = level->getPeriodicBoundaries();

            IntVector low, high;
            level->findCellIndexRange(low, high);
            IntVector range = high-low;

            int periodic[3];
            periodic[0] = periodic_vector.x() * range.x();
            periodic[1] = periodic_vector.y() * range.y();
            periodic[2] = periodic_vector.z() * range.z();
            HYPRE_StructGridSetPeriodic(grid, periodic);
          }

          // Assemble the grid
          HYPRE_StructGridAssemble(grid);
        }

        //__________________________________
        // Create the stencil
        HYPRE_StructStencil stencil;
        if ( timeStep == 1 || do_setup || recompute) {
          cout_assembling << "hyper stencil setup" << std::endl;
          if( m_params->getSymmetric()){

            HYPRE_StructStencilCreate(3, 4, &stencil);
            int offsets[4][3] = {{0,0,0},
              {-1,0,0},
              {0,-1,0},
              {0,0,-1}};
            for(int i=0;i<4;i++) {
              HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
            }

          } else {

            HYPRE_StructStencilCreate(3, 7, &stencil);
            int offsets[7][3] = {{0,0,0},
              {1,0,0}, {-1,0,0},
              {0,1,0}, {0,-1,0},
              {0,0,1}, {0,0,-1}};

            for(int i=0;i<7;i++){
              HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
            }
          }
        }

        //__________________________________
        // Create the matrix
        HYPRE_StructMatrix* HA = hypre_solver_s->HA_p;

        if (timeStep == 1 || recompute) {
          cout_assembling << "hyper matrix setup" << std::endl;
          HYPRE_StructMatrixCreate( pg->getComm(), grid, stencil, HA );
          HYPRE_StructMatrixSetSymmetric( *HA, m_params->getSymmetric() );
          int ghost[] = {1,1,1,1,1,1};
          HYPRE_StructMatrixSetNumGhost( *HA, ghost );
          HYPRE_StructMatrixInitialize( *HA );
        }
        else if (do_setup) {
          cout_assembling << "hyper matrix destroy + setup" << std::endl;
          HYPRE_StructMatrixDestroy( *HA );
          HYPRE_StructMatrixCreate(pg->getComm(), grid, stencil, HA );
          HYPRE_StructMatrixSetSymmetric( *HA, m_params->getSymmetric() );
          int ghost[] = {1,1,1,1,1,1};
          HYPRE_StructMatrixSetNumGhost( *HA, ghost );
          HYPRE_StructMatrixInitialize(*HA );
        }

        // setup the coefficient matrix ONLY on the first timeStep, if
        // we are doing a recompute, or if we set setupFrequency != 0,
        // or if UpdateCoefFrequency != 0
        if (timeStep == 1 || recompute || do_setup || updateCoefs) {
          cout_assembling << "hyper matrix update coefficients" << std::endl;
          for(int p=0;p<patches->size();p++) {
            const Patch* patch = patches->get(p);
            printTask( patches, patch, cout_doing, "HypreSolver:solve: Create Matrix" );

            //__________________________________
            // Get A matrix from the DW
            typename GridVarType::symmetric_matrix_type AStencil4;
            typename GridVarType::matrix_type A;

            if ( m_params->getUseStencil4() ){
              A_dw->get( AStencil4, m_A_label, matl, patch, Ghost::None, 0);
            } else {
              A_dw->get( A, m_A_label, matl, patch, Ghost::None, 0);
            }

            IntVector l;
            IntVector h;
            getPatchExtents( patch, l, h );

            //__________________________________
            // Feed it to Hypre
            if( m_params->getSymmetric()){

              double* values = scinew double[(h.x()-l.x())*4];
              int stencil_indices[] = {0,1,2,3};

              // use stencil4 as coefficient matrix. NOTE: This should be templated
              // on the stencil type. This workaround is to get things moving
              // until we convince component developers to move to stencil4. You must
              // set m_params->setUseStencil4(true) when you setup your linear solver
              // if you want to use stencil4. You must also provide a matrix of type
              // stencil4 otherwise this will crash.
              if ( m_params->getUseStencil4()) {

                for(int z=l.z();z<h.z();z++){
                  for(int y=l.y();y<h.y();y++){

                    const Stencil4* AA = &AStencil4[IntVector(l.x(), y, z)];
                    double* p = values;

                    for(int x=l.x();x<h.x();x++){
                      *p++ = AA->p;
                      *p++ = AA->w;
                      *p++ = AA->s;
                      *p++ = AA->b;
                      AA++;
                    }
                    IntVector ll(l.x(), y, z);
                    IntVector hh(h.x()-1, y, z);
                    HYPRE_StructMatrixSetBoxValues(*HA,
                                                   ll.get_pointer(), hh.get_pointer(),
                                                   4, stencil_indices, values);

                  } // y loop
                }  // z loop

              } else { // use stencil7

                for(int z=l.z();z<h.z();z++){
                  for(int y=l.y();y<h.y();y++){

                    const Stencil7* AA = &A[IntVector(l.x(), y, z)];
                    double* p = values;

                    for(int x=l.x();x<h.x();x++){
                      *p++ = AA->p;
                      *p++ = AA->w;
                      *p++ = AA->s;
                      *p++ = AA->b;
                      AA++;
                    }
                    IntVector ll(l.x(), y, z);
                    IntVector hh(h.x()-1, y, z);
                    HYPRE_StructMatrixSetBoxValues(*HA,
                                                   ll.get_pointer(), hh.get_pointer(),
                                                   4, stencil_indices, values);

                  } // y loop
                }  // z loop
              }
              delete[] values;
            } else {
              double* values = scinew double[(h.x()-l.x())*7];
              int stencil_indices[] = {0,1,2,3,4,5,6};

              for(int z=l.z();z<h.z();z++){
                for(int y=l.y();y<h.y();y++){

                  const Stencil7* AA = &A[IntVector(l.x(), y, z)];
                  double* p = values;

                  for(int x=l.x();x<h.x();x++){
                    *p++ = AA->p;
                    *p++ = AA->e;
                    *p++ = AA->w;
                    *p++ = AA->n;
                    *p++ = AA->s;
                    *p++ = AA->t;
                    *p++ = AA->b;
                    AA++;
                  }

                  IntVector ll(l.x(), y, z);
                  IntVector hh(h.x()-1, y, z);
                  HYPRE_StructMatrixSetBoxValues(*HA,
                                                 ll.get_pointer(), hh.get_pointer(),
                                                 7, stencil_indices,
                                                 values);
                }  // y loop
              } // z loop
              delete[] values;
            }
          }
          if (timeStep == 1 || recompute || do_setup){
            cout_assembling << "hyper matrix assemble" << std::endl;
            HYPRE_StructMatrixAssemble(*HA);
          }
        }

        //__________________________________
        // Create the RHS
        HYPRE_StructVector HB;
        HB = createPopulateHypreVector(  timeStep, recompute, do_setup, pg, grid, patches, matl, m_b_label, b_dw, hypre_solver_s->HB_p);

        //__________________________________
        // Create the solution vector
        HYPRE_StructVector HX;
        HX = createPopulateHypreVector(  timeStep, recompute, do_setup, pg, grid, patches, matl, m_guess_label, guess_dw, hypre_solver_s->HX_p);

        hypre_EndTiming( m_tMatVecSetup );

        //__________________________________
        Timers::Simple solve_timer;
        solve_timer.start();

        hypre_BeginTiming(m_tSolveOnly);

        int num_iterations;
        double final_res_norm;

        //______________________________________________________________________
        // Solve the system
        switch( hypre_solver_s->solver_type ){
        //______________________________________________________________________
        //
        case jacobi: {

          HYPRE_StructSolver * solver = hypre_solver_s->solver_p;

          if ( do_setup ){
            HYPRE_StructJacobiDestroy(*solver);
          }

          if ( timeStep == 1 || recompute || do_setup ) {
            HYPRE_StructJacobiCreate(pg->getComm(),solver);

            HYPRE_StructJacobiSetTol    (*solver, m_params->tolerance);
            HYPRE_StructJacobiSetMaxIter(*solver, m_params->maxiterations);

            HYPRE_StructJacobiSetup(*solver,*HA,HB,HX);
          }

          HYPRE_StructJacobiSolve(*solver,*HA,HB,HX);

          HYPRE_StructJacobiGetNumIterations(*solver, &num_iterations);
          HYPRE_StructJacobiGetFinalRelativeResidualNorm(*solver, &final_res_norm);
          break;
        }
        //______________________________________________________________________
        //
        case pfmg: {

          HYPRE_StructSolver * solver = hypre_solver_s->solver_p;

          if ( do_setup ){
            HYPRE_StructPFMGDestroy(*solver);
          }

          if ( timeStep == 1 || recompute || do_setup ) {
            HYPRE_StructPFMGCreate(pg->getComm(),solver);

            HYPRE_StructPFMGSetTol         (*solver, m_params->tolerance);
            HYPRE_StructPFMGSetMaxIter     (*solver, m_params->maxiterations);
            HYPRE_StructPFMGSetMaxLevels   (*solver, m_params->maxlevels);
            HYPRE_StructPFMGSetRelChange   (*solver, 0);
            HYPRE_StructPFMGSetRelaxType   (*solver, m_params->relax_type); // Jacobi = 0; weighted Jacobi = 1; symm red-black GS = 2; nonsymm red-black GS = 3
            HYPRE_StructPFMGSetJacobiWeight(*solver, m_params->jacobiweight);
            HYPRE_StructPFMGSetRAPType     (*solver, m_params->raptype); // Galerkin = 0; non-Galerkin = 1
            HYPRE_StructPFMGSetNumPreRelax (*solver, m_params->npre);
            HYPRE_StructPFMGSetNumPostRelax(*solver, m_params->npost);
            HYPRE_StructPFMGSetSkipRelax   (*solver, m_params->skip);
            HYPRE_StructPFMGSetLogging     (*solver, m_params->logging);

            HYPRE_StructPFMGSetup(*solver,*HA,HB,HX);
          }

          HYPRE_StructPFMGSolve(*solver,*HA,HB,HX);

          HYPRE_StructPFMGGetNumIterations(*solver, &num_iterations);
          HYPRE_StructPFMGGetFinalRelativeResidualNorm(*solver, &final_res_norm);
          break;
        }
        //__________________________________
        // use symmetric SMG
        case smg: {

          HYPRE_StructSolver * solver = hypre_solver_s->solver_p;

          if ( do_setup ){
            HYPRE_StructSMGDestroy(*solver);
          }

          if ( timeStep == 1 || recompute || do_setup ) {
            HYPRE_StructSMGCreate(pg->getComm(), solver);

            HYPRE_StructSMGSetMemoryUse   (*solver, 0);
            HYPRE_StructSMGSetTol         (*solver, m_params->tolerance);
            HYPRE_StructSMGSetMaxIter     (*solver, m_params->maxiterations);
            HYPRE_StructSMGSetRelChange   (*solver, 0);
            HYPRE_StructSMGSetNumPreRelax (*solver, m_params->npre);
            HYPRE_StructSMGSetNumPostRelax(*solver, m_params->npost);
            HYPRE_StructSMGSetLogging     (*solver, m_params->logging);

            HYPRE_StructSMGSetup(*solver,*HA,HB,HX);
          }

          HYPRE_StructSMGSolve(*solver,*HA,HB,HX);

          HYPRE_StructSMGGetNumIterations(*solver, &num_iterations);
          HYPRE_StructSMGGetFinalRelativeResidualNorm(*solver, &final_res_norm);
          break;
        }
        //______________________________________________________________________
        //
        case cycred: {

          HYPRE_StructSolver * solver = hypre_solver_s->solver_p;

          if ( do_setup ){
            HYPRE_StructCycRedDestroy(*solver);
          }
          if ( timeStep == 1 || recompute || do_setup ) {
            HYPRE_StructCycRedCreate(pg->getComm(),solver);

            HYPRE_StructCycRedSetup(*solver,*HA,HB,HX);
          }

          HYPRE_StructCycRedSolve(*solver,*HA,HB,HX);

          num_iterations=0;
          final_res_norm=0;
          break;
        }
        //______________________________________________________________________
        //
        case pcg: {

          HYPRE_StructSolver * solver         = hypre_solver_s->solver_p;
          HYPRE_StructSolver * precond_solver = hypre_solver_s->precond_solver_p;
          HYPRE_Solver       * krylov         = (HYPRE_Solver*)solver;

          if( do_setup ){
            destroyPrecond(hypre_solver_s, *precond_solver);
            HYPRE_StructPCGDestroy(*solver);
          }

          if ( timeStep == 1 || recompute || do_setup ) {
            HYPRE_StructPCGCreate(pg->getComm(),solver);

            HYPRE_PtrToStructSolverFcn precond;
            HYPRE_PtrToStructSolverFcn precond_setup;

            setupPrecond(pg, precond, precond_setup, hypre_solver_s, *precond_solver);
            if (precond) HYPRE_StructPCGSetPrecond(*solver, precond, precond_setup, *precond_solver);

            HYPRE_StructPCGSetTol           (*solver, m_params->tolerance);
            HYPRE_PCGSetAbsoluteTol         (*krylov, m_params->absolutetol);
            // HYPRE_PCGSetResidualTol      (*krylov, obsolete);
            HYPRE_PCGSetAbsoluteTolFactor   (*krylov, m_params->abstolf);
            HYPRE_PCGSetConvergenceFactorTol(*krylov, m_params->ctolf);
            // HYPRE_PCGSetStopCrit         (*krylov, obsolete);
            HYPRE_StructPCGSetMaxIter       (*solver, m_params->maxiterations);
            HYPRE_StructPCGSetTwoNorm       (*solver, 1);
            HYPRE_StructPCGSetRelChange     (*solver, 0);
            HYPRE_PCGSetRecomputeResidual   (*krylov, m_params->recomputeresidualperiod > 0 ? 1 : 0);
            HYPRE_PCGSetRecomputeResidualP  (*krylov, m_params->recomputeresidualperiod);
            HYPRE_StructPCGSetLogging       (*solver, m_params->logging);

            HYPRE_StructPCGSetup(*solver,*HA,HB,HX);
          }

          HYPRE_StructPCGSolve(*solver,*HA,HB,HX);

          HYPRE_StructPCGGetNumIterations(*solver, &num_iterations);
          HYPRE_StructPCGGetFinalRelativeResidualNorm(*solver,&final_res_norm);
          break;
        }
        //______________________________________________________________________
        //
        case gmres: {

          HYPRE_StructSolver * solver         = hypre_solver_s->solver_p;
          HYPRE_StructSolver * precond_solver = hypre_solver_s->precond_solver_p;
          HYPRE_Solver       * krylov         = (HYPRE_Solver*)solver;

          if ( do_setup ){
            destroyPrecond(hypre_solver_s, *precond_solver);
            HYPRE_StructGMRESDestroy(*solver);
          }

          if ( timeStep == 1 || recompute || do_setup ) {
            HYPRE_StructGMRESCreate(pg->getComm(),solver);

            HYPRE_PtrToStructSolverFcn precond;
            HYPRE_PtrToStructSolverFcn precond_setup;

            setupPrecond(pg, precond, precond_setup, hypre_solver_s, *precond_solver);
            if (precond) HYPRE_StructGMRESSetPrecond(*solver, precond, precond_setup, *precond_solver);

            HYPRE_StructGMRESSetTol            (*solver, m_params->tolerance);
            HYPRE_StructGMRESSetAbsoluteTol    (*solver, m_params->absolutetol);
            HYPRE_GMRESSetConvergenceFactorTol (*krylov, m_params->ctolf);
            // HYPRE_GMRESSetStopCrit          (*krylov, obsolete);
            HYPRE_GMRESSetMinIter              (*krylov, m_params->miniterations);
            HYPRE_StructGMRESSetMaxIter        (*solver, m_params->maxiterations);
            HYPRE_StructGMRESSetKDim           (*solver, m_params->kdim);
            HYPRE_GMRESSetRelChange            (*krylov, 0);
            HYPRE_GMRESSetSkipRealResidualCheck(*krylov, m_params->skiprealrcheck);
            HYPRE_StructGMRESSetLogging        (*solver, m_params->logging);

            HYPRE_StructGMRESSetup(*solver,*HA,HB,HX);
          }

          HYPRE_StructGMRESSolve(*solver,*HA,HB,HX);

          HYPRE_StructGMRESGetNumIterations(*solver, &num_iterations);
          HYPRE_StructGMRESGetFinalRelativeResidualNorm(*solver, &final_res_norm);
          break;
        }
        //______________________________________________________________________
        //
        case flexgmres: {

          HYPRE_StructSolver * solver         = hypre_solver_s->solver_p;
          HYPRE_StructSolver * precond_solver = hypre_solver_s->precond_solver_p;
          HYPRE_Solver       * krylov         = (HYPRE_Solver*)solver;

          if ( do_setup ){
            destroyPrecond(hypre_solver_s, *precond_solver);
            HYPRE_StructFlexGMRESDestroy(*solver);
          }

          if ( timeStep == 1 || recompute || do_setup ) {
            HYPRE_StructFlexGMRESCreate(pg->getComm(),solver);

            HYPRE_PtrToStructSolverFcn precond;
            HYPRE_PtrToStructSolverFcn precond_setup;

            setupPrecond(pg, precond, precond_setup, hypre_solver_s, *precond_solver);
            if (precond) HYPRE_StructFlexGMRESSetPrecond(*solver, precond, precond_setup, *precond_solver);

            HYPRE_StructFlexGMRESSetTol           (*solver, m_params->tolerance);
            HYPRE_FlexGMRESSetAbsoluteTol         (*krylov, m_params->absolutetol);
            HYPRE_FlexGMRESSetConvergenceFactorTol(*krylov, m_params->ctolf);
            HYPRE_FlexGMRESSetMinIter             (*krylov, m_params->miniterations);
            HYPRE_StructFlexGMRESSetMaxIter       (*solver, m_params->maxiterations);
            HYPRE_StructFlexGMRESSetKDim          (*solver, m_params->kdim);
            HYPRE_StructFlexGMRESSetLogging       (*solver, m_params->logging);

            HYPRE_StructFlexGMRESSetup(*solver,*HA,HB,HX);
          }

          HYPRE_StructFlexGMRESSolve(*solver,*HA,HB,HX);

          HYPRE_StructFlexGMRESGetNumIterations(*solver, &num_iterations);
          HYPRE_StructFlexGMRESGetFinalRelativeResidualNorm(*solver, &final_res_norm);
          break;
        }
        //______________________________________________________________________
        //
        case lgmres: {

          HYPRE_StructSolver * solver         = hypre_solver_s->solver_p;
          HYPRE_StructSolver * precond_solver = hypre_solver_s->precond_solver_p;
          HYPRE_Solver       * krylov         = (HYPRE_Solver*)solver;

          if ( do_setup ){
            destroyPrecond(hypre_solver_s, *precond_solver);
            HYPRE_StructLGMRESDestroy(*solver);
          }
          if ( timeStep == 1 || recompute || do_setup ) {
            HYPRE_StructLGMRESCreate(pg->getComm(),solver);

            HYPRE_PtrToStructSolverFcn precond;
            HYPRE_PtrToStructSolverFcn precond_setup;

            setupPrecond(pg, precond, precond_setup, hypre_solver_s, *precond_solver);
            if (precond) HYPRE_StructLGMRESSetPrecond(*solver, precond, precond_setup, *precond_solver);

            HYPRE_StructLGMRESSetTol           (*solver, m_params->tolerance);
            HYPRE_LGMRESSetAbsoluteTol         (*krylov, m_params->absolutetol);
            HYPRE_LGMRESSetConvergenceFactorTol(*krylov, m_params->ctolf);
            HYPRE_LGMRESSetMinIter             (*krylov, m_params->miniterations);
            HYPRE_StructLGMRESSetMaxIter       (*solver, m_params->maxiterations);
            HYPRE_StructLGMRESSetKDim          (*solver, m_params->kdim);
            HYPRE_StructLGMRESSetAugDim        (*solver, m_params->augdim);
            HYPRE_StructLGMRESSetLogging       (*solver, m_params->logging);

            HYPRE_StructLGMRESSetup(*solver,*HA,HB,HX);
          }

          HYPRE_StructLGMRESSolve(*solver,*HA,HB,HX);

          HYPRE_StructLGMRESGetNumIterations(*solver, &num_iterations);
          HYPRE_StructLGMRESGetFinalRelativeResidualNorm(*solver, &final_res_norm);
          break;
        }
        //______________________________________________________________________
        //
        case bicgstab: {

          HYPRE_StructSolver * solver         = hypre_solver_s->solver_p;
          HYPRE_StructSolver * precond_solver = hypre_solver_s->precond_solver_p;
          HYPRE_Solver       * krylov         = (HYPRE_Solver*)solver;

          if ( do_setup ){
            destroyPrecond(hypre_solver_s, *precond_solver);
            HYPRE_StructBiCGSTABDestroy(*solver);
          }
          if ( timeStep == 1 || recompute || do_setup ) {
            HYPRE_StructBiCGSTABCreate(pg->getComm(),solver);

            HYPRE_PtrToStructSolverFcn precond;
            HYPRE_PtrToStructSolverFcn precond_setup;

            setupPrecond(pg, precond, precond_setup, hypre_solver_s, *precond_solver);
            if (precond) HYPRE_StructBiCGSTABSetPrecond    (*solver, precond, precond_setup, *precond_solver);

            HYPRE_StructBiCGSTABSetTol        (*solver, m_params->tolerance);
            HYPRE_StructBiCGSTABSetAbsoluteTol(*solver, m_params->absolutetol);
            HYPRE_BiCGSTABSetMinIter          (*krylov, m_params->miniterations);
            HYPRE_StructBiCGSTABSetMaxIter    (*solver, m_params->maxiterations);
            // HYPRE_BiCGSTABSetStopCrit      (*krylov, obsolete);
            HYPRE_StructBiCGSTABSetLogging    (*solver, m_params->logging);

            HYPRE_StructBiCGSTABSetup         (*solver,*HA,HB,HX);
          }

          HYPRE_StructBiCGSTABSolve(*solver,*HA,HB,HX);

          HYPRE_StructBiCGSTABGetNumIterations(*solver, &num_iterations);
          HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm(*solver, &final_res_norm);
          break;
        }
        //______________________________________________________________________
        //
        case hybrid: {

          HYPRE_StructSolver * solver         = hypre_solver_s->solver_p;
          HYPRE_StructSolver * precond_solver = hypre_solver_s->precond_solver_p;

          if ( do_setup ){
            destroyPrecond(hypre_solver_s, *precond_solver);
            HYPRE_StructHybridDestroy(*solver);
          }

          if ( timeStep == 1 || recompute || do_setup ) {
            HYPRE_StructHybridCreate(pg->getComm(), solver);

            HYPRE_PtrToStructSolverFcn precond;
            HYPRE_PtrToStructSolverFcn precond_setup;

            setupPrecond(pg, precond, precond_setup, hypre_solver_s, *precond_solver);
            if (precond) HYPRE_StructHybridSetPrecond(*solver, precond, precond_setup, *precond_solver);

            HYPRE_StructHybridSetTol               (*solver, m_params->tolerance);
            HYPRE_StructHybridSetConvergenceTol    (*solver, 0.90);
            HYPRE_StructHybridSetDSCGMaxIter       (*solver, 100);
            HYPRE_StructHybridSetPCGMaxIter        (*solver, m_params->maxiterations);
            HYPRE_StructHybridSetTwoNorm           (*solver, 1);
            // HYPRE_StructHybridSetStopCrit       (*solver, obsolete);
            HYPRE_StructHybridSetRelChange         (*solver, 0);
            HYPRE_StructHybridSetSolverType        (*solver, m_params->krylovsolver); // 0 = PCG; 1 = GMRES; 2: BiCGSTAB
            HYPRE_StructHybridSetRecomputeResidual (*solver, m_params->recomputeresidualperiod > 0 ? 1 : 0);
            HYPRE_StructHybridSetRecomputeResidualP(*solver, m_params->recomputeresidualperiod);
            HYPRE_StructHybridSetKDim              (*solver, m_params->kdim);
            HYPRE_StructHybridSetLogging           (*solver, m_params->logging);

            HYPRE_StructHybridSetup                (*solver,*HA,HB,HX);
          }

          HYPRE_StructHybridSolve(*solver,*HA,HB,HX);

          HYPRE_StructHybridGetNumIterations(*solver,&num_iterations);
          HYPRE_StructHybridGetFinalRelativeResidualNorm(*solver, &final_res_norm);
          break;
        }
        //______________________________________________________________________
        //
        case sparsemsg: {

          HYPRE_StructSolver * solver = hypre_solver_s->solver_p;
          if ( do_setup ){
            HYPRE_StructSparseMSGDestroy(*solver);
          }
          if ( timeStep == 1 || recompute || do_setup ) {
            HYPRE_StructSparseMSGCreate(pg->getComm(), solver);

            HYPRE_StructSparseMSGSetTol         (*solver, m_params->tolerance);
            HYPRE_StructSparseMSGSetMaxIter     (*solver, m_params->maxiterations);
            HYPRE_StructSparseMSGSetJump        (*solver, m_params->jump);
            HYPRE_StructSparseMSGSetRelChange   (*solver, 0);
            HYPRE_StructSparseMSGSetRelaxType   (*solver, m_params->relax_type);  // Jacobi = 0; weighted Jacobi = 1; symm red-black GS = 2; nonsymm red-black GS = 3
            HYPRE_StructSparseMSGSetJacobiWeight(*solver, m_params->jacobiweight);  // Jacobi = 0; weighted Jacobi = 1; symm red-black GS = 2; nonsymm red-black GS = 3
            HYPRE_StructSparseMSGSetNumPreRelax (*solver, m_params->npre);
            HYPRE_StructSparseMSGSetNumPostRelax(*solver, m_params->npost);
            HYPRE_StructSparseMSGSetNumFineRelax(*solver, m_params->nfine);
            HYPRE_StructSparseMSGSetLogging     (*solver, m_params->logging);

            HYPRE_StructSparseMSGSetup(*solver,*HA,HB,HX);
          }

          HYPRE_StructSparseMSGSolve(*solver,*HA,HB,HX);

          HYPRE_StructSparseMSGGetNumIterations(*solver, &num_iterations);
          HYPRE_StructSparseMSGGetFinalRelativeResidualNorm(*solver, &final_res_norm);
          break;
        }
        default:
          throw InternalError("Unknown solver type: "+ m_params->solvertype, __FILE__, __LINE__);
        }

        //______________________________________________________________________
        //
#ifdef PRINTSYSTEM
        //__________________________________
        //   Debugging
        vector<string> fname;
        m_params->getOutputFileName(fname);
        HYPRE_StructMatrixPrint( fname[0].c_str(), *HA, 0 );
        HYPRE_StructVectorPrint( fname[1].c_str(), *HB, 0 );
        HYPRE_StructVectorPrint( fname[2].c_str(), HX, 0 );
#endif

        printTask( patches, patches->empty() ? nullptr : patches->get(0), cout_doing, "HypreSolver:solve: testConvergence" );

        solve_timer.stop();
        hypre_EndTiming ( m_tSolveOnly );

        //__________________________________
        // Push the solution into Uintah data structure
        for(int p=0;p<patches->size();p++){
          const Patch* patch = patches->get(p);
          printTask( patches, patch, cout_doing, "HypreSolver:solve: copy solution" );

          IntVector l;
          IntVector h;
          getPatchExtents( patch, l, h );

          CellIterator iter(l, h);

          typename GridVarType::double_type Xnew;
          if( m_modifies_X ){
            new_dw->getModifiable(Xnew, m_X_label, matl, patch);
          }else{
            new_dw->allocateAndPut(Xnew, m_X_label, matl, patch);
          }

          // Get the solution back from hypre
          for(int z=l.z();z<h.z();z++){
            for(int y=l.y();y<h.y();y++){

              double* values = &Xnew[IntVector(l.x(), y, z)];
              IntVector ll(l.x(), y, z);
              IntVector hh(h.x()-1, y, z);

              HYPRE_StructVectorGetBoxValues(HX,
                  ll.get_pointer(), hh.get_pointer(),
                  values);
            }
          }
        }
        //__________________________________
        // clean up
         m_firstPassThrough  = false;
         hypre_solver_s->isRecomputeTimeStep  = false;

        if ( timeStep == 1 || do_setup || recompute ) {
          HYPRE_StructStencilDestroy(stencil);
          HYPRE_StructGridDestroy(grid);
        }

        hypre_EndTiming (m_tHypreAll);

        hypre_PrintTiming   ("Hypre Timings:", pg->getComm());
        hypre_FinalizeTiming( m_tMatVecSetup );
        hypre_FinalizeTiming( m_tSolveOnly );
        hypre_FinalizeTiming( m_tHypreAll );
        hypre_ClearTiming();

        timer.stop();

        if(pg->myRank() == 0) {

          cout << "Solve of " << m_X_label->getName()
               << " on level " << m_level->getIndex()
               << " completed in " << timer().seconds()
               << " s (solve only: " << solve_timer().seconds() << " s, ";

          if (timeStep > 2) {
            // alpha = 2/(N+1)
            // averaging window is 10 timeSteps.
            double alpha   = 2.0/(std::min( int(timeStep) - 2, 10) + 1);
            m_movingAverage = alpha*solve_timer().seconds() + (1-alpha) * m_movingAverage;

            cout << "mean: " <<  m_movingAverage << " s, ";
          }

          cout << num_iterations << " iterations, residual = "
               << final_res_norm << ")." << std::endl;
        }

        timer.reset( true );

        //__________________________________
        // Test for convergence failure

        if( final_res_norm > m_params->tolerance || std::isfinite(final_res_norm) == 0 ){
          if( m_params->getRecomputeTimeStepOnFailure() ){
            proc0cout << "  WARNING:  HypreSolver not converged in " << num_iterations
                      << " iterations, final residual= " << final_res_norm
                      << ", requesting the time step be recomputed.\n";

            new_dw->put( bool_or_vartype(true), VarLabel::find(abortTimeStep_name));
            new_dw->put( bool_or_vartype(true), VarLabel::find(recomputeTimeStep_name));
          } else {
            throw ConvergenceFailure("HypreSolver variable: "+ m_X_label->getName()+", solver: "+ m_params->solvertype+", preconditioner: "+ m_params->precondtype,
                                     num_iterations, final_res_norm,
                                     m_params->tolerance,__FILE__,__LINE__);
          }
        }
      }
    }

    //---------------------------------------------------------------------------------------------
    void
    setupPrecond( const ProcessorGroup              * pg
                 ,       HYPRE_PtrToStructSolverFcn & precond
                 ,       HYPRE_PtrToStructSolverFcn & pcsetup
                 ,struct hypre_solver_struct        * hypre_solver_s
                 ,       HYPRE_StructSolver         & precond_solver
                 )
    {
      switch( hypre_solver_s->precond_solver_type ){
      //__________________________________
      //  use diagonal scaling as preconditioner
      case none:{
        precond_solver = NULL;

        precond = NULL;
        pcsetup = NULL;
        break;
      }      //__________________________________
      //  use diagonal scaling as preconditioner
      case diagonal:{
        precond_solver = NULL;

        precond = HYPRE_StructDiagScale;
        pcsetup = HYPRE_StructDiagScaleSetup;
        break;
      }
      //__________________________________
      //  use two-step Jacobi as preconditioner
      case jacobi: {
        HYPRE_StructJacobiCreate(pg->getComm(),&precond_solver);

        HYPRE_StructJacobiSetTol      (precond_solver, m_params->precond_tolerance);
        HYPRE_StructJacobiSetMaxIter  (precond_solver, m_params->precond_maxiters);
        HYPRE_StructJacobiSetZeroGuess(precond_solver);

        precond = HYPRE_StructJacobiSolve;
        pcsetup = HYPRE_StructJacobiSetup;
        break;

      }
      //__________________________________
      // use symmetric PFMG as preconditioner
      case pfmg: {
        HYPRE_StructPFMGCreate(pg->getComm(),&precond_solver);

        HYPRE_StructPFMGSetTol         (precond_solver, m_params->precond_tolerance);
        HYPRE_StructPFMGSetMaxIter     (precond_solver, m_params->precond_maxiters);
        HYPRE_StructPFMGSetMaxLevels   (precond_solver, m_params->maxlevels);
        HYPRE_StructPFMGSetRelaxType   (precond_solver, m_params->relax_type); // Jacobi = 0; weighted Jacobi = 1; symm red-black GS = 2; nonsymm red-black GS = 3
        HYPRE_StructPFMGSetJacobiWeight(precond_solver, m_params->jacobiweight);
        HYPRE_StructPFMGSetRAPType     (precond_solver, m_params->raptype); // Galerkin = 0; non-Galerkin = 1
        HYPRE_StructPFMGSetNumPreRelax (precond_solver, m_params->npre);
        HYPRE_StructPFMGSetNumPostRelax(precond_solver, m_params->npost);
        HYPRE_StructPFMGSetSkipRelax   (precond_solver, m_params->skip);
        HYPRE_StructPFMGSetLogging     (precond_solver, 0);
        HYPRE_StructPFMGSetZeroGuess   (precond_solver);

        precond = HYPRE_StructPFMGSolve;
        pcsetup = HYPRE_StructPFMGSetup;
        break;
      }
      //__________________________________
      // use symmetric SMG as preconditioner
      case smg: {
        HYPRE_StructSMGCreate(pg->getComm(),&precond_solver);

        HYPRE_StructSMGSetMemoryUse   (precond_solver, 0);
        HYPRE_StructSMGSetTol         (precond_solver, m_params->precond_tolerance);
        HYPRE_StructSMGSetMaxIter     (precond_solver, m_params->precond_maxiters);
        HYPRE_StructSMGSetNumPreRelax (precond_solver, m_params->npre);
        HYPRE_StructSMGSetNumPostRelax(precond_solver, m_params->npost);
        HYPRE_StructSMGSetLogging     (precond_solver, 0);
        HYPRE_StructSMGSetZeroGuess   (precond_solver);

        precond = HYPRE_StructSMGSolve;
        pcsetup = HYPRE_StructSMGSetup;
        break;

      }
      //__________________________________
      //
      case pcg: {
        HYPRE_StructPCGCreate(pg->getComm(),&precond_solver);

        HYPRE_Solver precond_krylov = (HYPRE_Solver)&precond_solver;

        HYPRE_StructPCGSetTol           (precond_solver, m_params->precond_tolerance);
        HYPRE_PCGSetAbsoluteTol         (precond_krylov, m_params->precond_absolutetol);
        // HYPRE_PCGSetResidualTol      (precond_krylov, obsolete);
        HYPRE_PCGSetAbsoluteTolFactor   (precond_krylov, m_params->precond_abstolf);
        HYPRE_PCGSetConvergenceFactorTol(precond_krylov, m_params->precond_ctolf);
        // HYPRE_PCGSetStopCrit         (precond_krylov, obsolete);
        HYPRE_StructPCGSetMaxIter       (precond_solver, m_params->precond_maxiters);
        HYPRE_StructPCGSetTwoNorm       (precond_solver, 1);
        HYPRE_PCGSetRecomputeResidual   (precond_krylov, m_params->recomputeresidualperiod > 0 ? 1 : 0);
        HYPRE_PCGSetRecomputeResidualP  (precond_krylov, m_params->recomputeresidualperiod);
        HYPRE_StructPCGSetLogging       (precond_solver, 0);

        precond = HYPRE_StructPCGSolve;
        pcsetup = HYPRE_StructPCGSetup;
        break;
      }
      //__________________________________
      //
      case gmres: {
        HYPRE_StructGMRESCreate(pg->getComm(),&precond_solver);

        HYPRE_Solver precond_krylov = (HYPRE_Solver)&precond_solver;

        HYPRE_StructGMRESSetTol            (precond_solver, m_params->precond_tolerance);
        HYPRE_StructGMRESSetAbsoluteTol    (precond_solver, m_params->precond_absolutetol);
        HYPRE_GMRESSetConvergenceFactorTol (precond_krylov, m_params->precond_ctolf);
        // HYPRE_GMRESSetStopCrit          (precond_krylov, obsolete);
        HYPRE_GMRESSetMinIter              (precond_krylov, m_params->precond_miniters);
        HYPRE_StructGMRESSetMaxIter        (precond_solver, m_params->precond_maxiters);
        HYPRE_StructGMRESSetKDim           (precond_solver, m_params->kdim);
        HYPRE_GMRESSetSkipRealResidualCheck(precond_krylov, m_params->skiprealrcheck);
        HYPRE_StructGMRESSetLogging        (precond_solver, 0);

        precond = HYPRE_StructGMRESSolve;
        pcsetup = HYPRE_StructGMRESSetup;
        break;
      }
      //__________________________________
      //
      case flexgmres: {
        HYPRE_StructGMRESCreate(pg->getComm(),&precond_solver);

        HYPRE_Solver precond_krylov = (HYPRE_Solver)&precond_solver;

        HYPRE_StructGMRESSetTol               (precond_solver, m_params->precond_tolerance);
        HYPRE_FlexGMRESSetAbsoluteTol         (precond_krylov, m_params->precond_absolutetol);
        HYPRE_FlexGMRESSetConvergenceFactorTol(precond_krylov, m_params->precond_ctolf);
        HYPRE_FlexGMRESSetMinIter             (precond_krylov, m_params->precond_miniters);
        HYPRE_StructGMRESSetMaxIter           (precond_solver, m_params->precond_maxiters);
        HYPRE_StructGMRESSetKDim              (precond_solver, m_params->kdim);
        HYPRE_StructGMRESSetLogging           (precond_solver, 0);

        precond = HYPRE_StructFlexGMRESSolve;
        pcsetup = HYPRE_StructFlexGMRESSetup;
        break;
      }
      //__________________________________
      //
      case lgmres: {
        HYPRE_StructLGMRESCreate(pg->getComm(),&precond_solver);

        HYPRE_Solver precond_krylov = (HYPRE_Solver)&precond_solver;

        HYPRE_StructLGMRESSetTol           (precond_solver, m_params->precond_tolerance);
        HYPRE_LGMRESSetAbsoluteTol         (precond_krylov, m_params->precond_absolutetol);
        HYPRE_LGMRESSetConvergenceFactorTol(precond_krylov, m_params->precond_ctolf);
        HYPRE_LGMRESSetMinIter             (precond_krylov, m_params->precond_miniters);
        HYPRE_StructLGMRESSetMaxIter       (precond_solver, m_params->precond_maxiters);
        HYPRE_StructLGMRESSetKDim          (precond_solver, m_params->kdim);
        HYPRE_StructLGMRESSetAugDim        (precond_solver, m_params->augdim);
        HYPRE_StructLGMRESSetLogging       (precond_solver, 0);

        precond = HYPRE_StructLGMRESSolve;
        pcsetup = HYPRE_StructLGMRESSetup;
        break;
      }
      //__________________________________
      //
      case bicgstab: {
        HYPRE_StructBiCGSTABCreate(pg->getComm(),&precond_solver);

        HYPRE_Solver precond_krylov = (HYPRE_Solver)&precond_solver;

        HYPRE_StructBiCGSTABSetTol        (precond_solver, m_params->precond_tolerance);
        HYPRE_StructBiCGSTABSetAbsoluteTol(precond_solver, m_params->precond_absolutetol);
        HYPRE_BiCGSTABSetMinIter          (precond_krylov, m_params->precond_miniters);
        HYPRE_StructBiCGSTABSetMaxIter    (precond_solver, m_params->precond_maxiters);
        // HYPRE_BiCGSTABSetStopCrit      (precond_solver, obsolete);
        HYPRE_StructBiCGSTABSetLogging    (precond_solver, 0);

        precond = HYPRE_StructBiCGSTABSolve;
        pcsetup = HYPRE_StructBiCGSTABSetup;
        break;
      }
      //__________________________________
      //  use symmetric SparseMSG as preconditioner
      case sparsemsg:{
        HYPRE_StructSparseMSGCreate(pg->getComm(),&precond_solver);

        HYPRE_StructSparseMSGSetTol         (precond_solver, m_params->precond_tolerance);
        HYPRE_StructSparseMSGSetMaxIter     (precond_solver, m_params->precond_maxiters);
        HYPRE_StructSparseMSGSetJump        (precond_solver, m_params->jump);
        HYPRE_StructSparseMSGSetRelaxType   (precond_solver, m_params->relax_type);
        HYPRE_StructSparseMSGSetNumPreRelax (precond_solver, m_params->npre);
        HYPRE_StructSparseMSGSetNumPostRelax(precond_solver, m_params->npost);
        HYPRE_StructSparseMSGSetNumFineRelax(precond_solver, m_params->nfine);
        HYPRE_StructSparseMSGSetLogging     (precond_solver, 0);
        HYPRE_StructSparseMSGSetZeroGuess   (precond_solver);

        precond = HYPRE_StructSparseMSGSolve;
        pcsetup = HYPRE_StructSparseMSGSetup;
        break;

      }
      default:
        // This should have been caught in readParameters...
        throw InternalError("Unknown preconditionertype: "+ m_params->precondtype, __FILE__, __LINE__);
      }
    }

    //---------------------------------------------------------------------------------------------
    void
    destroyPrecond(struct hypre_solver_struct * hypre_solver_s
                   ,       HYPRE_StructSolver  & precond_solver )
    {

      switch( hypre_solver_s->precond_solver_type ){

      case none:
      case diagonal: {
        // do nothing
        break;
      }
      case jacobi: {
        HYPRE_StructJacobiDestroy( precond_solver );
        break;
      }
      case pfmg: {
        HYPRE_StructPFMGDestroy( precond_solver );
        break;
      }
      case smg: {
        HYPRE_StructSMGDestroy( precond_solver );
        break;
      }
      case pcg: {
        HYPRE_StructPCGDestroy( precond_solver );
        break;
      }
      case gmres: {
        HYPRE_StructGMRESDestroy( precond_solver );
        break;
      }
      case flexgmres: {
        HYPRE_StructFlexGMRESDestroy( precond_solver );
        break;
      }
      case lgmres: {
        HYPRE_StructLGMRESDestroy( precond_solver );
        break;
      }
      case bicgstab: {
        HYPRE_StructBiCGSTABDestroy( precond_solver );
        break;
      }
      case sparsemsg: {
        HYPRE_StructSparseMSGDestroy( precond_solver );
        break;
      }
      default:
        // This should have been caught in readParameters...
        throw InternalError("Unknown preconditionertype: "+ m_params->precondtype, __FILE__, __LINE__);
      }
    }

    //---------------------------------------------------------------------------------------------

  private:

    const Level*       m_level;
    const MaterialSet* m_matlset;
    const VarLabel*    m_A_label;
    Task::WhichDW      m_which_A_dw;
    const VarLabel*    m_X_label;
    bool               m_modifies_X;
    const VarLabel*    m_b_label;
    Task::WhichDW      m_which_b_dw;
    const VarLabel*    m_guess_label;
    Task::WhichDW      m_which_guess_dw;
    const HypreParams* m_params;
    bool               m_isFirstSolve;

    const VarLabel*    m_timeStepLabel;
    SoleVariable<hypre_solver_structP> m_hypre_solverP;
    bool   m_firstPassThrough;
    double m_movingAverage;

    // hypre timers - note that these variables do NOT store timings - rather, each corresponds to
    // a different timer index that is managed by Hypre. To enable the use and reporting of these
    // hypre timings, #define HYPRE_TIMING in HypreSolver.h
    int m_tHypreAll;    // Tracks overall time spent in Hypre = matrix/vector setup & assembly + solve time.
    int m_tSolveOnly;   // Tracks time taken by hypre to solve the system of equations
    int m_tMatVecSetup; // Tracks the time taken by uintah/hypre to allocate and set matrix and vector box vaules

  }; // class HypreStencil7

  //==============================================================================
  //
  // HypreSolver2 Implementation
  //
  //==============================================================================

  HypreSolver2::HypreSolver2(const ProcessorGroup* myworld)
  : SolverCommon(myworld)
  {
    // Time Step
    m_timeStepLabel = VarLabel::create(timeStep_name, timeStep_vartype::getTypeDescription() );

    if ( m_hypre_solver_label.find ( 0 ) == m_hypre_solver_label.end() )
        m_hypre_solver_label[0] = VarLabel::create ( "hypre_solver_label_" + std::to_string ( 0 ),
                                                            SoleVariable<hypre_solver_structP>::getTypeDescription() );

    m_params = scinew HypreParams();

  }

  //---------------------------------------------------------------------------------------------

  HypreSolver2::~HypreSolver2()
  {
    VarLabel::destroy(m_timeStepLabel);
    for (auto &&l : m_hypre_solver_label) VarLabel::destroy(l.second);
    delete m_params;
  }

  //---------------------------------------------------------------------------------------------

  void HypreSolver2::readParameters(       ProblemSpecP & params_ps
                                   , const string       & varname)
  {
    bool found=false;
    if(params_ps){
      for( ProblemSpecP param_ps = params_ps->findBlock("Parameters"); param_ps != nullptr; param_ps = param_ps->findNextBlock("Parameters")) {

        string variable;
        if( param_ps->getAttribute("variable", variable) && variable != varname ) {
          continue;
        }

        int sFreq;
        int coefFreq;
        string str_solver;
        string str_precond;

        param_ps->getWithDefault ("solver",        str_solver,      "smg");
        param_ps->getWithDefault ("preconditioner",str_precond,"diagonal");

        param_ps->getWithDefault ("setupFrequency",     sFreq,                   1);
        param_ps->getWithDefault ("updateCoefFrequency",coefFreq,                1);
        param_ps->getWithDefault ("solveFrequency",     m_params->solveFrequency,1);

        param_ps->getWithDefault ("tolerance",              m_params->tolerance,          1.e-10);
        param_ps->getWithDefault ("absolutetol",            m_params->absolutetol,           0.0);
        param_ps->getWithDefault ("abstolf",                m_params->abstolf,               0.0);
        param_ps->getWithDefault ("ctolf",                  m_params->ctolf,                 0.0);
        param_ps->getWithDefault ("miniterations",          m_params->miniterations,           0);
        param_ps->getWithDefault ("maxiterations",          m_params->maxiterations,          75);
        param_ps->getWithDefault ("precond_tolerance",      m_params->precond_tolerance,       0);
        param_ps->getWithDefault ("precond_absolutetol",    m_params->precond_absolutetol,   0.0);
        param_ps->getWithDefault ("precond_abstolf",        m_params->precond_abstolf,       0.0);
        param_ps->getWithDefault ("precond_ctolf",          m_params->precond_ctolf,         0.0);
        param_ps->getWithDefault ("precond_miniters",       m_params->precond_miniters,        0);
        param_ps->getWithDefault ("precond_maxiters",       m_params->precond_maxiters,        1);
        param_ps->getWithDefault ("logging",                m_params->logging,                 0);
        param_ps->getWithDefault ("maxlevels",              m_params->maxlevels,               0);
        param_ps->getWithDefault ("raptype",                m_params->raptype,                 0);
        param_ps->getWithDefault ("skip",                   m_params->skip,                    0);
        param_ps->getWithDefault ("relax_type",             m_params->relax_type,              1);
        param_ps->getWithDefault ("npre",                   m_params->npre,                    1);
        param_ps->getWithDefault ("npost",                  m_params->npost,                   1);
        param_ps->getWithDefault ("nfine",                  m_params->nfine,                   1);
        param_ps->getWithDefault ("recomputeresidualperiod",m_params->recomputeresidualperiod, 0);
        param_ps->getWithDefault ("kdim",                   m_params->kdim,                    5);
        param_ps->getWithDefault ("skiprealrcheck",         m_params->skiprealrcheck,          0);
        param_ps->getWithDefault ("augdim",                 m_params->augdim,                  2);
        param_ps->getWithDefault ("krylovsolver",           m_params->krylovsolver,            1);
        param_ps->getWithDefault ("jump",                   m_params->jump,                    0);

        // change to lowercase
        m_params->solvertype  = string_tolower( str_solver );
        m_params->precondtype = string_tolower( str_precond );

        m_params->setSetupFrequency(sFreq);
        m_params->setUpdateCoefFrequency(coefFreq);
        // Options from the HYPRE_ref_manual 2.8
        // npre:   Number of relaxation sweeps before coarse grid correction
        // npost:  Number of relaxation sweeps after coarse grid correction
        // skip:   Skip relaxation on certain grids for isotropic
        //         problems. This can greatly improve effciency by eliminating
        //         unnecessary relaxations when the underlying problem is isotropic.
        // jump:   not in manual
        //
        // relax_type
        // 0 : Jacobi
        // 1 : Weighted Jacobi (default)
        // 2 : Red/Black Gauss-Seidel (symmetric: RB pre-relaxation, BR post-relaxation)
        // 3 : Red/Black Gauss-Seidel (nonsymmetric: RB pre- and post-relaxation)

        found=true;
      }
    }
    if(!found){
      m_params->solvertype    = "smg";
      m_params->precondtype   = "diagonal";
      m_params->tolerance     = 1.e-10;
      m_params->maxiterations = 75;
      m_params->npre    = 1;
      m_params->npost   = 1;
      m_params->skip    = 0;
      m_params->jump    = 0;
      m_params->logging = 0;
      m_params->setSetupFrequency(1);
      m_params->setUpdateCoefFrequency(1);
      m_params->solveFrequency = 1;
      m_params->relax_type = 1;
    }
  }

  //---------------------------------------------------------------------------------------------

  void HypreSolver2::scheduleInitialize( const LevelP      & level
                                       ,       SchedulerP  & sched
                                       , const MaterialSet * matls
                                       )
  {
    const bool isRestart = sched->isRestartInitTimestep();
    Task* task = scinew Task("HypreSolver2::initialize_hypre", this, &HypreSolver2::initialize, isRestart );

    task->setType(Task::OncePerProc);  // must run this task on every proc.  It's possible to have
                                       // no patches on this proc when scheduling

    task->computes(hypre_solver_label(level));

    LoadBalancer * lb = sched->getLoadBalancer();

    sched->addTask(task, lb->getPerProcessorPatchSet(level), matls);
  }

 //---------------------------------------------------------------------------------------------

  void HypreSolver2::scheduleRestartInitialize( const LevelP      & level
                                              ,       SchedulerP  & sched
                                              , const MaterialSet * matls
                                              )
  {
    const bool isRestart = sched->isRestartInitTimestep();
    Task* task = scinew Task("HypreSolver2::restartInitialize_hypre", this, &HypreSolver2::initialize, isRestart);

    task->setType(Task::OncePerProc);  // must run this task on every proc.  It's possible to have
                                       // no patches  on this proc when scheduling restarts with regridding

    task->computes(hypre_solver_label(level));

    LoadBalancer * lb = sched->getLoadBalancer();

    sched->addTask(task, lb->getPerProcessorPatchSet(level), matls);
  }

  //---------------------------------------------------------------------------------------------

  void HypreSolver2::allocateHypreMatrices( DataWarehouse * new_dw,
                                            const Level   * level,
                                            const bool isRecomputeTimeStep_in )
  {
    SoleVariable<hypre_solver_structP> hypre_solverP;
    hypre_solver_struct* hypre_struct = scinew hypre_solver_struct;

    hypre_struct->solver_p         = scinew HYPRE_StructSolver( nullptr );
    hypre_struct->precond_solver_p = scinew HYPRE_StructSolver( nullptr );
    hypre_struct->HA_p             = scinew HYPRE_StructMatrix( nullptr );
    hypre_struct->HX_p             = scinew HYPRE_StructVector( nullptr );
    hypre_struct->HB_p             = scinew HYPRE_StructVector( nullptr );
    hypre_struct->solver_type         = stringToSolverType( m_params->solvertype );
    hypre_struct->precond_solver_type = stringToSolverType( m_params->precondtype );

    hypre_struct->isRecomputeTimeStep = isRecomputeTimeStep_in;

    hypre_solverP.setData( hypre_struct );
    new_dw->put( hypre_solverP, hypre_solver_label(level) );
  }

  //---------------------------------------------------------------------------------------------

  void
  HypreSolver2::initialize( const ProcessorGroup *
                          , const PatchSubset    * patches
                          , const MaterialSubset *
                          ,       DataWarehouse  *
                          ,       DataWarehouse  * new_dw
                          , const bool  isRestart
                          )
  {
    if ( patches->empty() ) return;
    allocateHypreMatrices( new_dw, getLevel(patches), isRestart );
  }

  //---------------------------------------------------------------------------------------------

  void
  HypreSolver2::scheduleSolve( const LevelP           & level
                             ,       SchedulerP       & sched
                             , const MaterialSet      * matls
                             , const VarLabel         * A_label
                             ,       Task::WhichDW      which_A_dw
                             , const VarLabel         * x_label
                             ,       bool               modifies_X
                             , const VarLabel         * b_label
                             ,       Task::WhichDW      which_b_dw
                             , const VarLabel         * guess_label
                             ,       Task::WhichDW      which_guess_dw
                             ,       bool               isFirstSolve /* = true */
                             )
  {
    printSchedule(level, cout_doing, "HypreSolver:scheduleSolve");

    Task* task;
    // The extra handle arg ensures that the stencil7 object will get freed
    // when the task gets freed.  The downside is that the refcount gets
    // tweaked everytime solve is called.

    TypeDescription::Type domtype = A_label->typeDescription()->getType();
    ASSERTEQ(domtype, x_label->typeDescription()->getType());
    ASSERTEQ(domtype, b_label->typeDescription()->getType());

    //__________________________________
    // bulletproofing
    IntVector periodic = level->getPeriodicBoundaries();
    if(periodic != IntVector(0,0,0)){

      IntVector l,h;
      level->findCellIndexRange( l, h );
      IntVector range = (h - l ) * periodic;

      if( fmodf(range.x(),2) != 0  || fmodf(range.y(),2) != 0 || fmodf(range.z(),2) != 0 ) {
        ostringstream warn;
        warn << "\nINPUT FILE WARNING: hypre solver: \n"
             << "With periodic boundary conditions the resolution of your grid "<<range<<", in each periodic direction, must be as close to a power of 2 as possible (i.e. M x 2^n).\n";

        if (m_params->solvertype == "smg") {
          throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
        }
        else {
          proc0cout << warn.str();
        }
      }
    }

    switch(domtype){
    case TypeDescription::SFCXVariable:
      {
        HypreStencil7<SFCXTypes>* that = scinew HypreStencil7<SFCXTypes>(level.get_rep(), matls, A_label, which_A_dw, x_label, modifies_X, b_label, which_b_dw, guess_label, which_guess_dw, m_params, isFirstSolve);
        Handle<HypreStencil7<SFCXTypes> > handle = that;
        task = scinew Task("Hypre:Matrix solve (SFCX)", that, &HypreStencil7<SFCXTypes>::solve, hypre_solver_label(level), handle);
      }
      break;
    case TypeDescription::SFCYVariable:
      {
        HypreStencil7<SFCYTypes>* that = scinew HypreStencil7<SFCYTypes>(level.get_rep(), matls, A_label, which_A_dw, x_label, modifies_X, b_label, which_b_dw, guess_label, which_guess_dw, m_params, isFirstSolve);
        Handle<HypreStencil7<SFCYTypes> > handle = that;
        task = scinew Task("Hypre:Matrix solve (SFCY)", that, &HypreStencil7<SFCYTypes>::solve, hypre_solver_label(level), handle);
      }
      break;
    case TypeDescription::SFCZVariable:
      {
        HypreStencil7<SFCZTypes>* that = scinew HypreStencil7<SFCZTypes>(level.get_rep(), matls, A_label, which_A_dw, x_label, modifies_X, b_label, which_b_dw, guess_label, which_guess_dw, m_params, isFirstSolve);
        Handle<HypreStencil7<SFCZTypes> > handle = that;
        task = scinew Task("Hypre:Matrix solve (SFCZ)", that, &HypreStencil7<SFCZTypes>::solve, hypre_solver_label(level), handle);
      }
      break;
    case TypeDescription::CCVariable:
      {
        HypreStencil7<CCTypes>* that = scinew HypreStencil7<CCTypes>(level.get_rep(), matls, A_label, which_A_dw, x_label, modifies_X, b_label, which_b_dw, guess_label, which_guess_dw, m_params, isFirstSolve);
        Handle<HypreStencil7<CCTypes> > handle = that;
        task = scinew Task("Hypre:Matrix solve (CC)", that, &HypreStencil7<CCTypes>::solve, hypre_solver_label(level), handle);
      }
      break;
    case TypeDescription::NCVariable:
      {
        HypreStencil7<NCTypes>* that = scinew HypreStencil7<NCTypes>(level.get_rep(), matls, A_label, which_A_dw, x_label, modifies_X, b_label, which_b_dw, guess_label, which_guess_dw, m_params, isFirstSolve);
        Handle<HypreStencil7<NCTypes> > handle = that;
        task = scinew Task("Hypre:Matrix solve (NC)", that, &HypreStencil7<NCTypes>::solve, hypre_solver_label(level), handle);
      }
      break;
    default:
      throw InternalError("Unknown variable type in scheduleSolve", __FILE__, __LINE__);
    }

    //__________________________________
    //  Computes and requires

    // Matrix A
    task->requires(which_A_dw, A_label, Ghost::None, 0);

    // Solution X
    if(modifies_X){
      task->modifies( x_label );
    } else {
      task->computes( x_label );
    }

    // Initial Guess
    if(guess_label){
      task->requires(which_guess_dw, guess_label, Ghost::None, 0);
    }

    // RHS  B
    task->requires(which_b_dw, b_label, Ghost::None, 0);

    // timestep
    // it could come from old_dw or parentOldDw
    Task::WhichDW old_dw = m_params->getWhichOldDW();
    task->requires( old_dw, m_timeStepLabel );

    // solve struct
    if (isFirstSolve) {
      task->requires( Task::OldDW, hypre_solver_label(level) );
      task->computes( hypre_solver_label(level) );
    }  else {
      task->requires( Task::NewDW, hypre_solver_label(level) );
    }

    sched->overrideVariableBehavior(hypre_solver_label(level)->getName(),false,false,false,true,true);

    task->setType(Task::Hypre);

    if( m_params->getRecomputeTimeStepOnFailure() ){
      task->computes( VarLabel::find(abortTimeStep_name) );
      task->computes( VarLabel::find(recomputeTimeStep_name) );
    }

    LoadBalancer * lb = sched->getLoadBalancer();

    sched->addTask(task, lb->getPerProcessorPatchSet(level), matls);
  }

  //---------------------------------------------------------------------------------------------

  string HypreSolver2::getName(){
    return "hypre";
  }

  //---------------------------------------------------------------------------------------------
  //  Return the solver or preconditioner type
  SolverType HypreSolver2::stringToSolverType( std::string str )
  {
    if ( str == "none" ){
      return none;
    }
    if ( str == "diagonal" ){
      return diagonal;
    }
    else if ( str == "jacobi" ){
      return jacobi;
    }
    else if ( str == "pfmg" ){
      return pfmg;
    }
    else if ( str == "smg" ){
      return smg;
    }
    else if ( str == "cycred" ){
      return cycred;
    }
    else if ( str == "cg" || str == "pcg" ){
      return pcg;
    }
    else if ( str == "gmres" ){
      return gmres;
    }
    else if ( str == "flexgmres" ){
      return flexgmres;
    }
    else if ( str == "lgmres" ){
      return lgmres;
    }
    else if ( str == "bicgstab" ){
      return bicgstab;
    }
    else if ( str == "hybrid" ){
      return hybrid;
    }
    else if ( str == "sparsemsg" ){
      return sparsemsg;
    }
    else {
      throw InternalError("ERROR:  Unknown solver type: "+ str, __FILE__, __LINE__);
    }
  }
  //---------------------------------------------------------------------------------------------
} // end namespace Uintah

