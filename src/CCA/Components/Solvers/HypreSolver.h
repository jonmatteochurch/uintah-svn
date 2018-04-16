/*
 * The MIT License
 *
 * Copyright (c) 1997-2018 The University of Utah
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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSolver_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSolver_h

//#define HYPRE_TIMING

#include <CCA/Components/Solvers/SolverCommon.h>

#include <Core/Grid/SimulationState.h>
#include <Core/Util/Handle.h>
#include <Core/Util/RefCounted.h>

#include <HYPRE_struct_ls.h>
#include <HYPRE_krylov.h>

#include <iostream>

/**
 *  @class  HypreSolver2
 *  @author Steve Parker
 *  @author Todd Haraman
 *  @author John Schmidt
 *  @author Tony Saad
 *  @author James Sutherland
 *  @author Oren Livne
 *  @brief  Uintah hypre solver interface.
 *  Allows the solution of a linear system of the form \[ \mathbf{A} \mathbf{x} = \mathbf{b}\] where \[\mathbf{A}\] is
 *  stencil7 matrix.
 *
 */

namespace Uintah {

  //______________________________________________________________________
  //
  class HypreSolver2Params : public SolverParameters {
  public:
    HypreSolver2Params(){}
    
    ~HypreSolver2Params() {}
    
    // Parameters common for all Hypre Solvers
    std::string solvertype;         // String corresponding to solver type
    std::string precondtype;        // String corresponding to preconditioner type
    double      tolerance;          // Residual tolerance for solver
    int         maxiterations;      // Maximum # iterations allowed
    int         logging;            // Log Hypre solver (using Hypre options)
    bool        restart;            // Allow solver to restart if not converged
    int         solveFrequency;     // Frequency for solving the linear system. timestep % solveFrequency
    int         relax_type;         // relaxation type
    
    // SMG parameters
    int    npre;               // # pre relaxations for Hypre SMG solver
    int    npost;              // # post relaxations for Hypre SMG solver
    
    // PFMG parameters
    int    skip;               // Hypre PFMG parameter
    
    // SparseMSG parameters
    int    jump;               // Hypre Sparse MSG parameter
  };

  //______________________________________________________________________
  //
  enum SolverType {
    smg,
    pfmg,
    sparsemsg,
    pcg,
    hybrid,
    gmres,
    jacobi,
    diagonal
  };

  //______________________________________________________________________
  //
  struct hypre_solver_struct : public RefCounted {
    bool                 created_solver;
    bool                 created_precond_solver;
    SolverType           solver_type;
    SolverType           precond_solver_type;
    HYPRE_StructSolver * solver;
    HYPRE_StructSolver * precond_solver;
    HYPRE_StructMatrix * HA;
    HYPRE_StructVector * HB;
    HYPRE_StructVector * HX;
    
    //__________________________________
    //
    hypre_solver_struct() {
      created_solver         = false;
      created_precond_solver = false;
      solver_type            = smg;
      precond_solver_type    = diagonal;
      solver                 = 0;
      precond_solver         = 0;
      HA = 0;
      HB = 0;
      HX = 0;
    };
    
    //__________________________________
    //
    virtual ~hypre_solver_struct() {
      if (created_solver) {
        HYPRE_StructMatrixDestroy( *HA );
        HYPRE_StructVectorDestroy( *HB );
        HYPRE_StructVectorDestroy( *HX );
      }
      if (created_solver)
        switch (solver_type) {
        case smg:
          HYPRE_StructSMGDestroy(*solver);
          break;
        case pfmg:
          HYPRE_StructPFMGDestroy(*solver);
          break;
        case sparsemsg:
          HYPRE_StructSparseMSGDestroy(*solver);
          break;
        case pcg:
          HYPRE_StructPCGDestroy(*solver);
          break;
        case gmres:
          HYPRE_StructGMRESDestroy(*solver);
          break;
        case jacobi:
          HYPRE_StructJacobiDestroy(*solver);
          break;
        default:
          throw InternalError( "HypreSolver given a bad solver type!", 
                               __FILE__, __LINE__ );
        }

      if (created_precond_solver)
        switch (precond_solver_type) {
        case smg:
          HYPRE_StructSMGDestroy(*precond_solver);
          break;
        case pfmg:
          HYPRE_StructPFMGDestroy(*precond_solver);
          break;
        case sparsemsg:
          HYPRE_StructSparseMSGDestroy(*precond_solver);
          break;
        case pcg:
          HYPRE_StructPCGDestroy(*precond_solver);
          break;
        case gmres:
          HYPRE_StructGMRESDestroy(*precond_solver);
          break;
        case jacobi:
          HYPRE_StructJacobiDestroy(*precond_solver);
          break;
        default:
          throw InternalError("HypreSolver given a bad solver type!", 
                              __FILE__, __LINE__);
      }

      if (HA) {
        delete HA;  
        HA = 0;
      }
      if (HB){
        delete HB;  
        HB = 0;
      }
      if (HX) {
        delete HX;  
        HX = 0;
      }
      if (solver) {
        delete solver;
        solver = 0;
      }
      if (precond_solver) {
        delete precond_solver;
        precond_solver = 0;
      }
    };
  };

  typedef Handle<hypre_solver_struct> hypre_solver_structP;
  
  //______________________________________________________________________
  //
  class HypreSolver2 : public SolverCommon {
  public:
    HypreSolver2(const ProcessorGroup* myworld);
    virtual ~HypreSolver2();

    virtual SolverParameters* readParameters(       ProblemSpecP & params,
                                              const std::string  & name  );

    /**
     *  @brief Schedules the solution of the linear system \[ \mathbf{A} \mathbf{x} = \mathbf{b}\].
     *
     *  @param level A reference to the level on which the system is solved.
     *  @param sched A reference to the Uintah scheduler.
     *  @param matls A pointer to the MaterialSet.
     *  @param A Varlabel of the coefficient matrix \[\mathbf{A}\]
     *  @param which_A_dw The datawarehouse in which the coefficient matrix lives.
     *  @param x The varlabel of the solutio1n vector.
     *  @param modifies_x A boolean that specifies the behaviour of the task 
                          associated with the ScheduleSolve. If set to true,
                          then the task will only modify x. Otherwise, it will
                          compute x. This is a key option when you are computing
                          x in another place.
     * @param b The VarLabel of the right hand side vector.
     * @param which_b_dw Specifies the datawarehouse in which b lives.
     * @param guess VarLabel of the initial guess.
     * @param guess_dw Specifies the datawarehouse of the initial guess.
     * @param params Specifies the solver parameters usually parsed from the input file.
     *
     */    
    virtual void scheduleSolve( const LevelP           & level_in,
                                      SchedulerP       & sched_in,
                                const MaterialSet      * matls_in,
                                const VarLabel         * A_in,
                                      Task::WhichDW      which_A_dw_in,  
                                const VarLabel         * x_in,
                                      bool               modifies_x_in,
                                const VarLabel         * b_in,
                                      Task::WhichDW      which_b_dw_in,  
                                const VarLabel         * guess_in,
                                      Task::WhichDW      which_guess_dw_in,
                                const SolverParameters * params_in,
                                      bool               isFirstSolve_in = true );

    virtual void scheduleInitialize( const LevelP      & level,
                                           SchedulerP  & sched,
                                     const MaterialSet * matls );

    virtual std::string getName();

    void allocateHypreMatrices( DataWarehouse * new_dw );

  private:
    void initialize( const ProcessorGroup *,
                     const PatchSubset    * patches,
                     const MaterialSubset * matls,
                           DataWarehouse  * old_dw,
                           DataWarehouse  * new_dw );

    const VarLabel * m_timeStepLabel;
    const VarLabel * hypre_solver_label;
  };
}

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSolver_h
