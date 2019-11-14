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

#ifndef Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructImplementation_h
#define Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructImplementation_h

#include <CCA/Components/Solvers/HypreSStruct/MatrixIndex.h>
#include <CCA/Components/Solvers/HypreSStruct/AdditionalEntries.h>
#include <CCA/Components/Solvers/HypreSStruct/GlobalData.h>
#include <CCA/Components/Solvers/HypreSStruct/PartData.h>
#include <CCA/Components/Solvers/HypreSStruct/ExtraEntries.h>
#include <CCA/Components/Solvers/HypreSStruct/SStructInterface.h>

#include <HYPRE_sstruct_ls.h>

#include <numeric>

#if 1
#define HYPRE(fn) HYPRE_##fn
#else
#define HYPREDBG_NDIM 2
#include "hypre_dbg.hpp"
#define HYPRE(fn) HYPREDBG_##fn
#endif


namespace Uintah
{
namespace HypreSStruct
{

template < int DIM, int C2F >
class SStructImplementation
    : public SStructInterface
{
protected:

    GlobalDataP gdata;  // global data  `
    PartDataP * pdatas; // local data per part/level
    int  * plevels;
    int ( *prefinements ) [HYPRE_MAXDIM];

    bool                   m_restart;

    HYPRE_SStructGrid      grid;
    HYPRE_SStructStencil   stencil;
    HYPRE_SStructGraph     graph;
    HYPRE_SStructSolver    solver;
    HYPRE_SStructMatrix    A;
    HYPRE_SStructVector    b;
    HYPRE_SStructVector    x;
    bool                   grid_initialized, stencil_initialized, graph_initialized;
    bool                   solver_initialized;
    bool                   A_initialized, b_initialized, x_initialized;
    bool                   guess_updated;

    ExtraEntries<StnIndex, DIM> ** extra_stn_entries;
    ExtraEntries<AddIndex, DIM> ** extra_add_entries;

protected:

    SStructImplementation (
        const GlobalDataP & global_data
    ) : gdata ( global_data ),
        pdatas ( nullptr ),
        plevels ( nullptr ),
        prefinements ( nullptr ),
        m_restart ( true ),
        grid_initialized ( false ),
        stencil_initialized ( false ),
        graph_initialized ( false ),
        solver_initialized ( false ),
        A_initialized ( false ),
        b_initialized ( false ),
        x_initialized ( false ),
        guess_updated ( false ),
        extra_stn_entries ( nullptr ),
        extra_add_entries ( nullptr )
    {
        const int & nvars = gdata->nVars();
        const int & nparts = gdata->nParts();

        pdatas = scinew PartDataP[nparts];
        plevels = scinew int[nparts];
        prefinements = scinew int[nparts][HYPRE_MAXDIM];

        if ( nparts > 1 )
        {
            extra_stn_entries = scinew ExtraEntries<StnIndex, DIM> * [nvars];
            extra_add_entries = scinew ExtraEntries<AddIndex, DIM> * [nvars];
            for ( int i = 0; i < nvars; ++i )
            {
                extra_stn_entries[i] = scinew ExtraEntries<StnIndex, DIM> [ nparts ];
                extra_add_entries[i] = scinew ExtraEntries<AddIndex, DIM> [ nparts ];
            }
        }
    }

    virtual ~SStructImplementation()
    {
        if ( x_initialized ) HYPRE ( SStructVectorDestroy ) ( x );
        if ( b_initialized ) HYPRE ( SStructVectorDestroy ) ( b );
        if ( A_initialized ) HYPRE ( SStructMatrixDestroy ) ( A );
        if ( graph_initialized ) HYPRE ( SStructGraphDestroy ) ( graph );
        if ( stencil_initialized ) HYPRE ( SStructStencilDestroy ) ( stencil );
        if ( grid_initialized ) HYPRE ( SStructGridDestroy ) ( grid );

        if ( extra_stn_entries )
        {
            for ( int i = 0; i < gdata->nVars(); ++i )
            {
                delete[] extra_stn_entries[i];
                delete[] extra_add_entries[i];
            }
            delete[] extra_stn_entries;
            delete[] extra_add_entries;
        }

        delete[] pdatas;
        delete[] plevels;
        delete[] prefinements;
    }

    virtual void solverFinalize() = 0;

    /// connection from
    void
    add_c2f_refinement_entry
    (
        const int & var,
        const int & coarse_part,
        const int & coarse_box,
        const IntVector & coarse_index,
        const IntVector & fine_index
    )
    {
        StnIndex stn { coarse_box, coarse_index, 0 };
        AddIndex ind { -1, coarse_index, coarse_part + 1, fine_index };

        int * index    = std::get<1> ( ind ).get_pointer();
        int & to_part  = std::get<2> ( ind );
        int * to_index = std::get<3> ( ind ).get_pointer();

        if ( extra_add_entries[var][coarse_part].emplace_back ( ind, stn, 0. ) )
            HYPRE ( SStructGraphAddEntries ) ( graph, coarse_part, index, var, to_part, to_index, var );
    }

    template<int C = C2F, typename std::enable_if<C == 0, int>::type = 0>
    void
    add_c2f_interface_entry
    (
        const int & var,
        const int & coarse_part,
        const int & coarse_box,
        const IntVector & coarse_index,
        const int & coarse_rank,
        const IntVector & fine_index,
        const int & face_dir,
        const int & face_sgn
    )
    {
        StnIndex stn { coarse_box, coarse_index, coarse_rank };
        AddIndex ind { -1, coarse_index, coarse_part + 1, fine_index };

        int * index    = std::get<1> ( ind ).get_pointer();
        int & to_part  = std::get<2> ( ind );
        int * to_index = std::get<3> ( ind ).get_pointer();

        if ( extra_add_entries[var][coarse_part].emplace_back ( ind, stn, 1. ) )
            HYPRE ( SStructGraphAddEntries ) ( graph, coarse_part, index, var, to_part, to_index, var );

        extra_stn_entries[var][coarse_part].emplace_back ( stn, stn, 0 );
    }

    template<int C = C2F, typename std::enable_if<C == 1, int>::type = 0>
    void
    add_c2f_interface_entry
    (
        const int & var,
        const int & coarse_part,
        const int & coarse_box,
        const IntVector & coarse_index,
        const int & coarse_rank,
        const IntVector & fine_index,
        const int & face_dir,
        const int & face_sgn
    )
    {
        StnIndex stn { coarse_box, coarse_index, coarse_rank };
        AddIndex ind { -1, coarse_index, coarse_part + 1, fine_index };

        int * index    = std::get<1> ( ind ).get_pointer();
        int & to_part  = std::get<2> ( ind );
        int * to_index = std::get<3> ( ind ).get_pointer();

        int ( & refinement ) [HYPRE_MAXDIM] = prefinements[to_part];
        const double & weight = pdatas[to_part]->pRefinementWeight();

        int lower[3] = { fine_index[0], fine_index[1], fine_index[2] };
        if ( face_sgn > 0 ) lower[face_dir] -= refinement[face_dir] - 1;
        int upper[3] = { lower[0] + refinement[0], lower[1] + refinement[1], lower[2] + refinement[2] };

        for ( to_index[2] = lower[2]; to_index[2] < upper[2]; ++to_index[2] )
            for ( to_index[1] = lower[1]; to_index[1] < upper[1]; ++to_index[1] )
                for ( to_index[0] = lower[0]; to_index[0] < upper[0]; ++to_index[0] )
                {
                    if ( extra_add_entries[var][coarse_part].emplace_back ( ind, stn, weight ) )
                        HYPRE ( SStructGraphAddEntries ) ( graph, coarse_part, index, var, to_part, to_index, var );
                }
        extra_stn_entries[var][coarse_part].emplace_back ( stn, stn, 0 );
    }

    void
    add_f2f_interface_entry
    (
        const int & var,
        const AddIndex & add,
        const IntVector & refined_index,
        const double & weight
    )
    {
        AddIndex ind { -1, std::get<1> ( add ), std::get<2> ( add ), refined_index };

        int * index    = std::get<1> ( ind ).get_pointer();
        int & part     = std::get<2> ( ind );
        int * to_index = std::get<3> ( ind ).get_pointer();

        // check if stencil entry
        for ( int rank = 0; rank < SStructStencil<DIM>::size; ++rank )
        {
            bool is_stencil = true;
            for ( int d = 0; d < DIM; ++d )
            {
                if ( to_index[d] - index[d] != SStructStencil<DIM>::offsets[rank][d] )
                {
                    is_stencil = false;
                    break;
                }
            }

            if ( is_stencil )
            {
                StnIndex stn { std::get<0> (add), refined_index, rank };
                // add entry to existing stencil
                extra_stn_entries[var][part].emplace_back ( stn, stn, 1. );
                extra_stn_entries[var][part].emplace_back ( stn, add, weight );
                return;
            }
        }

        // otherwise is additional entry
        if ( extra_add_entries[var][part].emplace_back ( ind, add, weight ) )
            HYPRE ( SStructGraphAddEntries ) ( graph, part, index, var, part, to_index, var );
    }

    template<int C = C2F, typename std::enable_if<C == 0, int>::type = 0>
    void
    add_f2f_interface_entry
    (
        const int & var,
        const int & fine_part,
        const int & fine_box,
        const IntVector & fine_index,
        const IntVector & coarse_index,
        const IntVector & refined_index
    )
    {
        AddIndex add { fine_box, fine_index, fine_part - 1, coarse_index };
        add_f2f_interface_entry ( var, add, refined_index, 1. );
    }

    template<int C = C2F, typename std::enable_if<C == 1, int>::type = 0>
    void
    add_f2f_interface_entry
    (
        const int & var,
        const int & fine_part,
        const int & fine_box,
        const IntVector & fine_index,
        const IntVector & coarse_index,
        const IntVector & refined_index
    )
    {
        AddIndex add { fine_box, fine_index, fine_part - 1, coarse_index };

        const int ( & refinement ) [HYPRE_MAXDIM] = prefinements[fine_part];
        const double & weight = pdatas[fine_part]->pRefinementWeight();

        int lower[3] = { refined_index[0], refined_index[1], refined_index[2] };
        int upper[3] = { lower[0] + refinement[0], lower[1] + refinement[1], lower[2] + refinement[2] };

        IntVector to_index;
        for ( to_index[2] = lower[2]; to_index[2] < upper[2]; ++to_index[2] )
            for ( to_index[1] = lower[1]; to_index[1] < upper[1]; ++to_index[1] )
                for ( to_index[0] = lower[0]; to_index[0] < upper[0]; ++to_index[0] )
                    add_f2f_interface_entry ( var, add, to_index, weight );
    }

    virtual void
    add_f2c_interface_entry
    (
        const int & var,
        const int & fine_part,
        const int & fine_box,
        const IntVector & fine_index,
        const IntVector & coarse_index
    )
    {
        AddIndex add { fine_box, fine_index, fine_part - 1, coarse_index };

        int * index    = std::get<1> ( add ).get_pointer();
        int & to_part  = std::get<2> ( add );
        int * to_index = std::get<3> ( add ).get_pointer();

        if ( extra_add_entries[var][fine_part].emplace_back ( add, add, 1. ) )
            HYPRE ( SStructGraphAddEntries ) ( graph, fine_part, index, var, to_part, to_index, var );
    };

    std::vector<double>
    get_stencil_values (
        const PartDataP & pdata,
        const int & box,
        const CCVariable<Stencil7> & stencil_entries
    )
    {
        std::vector<double> values ( pdata->boxSize ( box ) * SStructStencil<DIM>::size );

        std::vector<double>::iterator it = values.begin();
        for ( int k = pdata->iLowerD ( box, 2 ); k <= pdata->iUpperD ( box, 2 ); ++k )
            for ( int j = pdata->iLowerD ( box, 1 ); j <= pdata->iUpperD ( box, 1 ); ++j )
                for ( int i = pdata->iLowerD ( box, 0 ); i <= pdata->iUpperD ( box, 0 ); ++i )
                {
                    *it++ = stencil_entries[IntVector ( i, j, k )].p;
                    for ( int e = 0; e < 2 * DIM; ++e )
                        *it++ = stencil_entries[IntVector ( i, j, k )][e];
                }
        return values;
    }

public:

    virtual const GlobalDataP &
    globalData (
    ) const override
    {
        return gdata;
    }

    virtual const PartDataP &
    partData (
        const int & part
    ) const override
    {
        return pdatas[part];
    }

    virtual bool
    isPartSet (
        const int & part
    ) const override
    {
        return pdatas[part];
    }

    virtual void
    setPart (
        const int & part,
        const int & npatches,
        const int & nboxes,
        const int low[HYPRE_MAXDIM],
        const int high[HYPRE_MAXDIM],
        const int refinement[HYPRE_MAXDIM],
        const int periodic[HYPRE_MAXDIM]
    ) override
    {
        gdata->setPart ( part, npatches );
        pdatas[part] = scinew PartData ( nboxes, low, high, refinement, periodic );
        plevels[part] = part;
        std::copy_n ( refinement, HYPRE_MAXDIM, prefinements[part] );
    }

    virtual void
    addBox (
        const int & part,
        const int & id,
        const int ilower[HYPRE_MAXDIM],
        const int iupper[HYPRE_MAXDIM],
        const bool interface[2 * HYPRE_MAXDIM]
    ) override
    {
        pdatas[part]->addBox ( id, ilower, iupper, interface );
    }

    virtual void
    gridInitialize (
        const MPI_Comm & comm
    ) override
    {
        ASSERT ( !grid_initialized );

        HYPRE ( SStructGridCreate ) ( comm, DIM, gdata->nParts(), &grid );
        for ( int part = 0; part < gdata->nParts(); ++part )
        {
            for ( int box = 0; box < gdata->nBoxes ( part ); ++box )
                HYPRE ( SStructGridSetExtents ) ( grid, part, const_cast<int *> ( pdatas[part]->iLower ( box ) ), const_cast<int *> ( pdatas[part]->iUpper ( box ) ) );
            HYPRE ( SStructGridSetVariables ) ( grid, part, gdata->nVars(), const_cast<HYPRE_SStructVariable *> ( gdata->varTypes() ) );
            HYPRE ( SStructGridSetPeriodic ) ( grid, part, const_cast<int *> ( pdatas[part]->periodic() ) );
        }
        HYPRE ( SStructGridAssemble ) ( grid );
        grid_initialized = true;
    }

    virtual void
    stencilInitialize (
        const MPI_Comm & /*comm*/
    ) override
    {
        ASSERT ( !stencil_initialized );
        HYPRE ( SStructStencilCreate ) ( DIM, SStructStencil<DIM>::size, &stencil );
        for ( int var = 0; var < gdata->nVars(); ++var )
            for ( int entry = 0; entry < SStructStencil<DIM>::size; ++entry )
                HYPRE ( SStructStencilSetEntry ) ( stencil, entry, SStructStencil<DIM>::offsets[entry], var );

        stencil_initialized = true;
    }

    virtual void
    graphInitialize (
        const MPI_Comm & comm,
        const GridP & grd,
        AdditionalEntries ** ** additional_entries
    ) override
    {
        ASSERT ( !graph_initialized );
        HYPRE ( SStructGraphCreate ) ( comm, grid, &graph );
        HYPRE ( SStructGraphSetObjectType ) ( graph, HYPRE_SSTRUCT );

        /* set stencils */
        for ( int part = 0; part < gdata->nParts(); ++part )
            for ( int var = 0; var < gdata->nVars(); ++var )
                HYPRE ( SStructGraphSetStencil ) ( graph, part, var, stencil );

        if ( additional_entries )
        {
            /* add coarse to fine entries */
            for ( int fine_part = 1; fine_part < gdata->nParts(); ++fine_part )
            {
                int coarse_part = fine_part - 1;
                const int ( & refinement ) [HYPRE_MAXDIM] = prefinements[fine_part];

                /* add zero connection to fine boxes from the corresponding
                 * forming the coarse operators since it reassigns coarse boxs
                 * to the processors owning their refined fine boxes */
                for ( int coarse_box = 0; coarse_box < gdata->nBoxes ( coarse_part ); ++coarse_box )
                {
                    const Patch * coarse_patch = grd->getPatchByID ( pdatas[coarse_part]->patch ( coarse_box ), coarse_part );
                    const Level * coarse_level = coarse_patch->getLevel();
                    const Level * fine_level = coarse_level->getFinerLevel().get_rep();

                    // compute refined range for coarse_box
                    IntVector fine_lower, fine_upper;
                    for ( size_t d = 0; d < HYPRE_MAXDIM; ++d )
                    {
                        fine_lower[d] = pdatas[coarse_part]->iLowerD ( coarse_box, d ) * refinement[d];
                        fine_upper[d] = ( pdatas[coarse_part]->iUpperD ( coarse_box, d ) + 1 ) * refinement[d];
                    }

                    // select fine_patches in refined coarse_box range
                    std::vector<const Patch *> fine_patches;
                    fine_level->selectPatches ( fine_lower, fine_upper, fine_patches );

                    for ( const Patch * fine_patch : fine_patches )
                    {
                        // create connection to fine_patch low fine_index from
                        // its corresponding coarse_index
                        IntVector fine_index = fine_patch->getCellLowIndex();
                        IntVector coarse_index;
                        for ( size_t d = 0; d < HYPRE_MAXDIM; ++d )
                            coarse_index[d] = fine_index[d] / refinement[d] - ( fine_index[d] < 0 );

                        for ( int var = 0; var < gdata->nVars(); ++var )
                            add_c2f_refinement_entry ( var, coarse_part, coarse_box, coarse_index, fine_index );
                    }

                    /* add addtitional coarse to fine entries */

                    // stretch refined range for coarse_box to include ghosts
                    for ( int d = 0; d < DIM; ++d )
                    {
                        fine_lower[d] -= refinement[d];
                        fine_upper[d] += refinement[d];
                    }

                    // select fine_patches in refined coarse_box range
                    fine_level->selectPatches ( fine_lower, fine_upper, fine_patches );

                    for ( const Patch * fine_patch : fine_patches )
                    {
                        IntVector fine_lower = fine_patch->getCellLowIndex();
                        IntVector fine_upper = fine_patch->getCellHighIndex();

                        // compute refined coarse region corresponding to fine_patch
                        IntVector coarse_lower, coarse_upper;
                        for ( size_t d = 0; d < HYPRE_MAXDIM; ++d )
                        {
                            fine_upper[d] -= 1;
                            coarse_lower[d] = fine_lower[d] / refinement[d] - ( fine_lower[d] < 0 );
                            coarse_upper[d] = fine_upper[d] / refinement[d] - ( fine_upper[d] < 0 );
                        }

                        for ( int face_dir = 0; face_dir < DIM; ++face_dir )
                            for ( int face_sgn :
                                    {
                                        -1, 1
                                    } )
                            {
                                int face = 2 * face_dir + ( face_sgn + 1 ) / 2;

                                // check if face is an F2C interface
                                if ( fine_patch->getBCType ( ( Patch::FaceType ) face ) == Patch::Coarse )
                                {
                                    // stencil rank that points to refined coarse
                                    // indices corresponding to face
                                    int coarse_entry = face - face_sgn;

                                    // compute face coarse and fine ranges
                                    IntVector ilower = coarse_lower;
                                    IntVector iupper = coarse_upper;
                                    if ( fine_patch->isVirtual() )
                                    {
                                        fine_lower = fine_patch->getRealPatch()->getCellLowIndex();
                                        fine_upper = fine_patch->getRealPatch()->getCellHighIndex();
                                    }
                                    else
                                    {
                                        fine_lower = fine_patch->getCellLowIndex();
                                        fine_upper = fine_patch->getCellHighIndex();
                                    }
                                    if ( face_sgn < 0 )
                                    {
                                        iupper[face_dir] = ilower[face_dir] -= 1;
                                    }
                                    else
                                    {
                                        ilower[face_dir] = iupper[face_dir] += 1;
                                        fine_lower[face_dir] = fine_upper[face_dir] - refinement[face_dir] + 1;
                                    }

                                    for ( int d = 0; d < DIM; ++d )
                                    {
                                        if ( ilower[d] < pdatas[coarse_part]->iLowerD ( coarse_box, d ) )
                                            ilower[d] = pdatas[coarse_part]->iLowerD ( coarse_box, d );
                                        if ( iupper[d] > pdatas[coarse_part]->iUpperD ( coarse_box, d ) )
                                            iupper[d] = pdatas[coarse_part]->iUpperD ( coarse_box, d );
                                    }

                                    // loop over face ranges
                                    IntVector coarse_index, fine_index;
                                    for ( coarse_index[0] = ilower[0], fine_index[0] = fine_lower[0]; coarse_index[0] <= iupper[0]; ++coarse_index[0], fine_index[0] += refinement[0] )
                                        for ( coarse_index[1] = ilower[1], fine_index[1] = fine_lower[1]; coarse_index[1] <= iupper[1]; ++coarse_index[1], fine_index[1] += refinement[1] )
                                            for ( coarse_index[2] = ilower[2], fine_index[2] = fine_lower[2]; coarse_index[2] <= iupper[2]; ++coarse_index[2], fine_index[2] += refinement[2] )
                                                for ( int var = 0; var < gdata->nVars(); ++var )
                                                    add_c2f_interface_entry ( var, coarse_part, coarse_box, coarse_index, coarse_entry + 1, fine_index, face_dir, face_sgn );
                                }
                            }
                    }
                }

                /* add fine to coarse entries */

                for ( int var = 0; var < gdata->nVars(); ++var )
                    for ( int fine_box = 0; fine_box < gdata->nBoxes ( fine_part ); ++fine_box )
                    {
                        for ( const auto & entry : *additional_entries[var][fine_part][fine_box] )
                        {
                            IntVector fine_index = std::get<0> ( entry.first );
                            const int & coarse_part = std::get<1> ( entry.first );
                            IntVector coarse_index = std::get<2> ( entry.first );

                            // check if coarse index in entry is refined
                            bool is_refined = false;
                            if ( fine_part - coarse_part == 1 )
                            {
                                IntVector to_fine_index;
                                for ( size_t i = 0; i < HYPRE_MAXDIM; ++i ) // till last dimension otherwise getPatchFromIndex fails!
                                    to_fine_index[i] = coarse_index[i] * refinement[i];

                                for ( int to_fine_box = 0; to_fine_box < gdata->nBoxes ( fine_part ); ++to_fine_box )
                                {
                                    for ( int d = 0; d < DIM; ++d )
                                    {
                                        is_refined = true;
                                        if ( to_fine_index[d] < pdatas[fine_part]->iLowerD ( to_fine_box, d ) ||
                                                to_fine_index[d] > pdatas[fine_part]->iUpperD ( to_fine_box, d ) )
                                        {
                                            is_refined = false;
                                            break;
                                        }
                                    }

                                    if ( is_refined )
                                    {
                                        add_f2f_interface_entry ( var, fine_part, fine_box, fine_index, coarse_index, to_fine_index );
                                        break;
                                    }
                                }
                            }

                            if ( !is_refined )
                            {
                                for ( int var = 0; var < gdata->nVars(); ++var )
                                    add_f2c_interface_entry ( var, fine_part, fine_box, fine_index, coarse_index );
                            }
                        }
                    }
            }
        }

        HYPRE ( SStructGraphAssemble ) ( graph );
        graph_initialized = true;
    }

    virtual void
    matrixInitialize (
        const MPI_Comm & comm
    ) override
    {
        ASSERT ( !A_initialized );
        HYPRE ( SStructMatrixCreate ) ( comm, graph, &A );
        HYPRE ( SStructMatrixSetObjectType ) ( A, HYPRE_SSTRUCT );
        HYPRE ( SStructMatrixInitialize ) ( A );
        A_initialized = true;
    }

    virtual void
    rhsInitialize (
        const MPI_Comm & comm
    ) override
    {
        ASSERT ( !b_initialized );
        HYPRE ( SStructVectorCreate ) ( comm, grid, &b );
        HYPRE ( SStructVectorSetObjectType ) ( b, HYPRE_SSTRUCT );
        HYPRE ( SStructVectorInitialize ) ( b );
        b_initialized = true;
    }

    virtual void
    solutionInitialize (
        const MPI_Comm & comm
    ) override
    {
        ASSERT ( !x_initialized );
        HYPRE ( SStructVectorCreate ) ( comm, grid, &x );
        HYPRE ( SStructVectorSetObjectType ) ( x, HYPRE_SSTRUCT );
        HYPRE ( SStructVectorInitialize ) ( x );
        x_initialized = true;
    }

    virtual void
    matrixUpdate (
        constCCVariable<Stencil7> ** * stencil_entries,
        AdditionalEntries ** ** additional_entries
    ) override
    {
        ASSERT ( A_initialized );

        for ( int part = 0; part < gdata->nParts(); ++part )
        {
            auto & pdata = pdatas[part];
            for ( int box = 0; box < gdata->nBoxes ( part ); ++box )
            {
                std::vector<int> stencil_indices ( SStructStencil<DIM>::size );
                std::iota ( std::begin ( stencil_indices ), std::end ( stencil_indices ), 0 );

                for ( int var = 0; var < gdata->nVars(); ++var )
                {
                    std::vector<double> values = get_stencil_values ( pdata, box, stencil_entries[var][part][box] );
                    HYPRE ( SStructMatrixSetBoxValues ) (
                        A, part, const_cast<int *> ( pdata->iLower ( box ) ), const_cast<int *> ( pdata->iUpper ( box ) ),
                        var, SStructStencil<DIM>::size, stencil_indices.data(), values.data() );
                }
            }
        }

        if ( extra_stn_entries && extra_add_entries )
        {
            for ( int var = 0; var < gdata->nVars(); ++var )
            {
                for ( int part = 0; part < gdata->nParts(); ++part )
                {
                    for ( const auto & extra : extra_stn_entries[var][part] ( stencil_entries[var][part], additional_entries[var][part] ) )
                        HYPRE ( SStructMatrixSetValues ) ( A, part, extra.index(), var, extra.nValues(), extra.entries(), extra.values() );
                    for ( const auto & extra : extra_add_entries[var][part] ( stencil_entries[var][part], additional_entries[var][part] ) )
                        HYPRE ( SStructMatrixSetValues ) ( A, part, extra.index(), var, extra.nValues(), extra.entries(), extra.values() );
                }
            }
        }

        if ( gdata->nParts() > 1 )
        {
            for ( int part = gdata->nParts() - 1; part > 0; part-- )
            {
                // THIS SHOULD WORK FOR ALL SOLVERS NOT ONLY FAC
                HYPRE ( SStructFACZeroCFSten ) ( A, grid, part, prefinements[part] );
                HYPRE ( SStructFACZeroFCSten ) ( A, grid, part );
                HYPRE ( SStructFACZeroAMRMatrixData ) ( A, part - 1, prefinements[part] );
            }
        }
    }

    virtual void
    rhsUpdate (
        constCCVariable<double> ** * rhs
    ) override
    {
        ASSERT ( b_initialized );

        for ( int part = 0; part < gdata->nParts(); ++part )
        {
            auto & pdata = pdatas[part];
            for ( int box = 0; box < gdata->nBoxes ( part ); ++box )
                for ( int k = pdata->iLowerD ( box, 2 ); k <= pdata->iUpperD ( box, 2 ); ++k )
                    for ( int j = pdata->iLowerD ( box, 1 ); j <= pdata->iUpperD ( box, 1 ); ++j )
                    {
                        IntVector ll ( pdata->iLowerD ( box, 0 ), j, k );
                        IntVector hh ( pdata->iUpperD ( box, 0 ), j, k );
                        for ( int var = 0; var < gdata->nVars(); ++var )
                        {
                            const double * vals = &rhs[var][part][box][IntVector ( pdata->iLowerD ( box, 0 ), j, k )];
                            HYPRE ( SStructVectorSetBoxValues ) ( b, part, ll.get_pointer(), hh.get_pointer(), var, const_cast<double *> ( vals ) );
                        }
                    }
        }
        if ( gdata->nParts() > 1 )
            HYPRE ( SStructFACZeroAMRVectorData ) ( b, plevels, prefinements );
    }

    virtual void
    guessUpdate (
        constCCVariable<double> ** * guess
    ) override
    {
        ASSERT ( x_initialized );
        guess_updated = true;

        for ( int part = 0; part < gdata->nParts(); ++part )
        {
            auto & pdata = pdatas[part];
            for ( int box = 0; box < gdata->nBoxes ( part ); ++box )
                for ( int k = pdata->iLowerD ( box, 2 ); k <= pdata->iUpperD ( box, 2 ); ++k )
                    for ( int j = pdata->iLowerD ( box, 1 ); j <= pdata->iUpperD ( box, 1 ); ++j )
                    {
                        IntVector ll ( pdata->iLowerD ( box, 0 ), j, k );
                        IntVector hh ( pdata->iUpperD ( box, 0 ), j, k );
                        for ( int var = 0; var < gdata->nVars(); ++var )
                        {
                            const double * vals = &guess[var][part][box][IntVector ( pdata->iLowerD ( box, 0 ), j, k )];
                            HYPRE ( SStructVectorSetBoxValues ) ( x, part, ll.get_pointer(), hh.get_pointer(), var, const_cast<double *> ( vals ) );
                        }
                    }
        }
        if ( gdata->nParts() > 1 )
            HYPRE ( SStructFACZeroAMRVectorData ) ( x, plevels, prefinements );
    }

    virtual void
    assemble (
    ) override
    {
        HYPRE ( SStructMatrixAssemble ) ( A );
        HYPRE ( SStructVectorAssemble ) ( b );
        HYPRE ( SStructVectorAssemble ) ( x );
    }

#ifdef PRINTSYSTEM
    virtual void
    printSystem (
        std::string * fname
    ) override
    {
        HYPRE_SStructMatrixPrint ( fname[0].c_str(), A, 0 );
        HYPRE_SStructVectorPrint ( fname[1].c_str(), b, 0 );
        HYPRE_SStructVectorPrint ( fname[2].c_str(), x, 0 );
    }
#endif

    virtual void
    getSolution (
        CCVariable<double> ** * solution
    ) override
    {
        HYPRE ( SStructVectorGather ) ( x );

        for ( int part = 0; part < gdata->nParts(); ++part )
        {
            auto & pdata = pdatas[part];
            for ( int box = 0; box < gdata->nBoxes ( part ); ++box )
                for ( int k = pdata->iLowerD ( box, 2 ); k <= pdata->iUpperD ( box, 2 ); ++k )
                    for ( int j = pdata->iLowerD ( box, 1 ); j <= pdata->iUpperD ( box, 1 ); ++j )
                    {
                        IntVector ll ( pdata->iLowerD ( box, 0 ), j, k );
                        IntVector hh ( pdata->iUpperD ( box, 0 ), j, k );
                        for ( int var = 0; var < gdata->nVars(); ++var )
                        {
                            double * vals = &solution[var][part][box][IntVector ( pdata->iLowerD ( box, 0 ), j, k )];
                            HYPRE ( SStructVectorGetBoxValues ) ( x, part, ll.get_pointer(), hh.get_pointer(), var, vals );
                        }
                    }
        }
    }

    virtual void
    finalize (
    ) override
    {
        if ( x_initialized ) HYPRE ( SStructVectorDestroy ) ( x );
        if ( b_initialized ) HYPRE ( SStructVectorDestroy ) ( b );
        if ( A_initialized ) HYPRE ( SStructMatrixDestroy ) ( A );
        if ( graph_initialized ) HYPRE ( SStructGraphDestroy ) ( graph );
        if ( stencil_initialized ) HYPRE ( SStructStencilDestroy ) ( stencil );
        if ( grid_initialized ) HYPRE ( SStructGridDestroy ) ( grid );
        solverFinalize();
        x_initialized = b_initialized = A_initialized = graph_initialized = stencil_initialized = grid_initialized = false;
    }

    virtual bool
    restart (
    ) override
    {
        if ( !m_restart ) return false;

        finalize();

        m_restart = false;
        return true;
    }
};

} // namespace HypreSStruct
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_HypreSStruct_SStructInterface_h




