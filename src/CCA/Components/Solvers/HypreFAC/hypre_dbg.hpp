#ifndef my_hypre_dbg
#define my_hypre_dbg

inline HYPRE_Int
HYPREDBG_SStructGridCreate ( MPI_Comm           comm,
                             HYPRE_Int          ndim,
                             HYPRE_Int          nparts,
                             HYPRE_SStructGrid * grid )
{
    std::cout << "HYPRE_SStructGridCreate ( comm, " << ndim << ", " << nparts << ", grid )" << std::endl;
    return HYPRE_SStructGridCreate ( comm, ndim, nparts, grid );
}

inline HYPRE_Int
HYPREDBG_SStructGridDestroy ( HYPRE_SStructGrid grid )
{
    std::cout << "HYPRE_SStructGridDestroy ( grid )" << std::endl;
    return HYPRE_SStructGridDestroy ( grid );
}

inline HYPRE_Int
HYPREDBG_SStructGridSetExtents ( HYPRE_SStructGrid  grid,
                                 HYPRE_Int          part,
                                 HYPRE_Int     *    ilower,
                                 HYPRE_Int     *    iupper )
{
#if NDIM == 1
    std::cout << "HYPRE_SStructGridSetExtents ( grid, " << part << ", [" << ilower[0] << "], [" << iupper[0] << "] )" << std::endl;
#elif NDIM == 2
    std::cout << "HYPRE_SStructGridSetExtents ( grid, " << part << ", [" << ilower[0] << "," << ilower[1] << "], [" << iupper[0] << "," << iupper[1] << "] )" << std::endl;
#elif NDIM == 3
    std::cout << "HYPRE_SStructGridSetExtents ( grid, " << part << ", [" << ilower[0] << "," << ilower[1] << "," << ilower[2] << "], [" << iupper[0] << "," << iupper[1] << "," << iupper[2] << "] )" << std::endl;
#endif
    return HYPRE_SStructGridSetExtents ( grid, part, ilower, iupper );
}

inline HYPRE_Int
HYPREDBG_SStructGridSetVariables ( HYPRE_SStructGrid       grid,
                                   HYPRE_Int               part,
                                   HYPRE_Int               nvars,
                                   HYPRE_SStructVariable * vartypes )
{
    std::cout << "HYPRE_SStructGridSetVariables ( grid, " << part << ", " << nvars << ", [";
    for ( int v = 0; v < nvars - 1; ++v )
        std::cout << vartypes[v] << ",";
    std::cout << vartypes[nvars - 1] << "] )" << std::endl;
    return HYPRE_SStructGridSetVariables ( grid, part, nvars, vartypes );
}

inline HYPRE_Int
HYPREDBG_SStructGridAssemble ( HYPRE_SStructGrid grid )
{
    std::cout << "HYPRE_SStructGridAssemble ( grid )" << std::endl;
    return HYPRE_SStructGridAssemble ( grid );
}


inline HYPRE_Int
HYPREDBG_SStructGridSetPeriodic ( HYPRE_SStructGrid grid,
                                  HYPRE_Int         part,
                                  HYPRE_Int       * periodic )
{
#if NDIM == 1
    std::cout << "HYPRE_SStructGridSetPeriodic ( grid, " << part << ", [" << periodic[0] << "] )" << std::endl;
#elif NDIM == 2
    std::cout << "HYPRE_SStructGridSetPeriodic ( grid, " << part << ", [" << periodic[0] << "," << periodic[1] << "] )" << std::endl;
#elif NDIM == 3
    std::cout << "HYPRE_SStructGridSetPeriodic ( grid, " << part << ", [" << periodic[0] << "," << periodic[1] << "," << periodic[3] << "] )" << std::endl;
#endif
    return HYPRE_SStructGridSetPeriodic ( grid, part, periodic );
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

inline HYPRE_Int
HYPREDBG_SStructStencilCreate ( HYPRE_Int             ndim,
                                HYPRE_Int             size,
                                HYPRE_SStructStencil * stencil )
{
    std::cout << "HYPRE_SStructStencilCreate ( " << ndim << ", " << size << ", stencil )" << std::endl;
    return HYPRE_SStructStencilCreate ( ndim, size, stencil );
}

inline HYPRE_Int
HYPREDBG_SStructStencilDestroy ( HYPRE_SStructStencil stencil )
{
    std::cout << "HYPRE_SStructStencilDestroy ( stencil )" << std::endl;
    return HYPRE_SStructStencilDestroy ( stencil );
}

inline HYPRE_Int
HYPREDBG_SStructStencilSetEntry ( HYPRE_SStructStencil  stencil,
                                  HYPRE_Int             entry,
                                  HYPRE_Int      *      offset,
                                  HYPRE_Int             var )
{
#if NDIM == 1
    std::cout << "HYPRE_SStructStencilSetEntry ( stencil, " << entry << ", [" << offset[0] << "], " << var << " )" << std::endl;
#elif NDIM == 2
    std::cout << "HYPRE_SStructStencilSetEntry ( stencil, " << entry << ", [" << offset[0] << "," << offset[1] << "], " << var << " )" << std::endl;
#elif NDIM == 3
    std::cout << "HYPRE_SStructStencilSetEntry ( stencil, " << entry << ", [" << offset[0] << "," << offset[1] << "," << offset[2] << "], " << var << " )" << std::endl;
#endif
    return HYPRE_SStructStencilSetEntry ( stencil, entry, offset, var );
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

inline HYPRE_Int
HYPREDBG_SStructGraphCreate ( MPI_Comm             comm,
                              HYPRE_SStructGrid    grid,
                              HYPRE_SStructGraph * graph )
{
    std::cout << "HYPRE_SStructGraphCreate ( comm, grid, graph )" << std::endl;
    return HYPRE_SStructGraphCreate ( comm, grid, graph );
}

inline HYPRE_Int
HYPREDBG_SStructGraphDestroy ( HYPRE_SStructGraph graph )
{
    std::cout << "HYPRE_SStructGraphDestroy ( graph )" << std::endl;
    return HYPRE_SStructGraphDestroy ( graph );
}

inline HYPRE_Int
HYPREDBG_SStructGraphSetStencil ( HYPRE_SStructGraph   graph,
                                  HYPRE_Int            part,
                                  HYPRE_Int            var,
                                  HYPRE_SStructStencil stencil )
{
    std::cout << "HYPRE_SStructGraphSetStencil ( graph, " << part << ", " << var << ", stencil )" << std::endl;
    return HYPRE_SStructGraphSetStencil ( graph, part, var, stencil );
}

inline HYPRE_Int
HYPREDBG_SStructGraphAddEntries ( HYPRE_SStructGraph   graph,
                                  HYPRE_Int            part,
                                  HYPRE_Int      *     index,
                                  HYPRE_Int            var,
                                  HYPRE_Int            to_part,
                                  HYPRE_Int      *     to_index,
                                  HYPRE_Int            to_var )
{
#if NDIM == 1
    std::cout << "HYPRE_SStructGraphAddEntries ( graph, " << part << ", [" << index[0] << "], " << var << ", " << to_part << ", [" << to_index[0] << "], " << to_var << " )" << std::endl;
#elif NDIM == 2
    std::cout << "HYPRE_SStructGraphAddEntries ( graph, " << part << ", [" << index[0] << "," << index[1] << "], " << var << ", " << to_part << ", [" << to_index[0] << "," << to_index[1] << "], " << to_var << " )" << std::endl;
#elif NDIM == 3
    std::cout << "HYPRE_SStructGraphAddEntries ( graph, " << part << ", [" << index[0] << "," << index[1] << "," << index[2] << "], " << var << ", " << to_part << ", [" << to_index[0] << "," << to_index[1] << "," << to_index[2] << "], " << to_var << " )" << std::endl;
#endif
    return HYPRE_SStructGraphAddEntries ( graph, part, index, var, to_part, to_index, to_var );
}

inline HYPRE_Int
HYPREDBG_SStructGraphAssemble ( HYPRE_SStructGraph graph )
{
    std::cout << "HYPRE_SStructGraphAssemble ( graph )" << std::endl;
    return HYPRE_SStructGraphAssemble ( graph );
}

inline HYPRE_Int
HYPREDBG_SStructGraphSetObjectType ( HYPRE_SStructGraph  graph,
                                     HYPRE_Int           type )
{
    std::cout << "HYPRE_SStructGraphSetObjectType ( graph, " << type << " )" << std::endl;
    return HYPRE_SStructGraphSetObjectType ( graph, type );
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

inline HYPRE_Int
HYPREDBG_SStructMatrixCreate ( MPI_Comm              comm,
                               HYPRE_SStructGraph    graph,
                               HYPRE_SStructMatrix * matrix )
{
    std::cout << "HYPRE_SStructMatrixCreate ( comm, graph, matrix )" << std::endl;
    return HYPRE_SStructMatrixCreate ( comm, graph, matrix );
}

inline HYPRE_Int
HYPREDBG_SStructMatrixDestroy ( HYPRE_SStructMatrix matrix )
{
    std::cout << "HYPRE_SStructMatrixDestroy ( matrix )" << std::endl;
    return HYPRE_SStructMatrixDestroy ( matrix );
}

inline HYPRE_Int
HYPREDBG_SStructMatrixInitialize ( HYPRE_SStructMatrix matrix )
{
    std::cout << "HYPRE_SStructMatrixInitialize ( matrix )" << std::endl;
    return HYPRE_SStructMatrixInitialize ( matrix );
}

inline HYPRE_Int
HYPREDBG_SStructMatrixSetValues ( HYPRE_SStructMatrix  matrix,
                                  HYPRE_Int            part,
                                  HYPRE_Int      *     index,
                                  HYPRE_Int            var,
                                  HYPRE_Int            nentries,
                                  HYPRE_Int      *     entries,
                                  HYPRE_Complex    *   values )
{
#if NDIM == 1
    std::cout << "HYPRE_SStructMatrixSetValues ( matrix, " << part << ", [" << index[0] << "], " << var << ", " << nentries << ", [";
#elif NDIM == 2
    std::cout << "HYPRE_SStructMatrixSetValues ( matrix, " << part << ", [" << index[0] << ","  << index[1] << "], " << var << ", " << nentries << ", [";
#elif NDIM == 3
    std::cout << "HYPRE_SStructMatrixSetValues ( matrix, " << part << ", [" << index[0] << ","  << index[1] << ","  << index[2] << "], " << var << ", " << nentries << ", [";
#endif
    for ( int e = 0; e < nentries - 1; ++e )
        std::cout << entries[e] << ",";
    std::cout << entries[nentries - 1] << "], [";
    for ( int s = 0; s < nentries - 1; ++s )
        std::cout << values[s] << ",";
    std::cout << values[nentries - 1] << "] )" << std::endl;
    return HYPRE_SStructMatrixSetValues ( matrix, part, index, var, nentries, entries, values );
}

inline HYPRE_Int
HYPREDBG_SStructMatrixGetValues ( HYPRE_SStructMatrix  matrix,
                                  HYPRE_Int            part,
                                  HYPRE_Int      *     index,
                                  HYPRE_Int            var,
                                  HYPRE_Int            nentries,
                                  HYPRE_Int      *     entries,
                                  HYPRE_Complex    *   values )
{
    HYPRE_Int res = HYPRE_SStructMatrixGetValues ( matrix, part, index, var, nentries, entries, values );
#if NDIM == 1
    std::cout << "HYPRE_SStructMatrixGetValues ( matrix, " << part << ", [" << index[0] << "], " << var << ", " << nentries << ", [";
#elif NDIM == 2
    std::cout << "HYPRE_SStructMatrixGetValues ( matrix, " << part << ", [" << index[0] << ","  << index[1] << "], " << var << ", " << nentries << ", [";
#elif NDIM == 3
    std::cout << "HYPRE_SStructMatrixGetValues ( matrix, " << part << ", [" << index[0] << ","  << index[1] << ","  << index[2] << "], " << var << ", " << nentries << ", [";
#endif
    for ( int e = 0; e < nentries - 1; ++e )
        std::cout << entries[e] << ",";
    std::cout << entries[nentries - 1] << "], [";
    for ( int s = 0; s < nentries - 1; ++s )
        std::cout << values[s] << ",";
    std::cout << values[nentries - 1] << "] )" << std::endl;
    return res;
}

inline HYPRE_Int
HYPREDBG_SStructMatrixSetBoxValues ( HYPRE_SStructMatrix  matrix,
                                     HYPRE_Int            part,
                                     HYPRE_Int      *     ilower,
                                     HYPRE_Int      *     iupper,
                                     HYPRE_Int            var,
                                     HYPRE_Int            nentries,
                                     HYPRE_Int      *     entries,
                                     HYPRE_Complex    *   values )
{
#if NDIM == 1
    int size = ( iupper[0] - ilower[0] + 1 ) * nentries;
    std::cout << "HYPRE_SStructMatrixSetBoxValues ( matrix, " << part << ", [" << ilower[0] << "], [" << iupper[0] << "], " << var << ", " << nentries << ", [";
#elif NDIM == 2
    int size = ( iupper[0] - ilower[0] + 1 ) * ( iupper[1] - ilower[1] + 1 ) * nentries;
    std::cout << "HYPRE_SStructMatrixSetBoxValues ( matrix, " << part << ", [" << ilower[0] << "," << ilower[1] << "], [" << iupper[0] << "," << iupper[1] << "], " << var << ", " << nentries << ", [";
#elif NDIM == 3
    int size = ( iupper[0] - ilower[0] + 1 ) * ( iupper[1] - ilower[1] + 1 ) * ( iupper[2] - ilower[2] + 1 ) * nentries;
    std::cout << "HYPRE_SStructMatrixSetBoxValues ( matrix, " << part << ", [" << ilower[0] << "," << ilower[1] << "," << ilower[2] << "], [" << iupper[0] << "," << iupper[1] << "," << iupper[2] << "], " << var << ", " << nentries << ", [";
#endif
    for ( int e = 0; e < nentries - 1; ++e )
        std::cout << entries[e] << ",";
    std::cout << entries[nentries - 1] << "], [";
    for ( int s = 0; s < size - 1; ++s )
        std::cout << values[s] << ",";
    std::cout << values[size - 1] << "] )" << std::endl;
    return HYPRE_SStructMatrixSetBoxValues ( matrix, part, ilower, iupper, var, nentries, entries, values );
}

inline HYPRE_Int
HYPREDBG_SStructMatrixGetBoxValues ( HYPRE_SStructMatrix  matrix,
                                     HYPRE_Int            part,
                                     HYPRE_Int      *     ilower,
                                     HYPRE_Int      *     iupper,
                                     HYPRE_Int            var,
                                     HYPRE_Int            nentries,
                                     HYPRE_Int      *     entries,
                                     HYPRE_Complex    *   values )
{
    HYPRE_Int res = HYPRE_SStructMatrixGetBoxValues ( matrix, part, ilower, iupper, var, nentries, entries, values );
#if NDIM == 1
    int size = ( iupper[0] - ilower[0] + 1 ) * nentries;
    std::cout << "HYPRE_SStructMatrixGetBoxValues ( matrix, " << part << ", [" << ilower[0] << "] , [" << iupper[0] << "], " << var << ", " << nentries << ", [";
#elif NDIM == 2
    int size = ( iupper[0] - ilower[0] + 1 ) * ( iupper[1] - ilower[1] + 1 ) * nentries;
    std::cout << "HYPRE_SStructMatrixGetBoxValues ( matrix, " << part << ", [" << ilower[0] << "," << ilower[1] << "] , [" << iupper[0] << "," << iupper[1] << "], " << var << ", " << nentries << ", [";
#elif NDIM == 3
    int size = ( iupper[0] - ilower[0] + 1 ) * ( iupper[1] - ilower[1] + 1 ) * ( iupper[2] - ilower[2] + 1 ) * nentries;
    std::cout << "HYPRE_SStructMatrixGetBoxValues ( matrix, " << part << ", [" << ilower[0] << "," << ilower[1] << "," << ilower[2] << "] , [" << iupper[0] << "," << iupper[1] << "," << iupper[2] << "], " << var << ", " << nentries << ", [";
#endif
    for ( int e = 0; e < nentries - 1; ++e )
        std::cout << entries[e] << ",";
    std::cout << entries[nentries - 1] << "], [";
    for ( int s = 0; s < size - 1; ++s )
        std::cout << values[s] << ",";
    std::cout << values[size - 1] << "] )" << std::endl;
    return res;
}

inline HYPRE_Int
HYPREDBG_SStructMatrixAssemble ( HYPRE_SStructMatrix matrix )
{
    std::cout << "HYPRE_SStructMatrixAssemble ( matrix )" << std::endl;
    return HYPRE_SStructMatrixAssemble ( matrix );
}

inline HYPRE_Int
HYPREDBG_SStructMatrixSetObjectType ( HYPRE_SStructMatrix  matrix,
                                      HYPRE_Int            type )
{
    std::cout << "HYPRE_SStructMatrixSetObjectType ( matrix, " << type  << " )" << std::endl;
    return HYPRE_SStructMatrixSetObjectType ( matrix, type );
}

inline HYPRE_Int
HYPREDBG_SStructVectorCreate ( MPI_Comm              comm,
                               HYPRE_SStructGrid     grid,
                               HYPRE_SStructVector * vector )
{
    std::cout << "HYPRE_SStructVectorCreate ( comm, grid, vector )" << std::endl;
    return HYPRE_SStructVectorCreate ( comm, grid, vector );
};

inline HYPRE_Int
HYPREDBG_SStructVectorDestroy ( HYPRE_SStructVector vector )
{
    std::cout << "HYPRE_SStructVectorDestroy ( vector )" << std::endl;
    return HYPRE_SStructVectorDestroy ( vector );
}

inline HYPRE_Int
HYPREDBG_SStructVectorInitialize ( HYPRE_SStructVector vector )
{
    std::cout << "HYPRE_SStructVectorInitialize ( vector )" << std::endl;
    return HYPRE_SStructVectorInitialize ( vector );
}

inline HYPRE_Int
HYPREDBG_SStructVectorSetValues ( HYPRE_SStructVector  vector,
                                  HYPRE_Int            part,
                                  HYPRE_Int      *     index,
                                  HYPRE_Int            var,
                                  HYPRE_Complex    *   value );

inline HYPRE_Int
HYPREDBG_SStructVectorGetValues ( HYPRE_SStructVector  vector,
                                  HYPRE_Int            part,
                                  HYPRE_Int      *     index,
                                  HYPRE_Int            var,
                                  HYPRE_Complex    *   value )
{
    HYPRE_Int res = HYPRE_SStructVectorGetValues ( vector, part, index, var, value );
#if NDIM == 1
    std::cout << "HYPREDBG_SStructVectorGetValues ( vector, " << part << ", [" << index[0] << "], " << var << ", " << *value << " )" << std::endl;
#elif NDIM == 2
    std::cout << "HYPREDBG_SStructVectorGetValues ( vector, " << part << ", [" << index[0] << ","  << index[1] << "], " << var << ", "  << *value << " )" << std::endl;
#elif NDIM == 3
    std::cout << "HYPREDBG_SStructVectorGetValues ( vector, " << part << ", [" << index[0] << ","  << index[1] << ","  << index[2] << "], " << var << ", " << *value << " )" << std::endl;
#endif
    return res;
}

inline HYPRE_Int
HYPREDBG_SStructVectorSetBoxValues ( HYPRE_SStructVector  vector,
                                     HYPRE_Int            part,
                                     HYPRE_Int      *     ilower,
                                     HYPRE_Int      *     iupper,
                                     HYPRE_Int            var,
                                     HYPRE_Complex    *   values )
{
#if NDIM == 1
    int size = ( iupper[0] - ilower[0] + 1 );
    std::cout << "HYPRE_SStructVectorSetBoxValues ( vector, " << part << ", [" << ilower[0] << "], [" << iupper[0] << "], " << var << ", [" ;
#elif NDIM == 2
    int size = ( iupper[0] - ilower[0] + 1 ) * ( iupper[1] - ilower[1] + 1 );
    std::cout << "HYPRE_SStructVectorSetBoxValues ( vector, " << part << ", [" << ilower[0] << "," << ilower[1] << "], [" << iupper[0] << "," << iupper[1] << "], " << var << ", [" ;
#elif NDIM == 3
    int size = ( iupper[0] - ilower[0] + 1 ) * ( iupper[1] - ilower[1] + 1 ) * ( iupper[2] - ilower[2] + 1 );
    std::cout << "HYPRE_SStructVectorSetBoxValues ( vector, " << part << ", [" << ilower[0] << "," << ilower[1] << "," << ilower[2] << "], [" << iupper[0] << "," << iupper[1] << "," << iupper[2] << "], " << var << ", [" ;
#endif
    for ( int s = 0; s < size - 1; ++s )
        std::cout << values[s] << ",";
    std::cout << values[size - 1] << "] )" << std::endl;
    return HYPRE_SStructVectorSetBoxValues ( vector, part, ilower, iupper, var, values );
}

inline HYPRE_Int
HYPREDBG_SStructVectorGetBoxValues ( HYPRE_SStructVector  vector,
                                     HYPRE_Int            part,
                                     HYPRE_Int      *     ilower,
                                     HYPRE_Int      *     iupper,
                                     HYPRE_Int            var,
                                     HYPRE_Complex    *   values )
{
    HYPRE_Int res = HYPRE_SStructVectorGetBoxValues ( vector, part, ilower, iupper, var, values );
#if NDIM == 1
    int size = ( iupper[0] - ilower[0] + 1 );
    std::cout << "HYPRE_SStructVectorGetBoxValues ( vector, " << part << ", [" << ilower[0] << "] , [" << iupper[0] << "] , " << 0 << ", [";
#elif NDIM == 2
    int size = ( iupper[0] - ilower[0] + 1 ) * ( iupper[1] - ilower[1] + 1 );
    std::cout << "HYPRE_SStructVectorGetBoxValues ( vector, " << part << ", [" << ilower[0] << "," << ilower[1] << "] , [" << iupper[0] << "," << iupper[1] << "] , " << 0 << ", [";
#elif NDIM == 3
    int size = ( iupper[0] - ilower[0] + 1 ) * ( iupper[1] - ilower[1] + 1 ) * ( iupper[2] - ilower[2] + 1 );
    std::cout << "HYPRE_SStructVectorGetBoxValues ( vector, " << part << ", [" << ilower[0] << "," << ilower[1] << "," << ilower[2] << "] , [" << iupper[0] << "," << iupper[1] << "," << iupper[2] << "] , " << 0 << ", [";
#endif
    for ( int s = 0; s < size - 1; ++s )
        std::cout << values[s] << ",";
    std::cout << values[size - 1] << "] )" << std::endl;
    return res;
}

inline HYPRE_Int
HYPREDBG_SStructVectorAssemble ( HYPRE_SStructVector vector )
{
    std::cout << "HYPRE_SStructVectorAssemble ( vector )" << std::endl;
    return HYPRE_SStructVectorAssemble ( vector );
}

inline HYPRE_Int
HYPREDBG_SStructVectorGather ( HYPRE_SStructVector vector )
{
    std::cout << "HYPRE_SStructVectorGather ( vector )" << std::endl;
    return HYPRE_SStructVectorGather ( vector );
}

inline HYPRE_Int
HYPREDBG_SStructVectorSetObjectType ( HYPRE_SStructVector  vector,
                                      HYPRE_Int            type )
{
    std::cout << "HYPRE_SStructVectorSetObjectType ( vector, " << type << " )" << std::endl;
    return HYPRE_SStructVectorSetObjectType ( vector, type );
}

/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------*/

inline HYPRE_Int
HYPREDBG_SStructFACCreate ( MPI_Comm             comm,
                            HYPRE_SStructSolver * solver )
{
    std::cout << "HYPRE_SStructFACCreate ( comm, solver )" << std::endl;
    return HYPRE_SStructFACCreate ( comm, solver );
}

inline HYPRE_Int
HYPREDBG_SStructFACDestroy2 ( HYPRE_SStructSolver solver )
{
    std::cout << "HYPRE_SStructFACDestroy2 ( solver )" << std::endl;
    return HYPRE_SStructFACDestroy2 ( solver );
}

inline HYPRE_Int
HYPREDBG_SStructFACSetup2 ( HYPRE_SStructSolver solver,
                            HYPRE_SStructMatrix A,
                            HYPRE_SStructVector b,
                            HYPRE_SStructVector x )
{
    std::cout << "HYPRE_SStructFACSetup2 ( solver, A, b, x );" << std::endl;
    return HYPRE_SStructFACSetup2 ( solver, A, b, x );
}


inline HYPRE_Int
HYPREDBG_SStructFACSolve3 ( HYPRE_SStructSolver solver,
                            HYPRE_SStructMatrix A,
                            HYPRE_SStructVector b,
                            HYPRE_SStructVector x )
{
    std::cout << "HYPRE_SStructFACSolve3 ( solver, A, b, x );" << std::endl;
    return HYPRE_SStructFACSolve3 ( solver, A, b, x );
}

inline HYPRE_Int
HYPREDBG_SStructFACSetPLevels ( HYPRE_SStructSolver solver,
                                HYPRE_Int           nparts,
                                HYPRE_Int     *     plevels )
{
    std::cout << "HYPRE_SStructFACSetPLevels ( solver, " << nparts << ", [";
    for ( int s = 0; s < nparts - 1; ++s )
        std::cout << plevels[s] << ",";
    std::cout << plevels[nparts - 1] << "] )" << std::endl;
    return HYPRE_SStructFACSetPLevels ( solver, nparts, plevels );
}

inline HYPRE_Int
HYPREDBG_SStructFACSetPRefinements ( HYPRE_SStructSolver  solver,
                                     HYPRE_Int            nparts,
                                     HYPRE_Int ( *rfactors ) [HYPRE_MAXDIM] )
{
    std::cout << "HYPRE_SStructFACSetPRefinements ( solver, " << nparts << ", [";
    for ( int s = 0; s < nparts - 1; ++s )
        std::cout << "[" << rfactors[s][0] << "," << rfactors[s][1] << "," << rfactors[s][2] << "],";
    std::cout << "[" << rfactors[nparts - 1][0] << "," << rfactors[nparts - 1][1] << "," << rfactors[nparts - 1][2] << "]] ) " << std::endl;
    return HYPRE_SStructFACSetPRefinements ( solver, nparts, rfactors );
}

inline HYPRE_Int
HYPREDBG_SStructFACZeroCFSten ( HYPRE_SStructMatrix  A,
                                HYPRE_SStructGrid    grid,
                                HYPRE_Int            part,
                                HYPRE_Int            rfactors[HYPRE_MAXDIM] )
{
    std::cout << "HYPRE_SStructFACZeroCFSten ( A, grid, " << part << ", [" << rfactors[0] << "," << rfactors[1] << "," << rfactors[2] << "] )" << std::endl;
    return HYPRE_SStructFACZeroCFSten ( A, grid, part, rfactors );
}

inline HYPRE_Int
HYPREDBG_SStructFACZeroFCSten ( HYPRE_SStructMatrix  A,
                                HYPRE_SStructGrid    grid,
                                HYPRE_Int            part )
{
    std::cout << "HYPRE_SStructFACZeroFCSten ( A, grid, " << part << " )" << std::endl;
    return HYPRE_SStructFACZeroFCSten ( A, grid, part );
}

inline HYPRE_Int
HYPREDBG_SStructFACZeroAMRMatrixData ( HYPRE_SStructMatrix  A,
                                       HYPRE_Int            part_crse,
                                       HYPRE_Int            rfactors[HYPRE_MAXDIM] )
{
    std::cout << "HYPRE_SStructFACZeroAMRMatrixData ( A, " << part_crse << ", [" << rfactors[0] << "," << rfactors[1] << "," << rfactors[2] << "] )" << std::endl;
    return HYPRE_SStructFACZeroAMRMatrixData ( A, part_crse, rfactors );
}

#define HYPREDBG_SStructFACZeroAMRVectorData(_1,_2,_3) hypredbg_SStructFACZeroAMRVectorData (_1, nparts, _2, _3)

inline HYPRE_Int
hypredbg_SStructFACZeroAMRVectorData ( HYPRE_SStructVector  b,
                                       HYPRE_Int            nparts,
                                       HYPRE_Int      *     plevels,
                                       HYPRE_Int ( *rfactors ) [HYPRE_MAXDIM] )
{
    std::cout << "HYPRE_SStructFACZeroAMRVectorData ( b, [";
    for ( int s = 0; s < nparts - 1; ++s )
        std::cout << plevels[s] << ",";
    std::cout << plevels[nparts - 1] << "], [";
    for ( int s = 0; s < nparts - 1; ++s )
        std::cout << "[" << rfactors[s][0] << "," << rfactors[s][1] << "," << rfactors[s][2] << "],";
    std::cout << "[" << rfactors[nparts - 1][0] << "," << rfactors[nparts - 1][1] << "," << rfactors[nparts - 1][2] << "]] ) " << std::endl;
    return HYPRE_SStructFACZeroAMRVectorData ( b, plevels, rfactors );

}


inline HYPRE_Int
HYPREDBG_SStructFACSetMaxLevels ( HYPRE_SStructSolver solver ,
                                  HYPRE_Int           max_levels )
{
    std::cout << "HYPRE_SStructFACSetMaxLevels ( solver, " << max_levels << " )" << std::endl;
    return HYPRE_SStructFACSetMaxLevels ( solver, max_levels );
}

inline HYPRE_Int
HYPREDBG_SStructFACSetTol ( HYPRE_SStructSolver solver,
                            HYPRE_Real          tol )
{
    std::cout << "HYPRE_SStructFACSetTol ( solver, " << tol << " );" << std::endl;
    return HYPRE_SStructFACSetTol ( solver, tol );
}

inline HYPRE_Int
HYPREDBG_SStructFACSetMaxIter ( HYPRE_SStructSolver solver,
                                HYPRE_Int           max_iter )
{
    std::cout << "HYPRE_SStructFACSetMaxIter ( solver, " << max_iter << " )" << std::endl;
    return HYPRE_SStructFACSetMaxIter ( solver, max_iter );
}

inline HYPRE_Int
HYPREDBG_SStructFACSetRelChange ( HYPRE_SStructSolver solver,
                                  HYPRE_Int           rel_change )
{
    std::cout << "HYPRE_SStructFACSetRelChange ( solver, " << rel_change << " )" << std::endl;
    return HYPRE_SStructFACSetRelChange ( solver, rel_change );
}

inline HYPRE_Int
HYPREDBG_SStructFACSetZeroGuess ( HYPRE_SStructSolver solver )
{
    std::cout << "HYPRE_SStructFACSetZeroGuess ( solver  )" << std::endl;
    return HYPRE_SStructFACSetZeroGuess ( solver );
}

inline HYPRE_Int
HYPREDBG_SStructFACSetNonZeroGuess ( HYPRE_SStructSolver solver )
{
    std::cout << "HYPRE_SStructFACSetNonZeroGuess ( solver )" << std::endl;
    return HYPRE_SStructFACSetNonZeroGuess ( solver );
}

inline HYPRE_Int
HYPREDBG_SStructFACSetRelaxType ( HYPRE_SStructSolver solver,
                                  HYPRE_Int           relax_type )
{
    std::cout << "HYPRE_SStructFACSetRelaxType ( solver, " << relax_type << " )" << std::endl;
    return HYPRE_SStructFACSetRelaxType ( solver, relax_type );
}

inline HYPRE_Int
HYPREDBG_SStructFACSetJacobiWeight ( HYPRE_SStructSolver solver,
                                     HYPRE_Real          weight );

inline HYPRE_Int
HYPREDBG_SStructFACSetNumPreRelax ( HYPRE_SStructSolver solver,
                                    HYPRE_Int           num_pre_relax )
{
    std::cout << "HYPRE_SStructFACSetNumPreRelax ( solver, " << num_pre_relax << " )" << std::endl;
    return HYPRE_SStructFACSetNumPreRelax ( solver, num_pre_relax );
}

inline HYPRE_Int
HYPREDBG_SStructFACSetNumPostRelax ( HYPRE_SStructSolver solver,
                                     HYPRE_Int           num_post_relax )
{
    std::cout << "HYPRE_SStructFACSetNumPostRelax ( solver, " << num_post_relax << " )" << std::endl;
    return HYPRE_SStructFACSetNumPostRelax ( solver, num_post_relax );
}

inline HYPRE_Int
HYPREDBG_SStructFACSetCoarseSolverType ( HYPRE_SStructSolver solver,
        HYPRE_Int           csolver_type )
{
    std::cout << "HYPRE_SStructFACSetCoarseSolverType ( solver, " << csolver_type << " )" << std::endl;
    return HYPRE_SStructFACSetCoarseSolverType ( solver, csolver_type );
}

inline HYPRE_Int
HYPREDBG_SStructFACSetLogging ( HYPRE_SStructSolver solver,
                                HYPRE_Int           logging )
{
    std::cout << "HYPRE_SStructFACSetLogging ( solver, " << logging << " )" << std::endl;
    return HYPRE_SStructFACSetLogging ( solver, logging );
}

inline HYPRE_Int
HYPREDBG_SStructFACGetNumIterations ( HYPRE_SStructSolver  solver,
                                      HYPRE_Int      *     num_iterations )
{
    HYPRE_Int res = HYPRE_SStructFACGetNumIterations ( solver, num_iterations );
    std::cout << "HYPRE_SStructFACGetNumIterations ( solver, " << *num_iterations << " )" << std::endl;
    return res;
}

inline HYPRE_Int
HYPREDBG_SStructFACGetFinalRelativeResidualNorm ( HYPRE_SStructSolver solver,
                                                  HYPRE_Real     *    norm )
{
    HYPRE_Int res = HYPRE_SStructFACGetFinalRelativeResidualNorm ( solver, norm );
    std::cout << "HYPRE_SStructFACGetFinalRelativeResidualNorm ( solver, " << *norm << " )" << std::endl;
    return res;
}

#endif



