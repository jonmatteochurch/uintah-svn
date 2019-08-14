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

/**
 * @file CCA/Components/PhaseField/DataTypes/SubProblems.h
 * @author Jon Matteo Church [j.m.church@leeds.ac.uk]
 * @date 2018/12
 */

#ifndef Packages_Uintah_CCA_Components_PhaseField_DataTypes_SubProblems_h
#define Packages_Uintah_CCA_Components_PhaseField_DataTypes_SubProblems_h

#include <CCA/Components/PhaseField/DataTypes/Problem.h>
#include <CCA/Components/PhaseField/BoundaryConditions/BCInterface.h>
#include <CCA/Components/PhaseField/Applications/Application.h>
// #include <numeric>

namespace Uintah
{
namespace PhaseField
{

/// Debug stream
extern Dout g_subproblems_dbg;

template < typename Problem, bool AMR > class Application;

/**
 * @brief PhaseField Problem Container Variable
 *
 * Stores the list of all Problem's in which a Patch is partitioned
 *
 * @tparam Problem type of PhaseField problem
 */
template < typename Problem > class SubProblems;

/**
 * @brief PhaseField Problem Container Variable (Problem implementation)
 *
 * @implements SubProblems < Problem >
 *
 * @tparam VAR type of variable representation
 * @tparam STN finite-difference stencil
 * @tparam Field list of type of fields (ScalarField < T > or VectorField < T, N >)
 */
template<VarType VAR, StnType STN, typename... Field>
class SubProblems < Problem<VAR, STN, Field...> >
    : public SubProblemsVariableBase

{
public: // STATIC MEMBERS

    /// Problem Dimension
    static constexpr DimType DIM = get_stn<STN>::dim;

    /// Problem number of faces
    static constexpr size_t NF = 2 * DIM;

    /// SubProblems TypeDescription::Type value
    static constexpr TypeDescription::Type TD = (VAR==CC) ? TypeDescription::CCSubProblems : TypeDescription::NCSubProblems;

private: // TYPES

    /**
     * @brief Container of data for MPI communications
     *
     * @remark SubProblems are designed to hold instances for accessing the DW
     * which need to be recreated on receiver. Only required data for their
     * instantiation is sent
     */
    struct MpiData
            : public RefCounted
    {
    private: // STATIC MEMBERS

        /// Register
        /// @remark Cause the (instantiated) type of this class to be registered with the Core/Disclosure/TypeDescription class when this class.
        /// @remark It is not used for anything else in the
        // program.
        static TypeDescription::Register registerMe;

        /// Unique instance of TypeDescription
        static TypeDescription * td;

    public: // MEMBERS

        /// face is flagged if any bc is applied to it
        std::array<bool, NF> m_flags;

        /// arrays of BC info for each face
        std::tuple < std::array< BCInfo<Field>, NF >... > m_bcs;

    public: // CONSTRUCTORS

        /**
         * @brief Default Constructor
         *
         * @remark Should be used only by SubProblems::allocate an empty stuct
         * for being filled with received MPI buffer
         */
        MpiData () = default;

        /**
         * @brief Constructor
         *
         * @param bcs list of arrays of BC info for each one of face of the patch for each one of the problem labels
         * @param flags array of flags to check if any bc is applied to each one of faces
         */
        MpiData (
            std::array<bool, NF> flags,
            std::array< BCInfo<Field>, NF>... bcs
        ) : m_flags ( flags ), m_bcs ( bcs... )
        {}

    private: // STATIC METHODS

        /**
         * @brief Get the offset of the SubProblems flag member
         *
         * @remark SubProblems is not POD and offsetof cannot be used
         *
         * @return the offset, in bytes, from the beginning of each object
         */
        static MPI_Aint
        get_flags_offset()
        {
            MpiData instance;
            MPI_Aint base, flags;
            MPI_Get_address ( &instance, &base );
            MPI_Get_address ( instance.m_flags.data(), &flags );
            return flags - base;
        }

        /**
         * @brief Get the offset of one the SubProblems bcs members
         *
         * @remark SubProblems is not POD and offsetof cannot be used
         *
         * @tparam I index of the bcs
         * @return the offset, in bytes, from the beginning of each object
         */
        template <size_t I>
        static MPI_Aint
        get_bcs_offset()
        {
            MpiData instance;
            MPI_Aint base, bcs;
            MPI_Get_address ( &instance, &base );
            MPI_Get_address ( std::get<I> ( instance.m_bcs ).data(), &bcs );
            return bcs - base;
        }

        /**
         * @brief Create and commit MPI_Datatype for BCInfo
         *
         * Define a new struct datatype for communicating BCInfo with MPI and commit it

         * @remark the reference to this private method is passed to the TypeDescription
         * constructor to enure a single instance of MPI_Datatype is created
         *
         * @tparam I list of boundary fields' indices
         * @param unused to allow template argument deduction
         * @return MPI_Datatype for the BCInfo struct instance
         */
        template <size_t... I>
        static MPI_Datatype
        make_mpitype (
            index_sequence<I...> _DOXYARG ( unused )
        )
        {
            MPI_Datatype datatype;
            static constexpr int count = 1 + sizeof... ( Field );

            MPI_Datatype types[count] =
            {
                PhaseField::getTypeDescription<bool>()->getMPIType(),
                PhaseField::getTypeDescription<BCInfo<Field>>()->getMPIType()...
            };

            int block_length[count] =
            {
                nfaces<bool>::value,
                nfaces<Field>::value...
            };

            MPI_Aint displacement[count] =
            {
                get_flags_offset(),
                get_bcs_offset<I>()...
            };

            Uintah::MPI::Type_create_struct ( count, block_length, displacement, types, &datatype );
            MPI::Type_commit ( &datatype );
            return datatype;
        }

        /**
         * @brief Create and commit MPI_Datatype for BCInfo
         *
         * Define a new struct datatype for communicating BCInfo with MPI and commit it

         * @remark the reference to this private method is passed to the TypeDescription
         * constructor to enure a single instance of MPI_Datatype is created
         *
         * @return MPI_Datatype for the BCInfo struct instance
         */
        static MPI_Datatype
        make_mpitype()
        {
            return make_mpitype ( index_sequence_for<Field...> {} );
        }

    public: // STATIC METHODS

        /**
         * @brief Get the TypeDescription for BCInfo
         *
         * Get the pointer to the single instance of TypeDescription for the BCInfo
         * type instance. If this is not yet initialized then a new instance is created.
         *
         * @return Pointer to the TypeDescription for the BCInfo struct instance
         */
        static const TypeDescription *
        getTypeDescription()
        {
            if ( !td )
            {
                td = scinew TypeDescription ( TypeDescription::Other, "SubProblems::MpiData", true, MpiData::make_mpitype );
            }
            return td;
        }
    };

    /**
     * @brief Container of Problem instances
     */
    struct ProblemList
            : public RefCounted
            , public std::list < Problem<VAR, STN, Field...> >
    {};

private: // STATIC MEMBERS

    /// Register
    /// @remark Cause the (instantiated) type of this class to be registered with the Core/Disclosure/TypeDescription class when this class.
    /// @remark It is not used for anything else in the
    // program.
    static TypeDescription::Register registerMe;

    /// Unique instance of TypeDescription
    static TypeDescription * td;

private: // MEMBERS

    /// Info required by addGhostRegion to Instantiate new SubProblems when MPI message is received
    const std::tuple< typename Field::label_type... > * m_boundary_labels; ///< List of labels for boundary variables
    const VarLabel * m_subproblems_label; ///< Label for subproblems in the DataWarehouse
    const int * m_material; ///< Index of material in the DataWarehouse
    const Patch * m_patch; ///< Grid patch where to instantiate subproblems

    /// List populated by addGhostRegion and required by setGhostRegion after MPI message is received
    std::list<Region> m_ghost_regions;

    /// Data for MPI communications
    Handle < MpiData > m_mpi_data;

    /// List of instantiated Problem's
    Handle < ProblemList > m_list;

private: // CONSTRUCTORS/DESTRUCTOR

    /**
     * @brief Default Constructor
     *
     * @remark used only by make for receiving MPI and by the DataWarehouse
     * contructor to create an empty instance to be populated
     */
    SubProblems (
    ) : m_boundary_labels ( nullptr ),
        m_subproblems_label ( nullptr ),
        m_material ( nullptr ),
        m_patch ( nullptr ),
        m_ghost_regions (),
        m_mpi_data ( nullptr ),
        m_list ( nullptr )
    {
    }

    /**
     * @brief Inner Patch Constructor
     *
     * Partition the given patch using the given flags array and bcs infos
     *
     * @param labels list of lables for each variable of the problem
     * @param subproblems_label label for subproblems in the DataWarehouse
     * @param material index in the DataWarehouse
     * @param patch grid patch
     * @param flags array of flags to check if any bc is applied to each one of faces
     * @param bcs list of arrays of BC info for each one of face of the patch for each one of the problem labels
     */
    SubProblems (
        const typename Field::label_type & ...  labels,
        const VarLabel * subproblems_label,
        int material,
        const Patch * patch,
        const std::array<bool, 2 * DIM> & flags,
        const std::array< BCInfo<Field>, get_dim<DIM>::face_end > & ... bcs
    ) : m_boundary_labels ( nullptr ),
        m_subproblems_label ( nullptr ),
        m_material ( nullptr ),
        m_patch ( nullptr ),
        m_ghost_regions (),
        m_mpi_data ( scinew MpiData ( flags, bcs... ) ),
        m_list ( scinew ProblemList () )
    {
        auto problems = BCInterface<VAR, STN>::template partition_patch <Field...> ( labels..., subproblems_label, material, patch, bcs..., flags );
        m_list->splice ( m_list->begin(), problems );

        // CHECK if patch->d_patchSate can replace flags
        for ( size_t f = 0; f < 2 * DIM; ++f )
        {
            bool flag = flags.at ( f );
            Patch::BCType state = patch->getBCType ( ( Patch::FaceType ) f );
            bool comp_flag = state == Patch::Coarse || state == Patch::None;
            std::stringstream msg;
            msg << "cannot replace flag at " << *patch << " " << Patch::getFaceName ( ( Patch::FaceType ) f ) << " flag (" << flag << ") state " << state << "(" << comp_flag << ")";
            ASSERTMSG ( flag == comp_flag, msg.str().c_str() );
        }
    }

    /**
     * @brief Middle Patch Constructor
     *
     * Call BCInterface::get_bcs to populate the given flags array and retrieve
     * the BCInfos required to call the Inner Patch Constructor
     *
     * @param labels list of lables for each variable of the problem
     * @param subproblems_label label for subproblems in the DataWarehouse
     * @param material index in the DataWarehouse
     * @param patch grid patch
     * @param c2f mapping between variable names and C2F condition for amr grids
     * @param[in,out] flags array of flags to check if any bc is applied to each one of faces
     */
    SubProblems (
        const typename Field::label_type & ...  labels,
        const VarLabel * subproblems_label,
        int material,
        const Patch * patch,
        const std::map<std::string, FC> * c2f,
        std::array<bool, 2 * DIM> flags
    ) : SubProblems ( labels..., subproblems_label, material, patch, flags, BCInterface<VAR, STN>::template get_bcs<Field> ( patch, 0, material, labels, flags, c2f )... )
    {
    }

    /**
     * @brief Copy constructor
     *
     * @param copy instance to be copied to this
     */
    SubProblems (
        const SubProblems & copy
    ) : m_boundary_labels ( copy.m_boundary_labels ),
        m_subproblems_label ( copy.m_subproblems_label ),
        m_material ( copy.m_material ? scinew int ( *copy.m_material ) : nullptr ),
        m_patch ( copy.m_patch ),
        m_ghost_regions ( copy.m_ghost_regions ),
        m_mpi_data ( copy.m_mpi_data ),
        m_list ( copy.m_list )
    {
        if ( copy.isForeign() ) this->setForeign();
        if ( copy.isValid() ) this->setValid();
    }

    /**
     * @brief Assignment operator
     *
     * @param copy instance to be copied to this
     * @return reference to this
     */
    SubProblems & operator= ( const SubProblems & copy )
    {
        m_boundary_labels = copy.m_boundary_labels;
        m_subproblems_label = copy.m_subproblems_label;
        m_material = copy.m_material ? scinew int ( *copy.m_material ) : nullptr;
        m_patch = copy.m_patch;
        m_ghost_regions = copy.m_ghost_regions;
        m_mpi_data = copy.m_mpi_data;
        m_list = copy.m_list;

        if ( copy.isForeign() ) this->setForeign();
        if ( copy.isValid() ) this->setValid();

        return *this;
    }

public:

    /**
     * @brief Patch Constructor
     *
     * Instantiate the list the views required to handle all given variables
     * over the given patch without retrieving data.
     *
     * @remark There is no need to recreate the views until the geometry is unchanged
     *
     * @param labels list of lables for each variable of the problem
     * @param subproblems_label label for subproblems in the DataWarehouse
     * @param material index in the DataWarehouse
     * @param patch grid patch
     * @param c2f mapping between variable names and C2F condition for amr grids
     */
    SubProblems (
        const typename Field::label_type & ...  labels,
        const VarLabel * subproblems_label,
        int material,
        const Patch * patch,
        const std::map<std::string, FC> * c2f = nullptr
    ) : SubProblems ( labels..., subproblems_label, material, patch, c2f, make_flags() )
    {}

    /**
     * @brief Constructor
     *
     * Retrieve the list the views required to handle all given variables
     * over the given patch without retrieving data.
     *
     * @remark There is no need to recreate the views until the geometry is unchanged
     *
     * @param dw DataWarehouse from which to retrieve the subproblem list
     * @param subproblems_label label for subproblems in the DataWarehouse
     * @param material index in the DataWarehouse
     * @param patch grid patch
     */
    SubProblems (
        DataWarehouse * dw,
        const VarLabel * subproblems_label,
        int material,
        const Patch * patch
    ) : SubProblems()
    {
        dw->get ( *this, subproblems_label, material, patch );
    };

    /**
     * @brief Destructor
     */
    ~SubProblems ()
    {
        delete m_material;
    }

private: // STATIC METHODS

    /**
     * @brief Create initial flags array
     *
     * @return default face flags
     */
    static std::array<bool, 2 * DIM> make_flags ()
    {
        std::array<bool, 2 * DIM> flags;
        flags.fill ( false );
        return flags;
    }

    /**
     * @brief Create an empty SubProblems instance
     *
     * @return new instance
     */
    // used in TypeDescription::createInstance() when OnDemandDataWarehouse::recvMPI needs to create an empty variable to fill with received buffer
    static Variable *
    maker()
    {
        return scinew SubProblems < Problem<VAR, STN, Field...> > ();
    }

public: // STATIC METHODS

    /**
    * @brief Get SubProblems TypeDescription
    *
    * @return corresponding TypeDescription sungle instance
    */
    static const TypeDescription *
    getTypeDescription()
    {
        if ( !td )
        {
            td = scinew TypeDescription ( TD, "SubProblems", &maker, MpiData::getTypeDescription() );
        }
        return td;
    }

private:

    /// Helper struct for the number of faces
    template<typename> using nfaces = std::integral_constant<int, 2 * DIM>;

    /**
     * @brief Get BufferInfo (internal indexed implementation)
     *
     * Retrieve the info required to serialize inner data into an MPI message
     *
     * @tparam I list of boundary fields' indices
     * @param unused to allow template argument deduction
     * @param[out] buffer requested information
     */
    template<size_t... I>
    void
    getMPIBuffer (
        index_sequence<I...> _DOXYARG ( unused ),
        BufferInfo & buffer
    ) const
    {
        MPI_Datatype datatype = getMPIType<MpiData>();
        void * startbuf = getBasePointer();
        buffer.add ( startbuf, 1, datatype, false );
    }

    /**
     * @brief Add to the list the Problem over the given ghost region
     *
     * @tparam I list of boundary field' indices
     * @param unused to allow template argument deduction
     * @param low lower bound of the ghost region to partition
     * @param high higher bound of the ghost region to partition
     */
    template<size_t... I>
    void
    setGhostRegion (
        index_sequence<I...> _DOXYARG ( unused ),
        const IntVector & low,
        const IntVector & high

    )
    {
        // this should by set when OnDemandDataWarehouse::recvMPI calls allocate
        ASSERT ( m_boundary_labels );

        DOUTR ( g_subproblems_dbg, "processing subproblem on ghost region " << low << high << " from patch " << *m_patch );
        auto problems = BCInterface<VAR, STN>::template partition_patch <Field...> ( std::get<I> ( *m_boundary_labels )..., m_subproblems_label, *m_material, m_patch, std::get<I> ( m_mpi_data->m_bcs )..., m_mpi_data->m_flags );
        auto p1 = problems.begin();
        while ( p1 != problems.end() )
        {
            if ( p1->restrict_range( low, high ) )
            {
                ++p1;
            }
            else
            {
                DOUTR ( g_subproblems_dbg,  "   not adding problem " << *p1 << " because out of ghost region limits" );
                p1 = problems.erase ( p1 );
            }
        }

        if ( m_list )
        {
            if ( !m_list->empty() )
            {
                auto p1 = problems.begin();
                auto p0 = m_list->begin();
                while ( p1 != problems.end() )
                {
                    if ( p0 == m_list->end() )
                    {
                        ++p1;
                        p0 = m_list->begin();
                        continue;
                    }
                    if ( p0->get_low() <= p1->get_low() && p1->get_high() <= p0->get_high() )
                    {
                        // no need to add to list
                        DOUTR ( g_subproblems_dbg,  "   not adding problem " << *p1 << " because " << *p0 << " already added" );
                        p1 = problems.erase ( p1 );
                        p0 = m_list->begin();
                    }
                    else if ( p1->get_low() <= p0->get_low() && p0->get_high() <= p1->get_high() )
                    {
                        // no need to add to list
                        DOUTR ( g_subproblems_dbg,  "   removing problem " << *p0 << " because " << *p1 << " is being added" );
                        p0 = m_list->erase ( p0 );
                    }
                    else
                    {
                        ++p0;
                    }
                }
            }
            if ( !problems.empty() )
            {
                DOUTR ( g_subproblems_dbg, "added more problems" );
                for ( const auto & p : problems )
                    DOUTR ( g_subproblems_dbg, "   added problem " << p );

                m_list->splice ( m_list->begin(), problems );
            }
        }
        else
        {
            m_list = scinew ProblemList ();
            m_list->splice ( m_list->begin(), problems );

            DOUTR ( g_subproblems_dbg, "added problems" );
            for ( const auto & p : *this )
                DOUTR ( g_subproblems_dbg, "   added problem " << p );
        }
    }

public: // VARIABLE METHODS

    /**
     * @brief Finizalize the initialization of the SubProblems list after an MPI
     * message has been received
     *
     * Once information about flags and bcs is received it is possible to create
     * the problems for handling the ghost regions specified through addGhostRegion
     * by OnDemandDataWarehouse::recvMPI
     *
     * @return If the initialization was successful
     */
    virtual bool
    onMPIReceived()
    override
    {
        for ( const auto & region : m_ghost_regions )
            setGhostRegion ( index_sequence_for<Field...> {}, region.getLow(), region.getHigh() );
        m_ghost_regions.clear();
        return true;
    }

    /**
    * @brief Get SubProblems TypeDescription
    *
    * @return corresponding TypeDescription sungle instance
    */
    virtual const TypeDescription *
    virtualGetTypeDescription()
    const override
    {
        return getTypeDescription();
    }

    /**
    * @brief Serialize SubProblems
    *
    * @remark Not Implemented
    *
    * @param out output stream
    * @param l low grid index of the region to output
    * @param h high grid index of the region to output
    * @param varnode xml variable node
    * @param outputDoubleAsFloat output precision
    */
    virtual void
    emitNormal (
        std::ostream & _DOXYARG ( out ),
        const IntVector & _DOXYARG ( l ),
        const IntVector & _DOXYARG ( h ),
        ProblemSpecP _DOXYARG ( varnode ),
        bool _DOXYARG ( outputDoubleAsFloat )
    )
    {
        SCI_THROW ( InternalError ( "Cannot yet write SubProblems!\n", __FILE__, __LINE__ ) );
    }

    /**
    * @brief Deserialize SubProblems
    *
    * @remark Not Implemented
    *
    * @param in input stream
    * @param swapBytes handle endiannes
    */
    virtual void
    readNormal (
        std::istream & _DOXYARG ( in ),
        bool _DOXYARG ( swapBytes )
    )
    {
        SCI_THROW ( InternalError ( "Cannot yet read SubProblems!\n", __FILE__, __LINE__ ) );
    }

    /**
     * @brief Get Info about variable memory usage
     *
     * @remark MPIData is the inner data
     *
     * @param[out] elems # of elements
     * @param[out] totsize # of bytes
     * @param[out] ptr pointer to the inner data
     */
    virtual void
    getSizeInfo (
        std::string & elems,
        unsigned long & totsize,
        void *& ptr
    ) const override
    {
        elems = "1";
        totsize = getDataSize();
        ptr = getBasePointer();
    }

    /**
     * @brief Set size info of the underlying data
     *
     * @remark this is for host-->device variable copy
     * @return inner data byte size
     */
    virtual size_t
    getDataSize()
    const override
    {
        return sizeof ( MpiData );
    }

    /**
     * @brief Byte copy the underlying data
     *
     * @param dst pointer to the memory location to copy to
     * @return if the copy is successful
     */
    virtual bool
    copyOut (
        void * dst
    ) const override
    {
        void * src = getBasePointer();
        size_t numBytes = getDataSize();
        void * retVal = std::memcpy ( dst, src, numBytes );
        return ( retVal == dst ) ? true : false;
    }

    /**
     * @brief Ponter copy
     *
     * Copy this instance to the one pointed by the given virtual base class reference
     *
     * @param copy reference to Variable to copy to
     */
    virtual void
    copyPointer (
        Variable & copy
    ) override
    {
        auto * c = dynamic_cast<const SubProblems < Problem<VAR, STN, Field...> > * > ( &copy );
        if ( !c )
        {
            SCI_THROW ( TypeMismatchException ( "Type mismatch in SubProblems variable", __FILE__, __LINE__ ) );
        }
        *this = *c;
    };

    /**
     * @brief Get the RefCounted interface
     *
     * @return pointer to inner data
     */
    virtual RefCounted *
    getRefCounted()
    override
    {
        return m_mpi_data.get_rep();
    };

public: // SUBPROBLEMSVARIABLEBASE METHODS

    /**
     * @brief Allocate the inner data
     *
     * used by OnDemandDataWarehouse::recvMPI to allocate an empty variable before
     * filling it with the received buffer
     *
     * @remark has to be called once even when data is requested from multiple patches/processes
     *
     * @param ai current application (has to be an implementation of PhaseField::Application)
     * @param subproblems_label label for subproblems in the DataWarehouse
     * @param material index in the DataWarehouse
     * @param patch grid patch
     */
    virtual void
    allocate (
        const ApplicationInterface * ai,
        const VarLabel * subproblems_label,
        int material,
        const Patch * patch
    ) override
    {
        auto application = dynamic_cast< const Application< Problem< VAR, STN, Field... >, false > * > ( ai );
        ASSERTMSG ( application, "PhaseField::Application required by SubProblems::allocate" );
        ASSERT ( application->getSubProblemsLabel() == subproblems_label );
        ASSERT ( !m_boundary_labels )
        ASSERT ( !m_subproblems_label )
        ASSERT ( !m_material && !m_patch && m_ghost_regions.empty() );
        ASSERT ( !m_mpi_data );
        m_boundary_labels = application->getBoundaryLabels();
        m_subproblems_label = subproblems_label;
        m_material = scinew int ( material );
        m_patch = patch;
        m_mpi_data = scinew MpiData();
    }

    /**
     * @brief Populate the ghost regions list
     *
     * used by OnDemandDataWarehouse::recvMPI to allocate an empty variable before
     * filling it with the received buffer
     *
     * @remark has to be called for each data request from foreign patches/processes
     *
     * @param low lower bound of the ghost region to partition
     * @param high higher bound of the ghost region to partition
     */
    virtual void
    addGhostRegion (
        const IntVector & low,
        const IntVector & high
    ) override
    {
        ASSERT ( m_boundary_labels );
        ASSERT ( m_subproblems_label )
        ASSERT ( m_material && m_patch );
        ASSERT ( m_mpi_data );

        for ( auto & region : m_ghost_regions )
            if ( region.getLow() <= low && high <= region.getHigh() )
            {
                // no need to add to list
                DOUTR ( g_subproblems_dbg, "not adding ghost region " << low << high << " because " << region << " already added" );
                return;
            }
            else if ( low <= region.getLow() && region.getHigh() <= high )
            {
                // replace this with other
                DOUTR ( g_subproblems_dbg, "modifying ghost region " << region << " because ghost region " << low << high << " includes it" );

                region.low() = low;
                region.high() = high;
                return;
            }

        m_ghost_regions.emplace_back ( low, high );
    }

    /**
     * @brief Create a new instance copy
     *
     * @return new SubProblems instance
     */
    virtual SubProblemsVariableBase *
    clone()
    const override
    {
        return scinew SubProblems < Problem<VAR, STN, Field...> > ( *this );
    };

    /**
     * @brief Get pointer to inner data
     *
     * @remark getRefCounted must be used by revcMPI
     *
     * @return pointer to the inner mpi data
     */
    virtual void *
    getBasePointer()
    const override
    {
        return const_cast<MpiData *> ( m_mpi_data.get_rep() );
    }

    /**
     * @brief Get BufferInfo
     *
     * Retrieve the info required to serialize inner data into an MPI message
     *
     * @param[out] buffer requested information
     */
    virtual void
    getMPIBuffer (
        BufferInfo & buffer
    ) const override
    {
        getMPIBuffer ( index_sequence_for<Field...> {}, buffer );
    }

public: // LIST

    /**
     * @brief Returns an iterator to the first element of the container
     *
     * @remark If the container is empty, the returned iterator will be equal to end()
     *
     * @return Iterator to the first element
     */
    typename std::list < Problem<VAR, STN, Field...> >::iterator
    begin() noexcept
    {
        return m_list->begin();
    }

    /**
     * @brief Returns an iterator to the first element of the container
     *
     * @remark If the container is empty, the returned iterator will be equal to end()
     *
     * @return Iterator to the first element
     */
    typename std::list < Problem<VAR, STN, Field...> >::const_iterator
    begin()
    const noexcept
    {
        return m_list->begin();
    }

    /**
     * @brief Returns an iterator to the element following the last element of the container
     *
     * @remark This element acts as a placeholder; attempting to access it results in undefined behavior
     *
     * @return Iterator to the element following the last element
     */
    typename std::list < Problem<VAR, STN, Field...> >::iterator
    end() noexcept
    {
        return m_list->end();
    }

    /**
     * @brief Returns an iterator to the element following the last element of the container
     *
     * @remark This element acts as a placeholder; attempting to access it results in undefined behavior
     *
     * @return Iterator to the element following the last element
     */
    typename std::list < Problem<VAR, STN, Field...> >::const_iterator
    end()
    const noexcept
    {
        return m_list->end();
    }

}; // struct SubProblems

template<VarType VAR, StnType STN, typename... Field> TypeDescription * SubProblems < Problem<VAR, STN, Field...> > ::td = nullptr;

template<VarType VAR, StnType STN, typename... Field> TypeDescription::Register SubProblems < Problem<VAR, STN, Field...> > ::registerMe ( getTypeDescription() );

template<VarType VAR, StnType STN, typename... Field> TypeDescription * SubProblems < Problem<VAR, STN, Field...> > ::MpiData::td = nullptr;

template<VarType VAR, StnType STN, typename... Field> TypeDescription::Register SubProblems < Problem<VAR, STN, Field...> > ::MpiData::registerMe ( getTypeDescription() );

} // namespace PhaseField

} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_PhaseField_DataTypes_SubProblems_h
