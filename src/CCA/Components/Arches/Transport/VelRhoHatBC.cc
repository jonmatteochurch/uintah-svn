#include <CCA/Components/Arches/Transport/VelRhoHatBC.h>
#include <CCA/Components/Arches/UPSHelper.h>

using namespace Uintah;
typedef ArchesFieldContainer AFC;

//--------------------------------------------------------------------------------------------------
VelRhoHatBC::VelRhoHatBC( std::string task_name, int matl_index ) :
AtomicTaskInterface( task_name, matl_index )
{
}

//--------------------------------------------------------------------------------------------------
VelRhoHatBC::~VelRhoHatBC()
{
}

//--------------------------------------------------------------------------------------------------
void VelRhoHatBC::problemSetup( ProblemSpecP& db ){

  using namespace ArchesCore;

  m_xmom = default_uMom_name;
  m_ymom = default_vMom_name;
  m_zmom = default_wMom_name;
  m_uVel = parse_ups_for_role( UVELOCITY_ROLE, db, default_uVel_name );
  m_vVel = parse_ups_for_role( VVELOCITY_ROLE, db, default_vVel_name );
  m_wVel = parse_ups_for_role( WVELOCITY_ROLE, db, default_wVel_name );
}

//--------------------------------------------------------------------------------------------------
void VelRhoHatBC::create_local_labels(){
}

//--------------------------------------------------------------------------------------------------
void VelRhoHatBC::register_timestep_eval( std::vector<AFC::VariableInformation>& variable_registry,
                                 const int time_substep, const bool pack_tasks ){
  register_variable( m_xmom, AFC::MODIFIES, variable_registry, m_task_name );
  register_variable( m_ymom, AFC::MODIFIES, variable_registry, m_task_name );
  register_variable( m_zmom, AFC::MODIFIES, variable_registry, m_task_name );
  register_variable( m_uVel, AFC::REQUIRES, 0, AFC::OLDDW, variable_registry, time_substep, m_task_name );
  register_variable( m_vVel, AFC::REQUIRES, 0, AFC::OLDDW, variable_registry, time_substep, m_task_name );
  register_variable( m_wVel, AFC::REQUIRES, 0, AFC::OLDDW, variable_registry, time_substep, m_task_name );



}

void VelRhoHatBC::eval( const Patch* patch, ArchesTaskInfoManager* tsk_info ){

  SFCXVariable<double>& xmom = tsk_info->get_field<SFCXVariable<double> >( m_xmom );
  SFCYVariable<double>& ymom = tsk_info->get_field<SFCYVariable<double> >( m_ymom );
  SFCZVariable<double>& zmom = tsk_info->get_field<SFCZVariable<double> >( m_zmom );


  const BndMapT& bc_info = m_bcHelper->get_boundary_information();
  const double possmall = 1e-16;

  for ( auto i_bc = bc_info.begin(); i_bc != bc_info.end(); i_bc++ ){

    const bool on_this_patch = i_bc->second.has_patch(patch->getID());
    if ( !on_this_patch ) continue;

    Uintah::ListOfCellsIterator& cell_iter = m_bcHelper->get_uintah_extra_bnd_mask( i_bc->second, patch->getID() );
    IntVector iDir = patch->faceDirection( i_bc->second.face );
    Patch::FaceType face = i_bc->second.face;
    BndTypeEnum my_type = i_bc->second.type;
    int sign = iDir[0] + iDir[1] + iDir[2];
    int bc_sign = 0;
    int move_to_face_value = ( sign < 1 ) ? 1 : 0;

    IntVector move_to_face(std::abs(iDir[0])*move_to_face_value,
                           std::abs(iDir[1])*move_to_face_value,
                           std::abs(iDir[2])*move_to_face_value);

    Array3<double>* var = NULL;

    if ( face == Patch::xminus || face == Patch::xplus ){
      var     = &xmom;
      constSFCXVariable<double>& old_var = tsk_info->get_field<constSFCXVariable<double> >( m_uVel );

      if ( my_type == OUTLET_BC ){
        bc_sign = 1.;
      } else if ( my_type == PRESSURE_BC ){
        bc_sign = -1.;
      }

      sign = bc_sign * sign;

      if ( my_type == OUTLET_BC || my_type == PRESSURE_BC ){
        // This applies the mostly in (pressure)/mostly out (outlet) boundary condition
        parallel_for(cell_iter.get_ref_to_iterator(),cell_iter.size(), [&] (const int i,const int j,const int k) {
          int i_f = i + move_to_face[0]; // cell on the face
          int j_f = j + move_to_face[1];
          int k_f = k + move_to_face[2];

          int im = i_f - iDir[0];// first interior cell
          int jm = j_f - iDir[1];
          int km = k_f - iDir[2];

          int ipp = i_f + iDir[0];// extra cell face in the last index (mostly outwardly position)
          int jpp = j_f + iDir[1];
          int kpp = k_f + iDir[2];

          if ( sign * old_var(i_f,j_f,k_f) > possmall ){
            // du/dx = 0
            (*var)(i_f,j_f,k_f)= (*var)(im,jm,km);
          } else {
            // shut off the hatted value to encourage the mostly-* condition
            (*var)(i_f,j_f,k_f) = 0.0;
          }

          (*var)(ipp,jpp,kpp) = (*var)(i_f,j_f,k_f);

        });
      }


    }  else if ( face == Patch::yminus || face == Patch::yplus ){
      var = &ymom;
      constSFCYVariable<double>& old_var = tsk_info->get_field<constSFCYVariable<double> >( m_vVel );


      if ( my_type == OUTLET_BC ){
        bc_sign = 1.;
      } else if ( my_type == PRESSURE_BC ){
        bc_sign = -1.;
      }

      sign = bc_sign * sign;

      if ( my_type == OUTLET_BC || my_type == PRESSURE_BC ){
        // This applies the mostly in (pressure)/mostly out (outlet) boundary condition
        parallel_for(cell_iter.get_ref_to_iterator(),cell_iter.size(), [&] (const int i,const int j,const int k) {
          int i_f = i + move_to_face[0]; // cell on the face
          int j_f = j + move_to_face[1];
          int k_f = k + move_to_face[2];

          int im = i_f - iDir[0];// first interior cell
          int jm = j_f - iDir[1];
          int km = k_f - iDir[2];

          int ipp = i_f + iDir[0];// extra cell face in the last index (mostly outwardly position)
          int jpp = j_f + iDir[1];
          int kpp = k_f + iDir[2];

          if ( sign * old_var(i_f,j_f,k_f) > possmall ){
            // du/dx = 0
            (*var)(i_f,j_f,k_f)= (*var)(im,jm,km);
          } else {
            // shut off the hatted value to encourage the mostly-* condition
            (*var)(i_f,j_f,k_f) = 0.0;
          }

          (*var)(ipp,jpp,kpp) = (*var)(i_f,j_f,k_f);

        });
      }

    } else {

      var = &zmom;
      constSFCZVariable<double>& old_var = tsk_info->get_field<constSFCZVariable<double> >( m_wVel );

      if ( my_type == OUTLET_BC ){
        bc_sign = 1.;
      } else if ( my_type == PRESSURE_BC ){
        bc_sign = -1.;
      }

      sign = bc_sign * sign;

      if ( my_type == OUTLET_BC || my_type == PRESSURE_BC ){
        // This applies the mostly in (pressure)/mostly out (outlet) boundary condition
        parallel_for(cell_iter.get_ref_to_iterator(),cell_iter.size(), [&] (const int i,const int j,const int k) {
          int i_f = i + move_to_face[0]; // cell on the face
          int j_f = j + move_to_face[1];
          int k_f = k + move_to_face[2];

          int im = i_f - iDir[0];// first interior cell
          int jm = j_f - iDir[1];
          int km = k_f - iDir[2];

          int ipp = i_f + iDir[0];// extra cell face in the last index (mostly outwardly position)
          int jpp = j_f + iDir[1];
          int kpp = k_f + iDir[2];

          if ( sign * old_var(i_f,j_f,k_f) > possmall ){
            // du/dx = 0
            (*var)(i_f,j_f,k_f)= (*var)(im,jm,km);
          } else {
            // shut off the hatted value to encourage the mostly-* condition
            (*var)(i_f,j_f,k_f) = 0.0;
          }

          (*var)(ipp,jpp,kpp) = (*var)(i_f,j_f,k_f);

        });
      }


    }

  }
}
