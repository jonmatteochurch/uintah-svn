#ifndef Uintah_Component_Arches_UPSHelper_h
#define Uintah_Component_Arches_UPSHelper_h

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>

namespace Uintah{ namespace ArchesCore {

  enum CFD_ROLE { UVELOCITY, VVELOCITY, WVELOCITY,
                  CCUVELOCITY, CCVVELOCITY, CCWVELOCITY,
                  PRESSURE, TEMPERATURE, ENTHALPY, DENSITY,
                  TOTAL_VISCOSITY };

  static inline CFD_ROLE role_string_to_enum( const std::string role ){

    if ( role == "uvelocity" ){
      return UVELOCITY;
    } else if ( role == "vvelocity" ){
      return VVELOCITY;
    } else if ( role == "wvelocity" ){
      return WVELOCITY;
    } else if ( role == "ccuvelocity" ){
      return CCUVELOCITY;
    } else if ( role == "ccvvelocity" ){
      return CCVVELOCITY;
    } else if ( role == "ccwvelocity" ){
      return CCWVELOCITY;
    } else if ( role == "pressure" ){
      return PRESSURE;
    } else if ( role == "temperature" ){
      return TEMPERATURE;
    } else if ( role == "enthalpy" ){
      return ENTHALPY;
    } else if ( role == "density" ){
      return DENSITY;
    } else if ( role == "total_viscosity" ){
      return TOTAL_VISCOSITY;
    } else {
      throw InvalidValue("Error: Cannot match role to CFD_ROLE enum. ", __FILE__, __LINE__ );
    }
  }

  static inline std::string role_enum_to_string( const CFD_ROLE role ){
    if ( role == UVELOCITY ){
      return "uvelocity";
    } else if ( role == VVELOCITY ){
      return "vvelocity";
    } else if ( role == WVELOCITY ){
      return "wvelocity";
    } else if ( role == CCUVELOCITY ){
      return "ccuvelocity";
    } else if ( role == CCVVELOCITY ){
      return "ccvvelocity";
    } else if ( role == CCWVELOCITY ){
      return "ccwvelocity";
    } else if ( role == PRESSURE ){
      return "pressure";
    } else if ( role == TEMPERATURE ){
      return "temperature";
    } else if ( role == ENTHALPY ){
      return "enthalpy";
    } else if ( role == DENSITY ){
      return "density";
    } else if ( role == TOTAL_VISCOSITY ){
      return "total_viscosity";
    } else {
      throw InvalidValue("Error: Role enum type not recognized.", __FILE__, __LINE__ );
    }
  }

  /** @brief Parse the VarID section in the UPS file for a specific CFD role **/
  static std::string parse_ups_for_role( CFD_ROLE role_enum, ProblemSpecP db, std::string mydefault="NotSet"  ){

    std::string role = role_enum_to_string( role_enum );

    ProblemSpecP db_varid = db->getRootNode()->findBlock("CFD")->findBlock("ARCHES")->findBlock("VarID");

    if ( db_varid ){
      for ( ProblemSpecP db_id = db_varid->findBlock("var"); db_id != nullptr; db_id = db_id->findNextBlock("var") ){

        std::string label="NotFound";
        std::string ups_role;

        db_id->getAttribute("label", label);
        db_id->getAttribute("role", ups_role);

        if ( ups_role == role ){
          return label;
        }
      }
    }

    return mydefault;

  }

  /** @brief Find a node with a specific attribute **/
  /** Use the *=> combination of characters to delineate between nodes. **/
  /** The * indicates the node name. **/
  /** Everything starts at the ARCHES node **/
  static ProblemSpecP inline find_node_with_att( ProblemSpecP& db, std::string start,
                                          std::string children_name,
                                          std::string att,
                                          std::string att_value ){

    ProblemSpecP db_arches = db->getRootNode()->findBlock("CFD")->findBlock("ARCHES");
    std::string delimiter = "=>";
    size_t pos = 0;
    std::vector<std::string> nodes;
    while ((pos = start.find(delimiter)) != std::string::npos) {
      std::string n = start.substr(0, pos);
      start.erase(0, pos+delimiter.length());
      nodes.push_back(n);
    }

    ProblemSpecP db_parent = db;
    for ( auto i = nodes.begin(); i != nodes.end(); i++){
      if ( db_parent->findBlock(*i)){
        db_parent = db_parent->findBlock(*i);
      } else{
        throw ProblemSetupException("Error: UPS node not found - "+*i, __FILE__, __LINE__);
      }
    }

    //once the parent is found, assume there are many children with the same
    // name. The requires that we search all potential children and
    // compare attributes to the one sought after.
    for ( ProblemSpecP db_child = db_parent->findBlock(children_name); db_child != nullptr;
          db_child = db_child->findNextBlock(children_name)){
      std::string found_att;
      db_child->getAttribute(att, found_att);
      if ( found_att == att_value ){
        return db_child;
      }
    }

    return nullptr;

  }

  // Defining default names for specific CFD variables.
  static std::string default_uVel_name{"uVelocity"};             // u-velocity, staggered
  static std::string default_vVel_name{"vVelocity"};             // v-velocity, staggered
  static std::string default_wVel_name{"wVelocity"};             // w-velocity, staggered
  static std::string default_uMom_name{"x-mom"};                 // x-momentum, staggered
  static std::string default_vMom_name{"y-mom"};                 // y-momentum, staggered
  static std::string default_wMom_name{"z-mom"};                 // z-momentum, staggered
  static std::string default_viscosity_name{"total_viscosity"};  // total viscosity (molecular + turb closure) - note: turb closure may or may not exist

}} // end Uintah::ArchesCore

#endif
