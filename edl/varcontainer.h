// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef VARCONTAINER_H
#define VARCONTAINER_H

#include "edl.h"

namespace EDL_NAMESPACE
{
class VarContainer;
}

#include "edl/list.h"
#include "edl/tmdimlist.h"
#include "edl/bmap.h"
#include "edl/tsdimvar.h"
#include "edl/tmdimvar.h"


namespace EDL_NAMESPACE
{

/**
 * This is a combination of a TMDimList with two symbolic dimensions
 * and an MObject. By default it is initialized with zero entries
 * and an increment of 100. VarContainer represents a three dimensional
 * list, where the first two inidices are symbolic (string). It is
 * possible to create referencing variables, which will be automatically
 * updated if the VarContainer expands or shrinks. This functionality is
 * implemented in the more general class TMDimList.
 */
class VarContainer : public List
{
public:

  typedef BMap<real, std::string, 2> map_t;

  /**
   * This is an automatically updated one dimensional variable, where the index
   * represents the point address.
   * The method Get1DVar is used to initialize such a variable.
   * The call
   * <pre>
   * VarContainer::var1_t res_rho = VC->Get1DVar("res", "rho");
   * </pre>
   * for example, will initialize the variable res_rho as an one dimensional
   * reference to the residual field (named "res") of the density variable
   * (named "rho"). VC in this case is of the type VarContainer*.
   * This variable can now ne used in the same manner a C-style
   * pointer would be used. The programmer does not have to take care about
   * updating it if the underlaying list expands or shrinks. 
   */
  typedef TSDimVar<real, std::string, 2, map_t> var1_t;

  /**
   * This is similiar to var1_t, but in this case representing
   * a two dimensional variable.
   * The call
   * <pre>
   * VarContainer::var2_t res = VC->Get2DVar("res");
   * </pre>
   * will initialize the variable res as a two dimensional
   * reference to the residual field (named "res"). Please note that if you use the
   * variable like a pointer, both indices are of the type size_t. That means once you
   * get the variable with the Get2DVar call, like described above, you use it like a C-style
   * array. If you have declared the following
   * <pre>
   * VarContainer *VC;
   * ...
   * VarContainer::var2_t res = VC->Get2DVar("res");
   * </pre>
   * you can write the following loop to set all entries of the "res" field to zero:
   * <pre>
   * FORALL(j, VC->) {
   *   for(size_t j = 0; j < VC->NumVars(); j++) {
   *     res[i][j] = 0;
   *   }
   * }
   * </pre>
   * (see also var1_t) 
   */
  typedef TMDimVar<real, std::string, 2, map_t, var1_t> var2_t;

#ifdef WITH_VAR3
  /** 
   * This is similiar to var1_t, but in this case representing
   * a three dimensional variable.
   * The call
   * <pre>
   * VarContainer::var3_t var = VC->Get3DVar();
   * </pre>
   * will initialize the variable var as a three dimensional
   * reference to all fields and variables of the VarContainer *VC.<br>
   * (see also var1_t and var2_t) 
   */
  typedef TMDimVar<real, std::string, 2, map_t, var2_t> var3_t;
#endif
  

protected:

  TMDimList<real, std::string, 2, map_t> *m_VarList;
  

public:

  VarContainer(bool init_list = true);
  VarContainer(List* a_master, std::string link_name = "__none");
  virtual ~VarContainer();
  void init(List* a_master, std::string link_name = "__none");

  /** 
     Add a variable entry. (2.index) to the container. 
     @param var_name the name of the new variable
  */
  void addVar(std::string var_name)
  {
    if (!m_VarList->indexExists(1, var_name)) m_VarList->addIndex(1, var_name);
  }

  /**
     Add a field entry. (1.index) to the container. 
     @param field_name the name of the new field
  */
  void addField(std::string field_name)
  {
    if (!m_VarList->indexExists(0, field_name)) m_VarList->addIndex(0, field_name);
  }

  /** 
     Get the number of variables 
     @return the number of variables 
  */
  size_t numVars() { return m_VarList->numSubIndices(1); }

  /**
     get the name of a variable 
     @param i the number of the variable
     @return the name of the variable
  */ 
  std::string varName(int i) { return m_VarList->subIndex(1, i); }

  /**
   * get the index of a variable.
   * @param name the name of the variable
   * @return the index of the variable
   */
  size_t varIndex(std::string var_name);

  /** 
     Get the number of fields
     @return the number of fields
  */
  size_t numFields() { return m_VarList->numSubIndices(0); }

  /**
   * get the name of a field.
   * @param i the number of the field
   * @return the name of the field
   */
  std::string fieldName(size_t i) { return m_VarList->subIndex(0, i); }

  /**
     Find out if a field is defined or not.
     @param field_name the name of the field to check
     @return true if this field exists
   */ 
  bool fieldDefined(std::string field_name);

  /**
     add new empty entries.
     @param n the number of new entries
  */
  void extendBy(size_t N)
  {
    for (size_t i = 0; i < N; i++) m_VarList->addEntry();
  }

  /**
   * Reset the VarContainer (delete all variables).
   */
  void reset();

  /**  Get a one dimensional variable.\\
       (see also the documentation for var1_t)
       @param i the field name of the variable
       @param j the variable name
       @return the one dimensional variable */ 
  var1_t get1DVar(std::string i, std::string j);

  /**  Get a two dimensional variable.\\
       (see also the documentation for var2_t)
       @param i the field name of the variable
       @return the two dimensional variable */ 
  var2_t get2DVar(std::string i);

#ifdef WITH_VAR3
  /**  Get a three dimensional variable.\\
       (see also the documentation for var3_t)
       @return the three dimensional variable */ 
  var3_t get3DVar();
#endif

};

} // namespace

#endif
























