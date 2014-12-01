// !!
// This is a part of MOUSE, a library for PDE's on unstructured grids
// Copyright (C) 1999 Oliver Gloth <oliver@vug.uni-duisburg.de>
// Institut fuer Verbrennung und Gasdynamik (Universitaet Duisburg)
// institute for combustion and gas dynamics (university of Duisburg)
// Thursday, 1 April, 1999 Duisburg, Germany
//
// please see http://www.vug.uni-duisburg.de/MOUSE for more information
// please send any questions or suggestions to mouse@www.vug.uni-duisburg.de
//  
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software 
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// !!

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
class VarContainer : public MObject, public List
{
public:
  typedef BMap<real, string, 2> map_t;

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
  typedef TSDimVar<real, string, 2, map_t> var1_t;

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
  typedef TMDimVar<real, string, 2, map_t, var1_t> var2_t;

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
  typedef TMDimVar<real, string, 2, map_t, var2_t> var3_t;
#endif
  
protected:
  TMDimList<real, string, 2, map_t> *varlist;
  
public:
  VarContainer(ConstructArg arg, bool init_list = true);
  VarContainer(ConstructArg arg, List* a_master);
  virtual ~VarContainer();
  void Init(List* a_master);

  /** 
     Add a variable entry. (2.index) to the container. 
     @param var_name the name of the new variable
  */
  void AddVar(string var_name) { 
    if (!varlist->IndexExists(1, var_name)) varlist->AddIndex(1, var_name); 
  };

  /**
     Add a field entry. (1.index) to the container. 
     @param field_name the name of the new field
  */
  void AddField(string field_name) {
    if (!varlist->IndexExists(0, field_name)) varlist->AddIndex(0, field_name);
  };

  /** 
     Get the number of variables 
     @return the number of variables 
  */
  size_t NumVars() { return varlist->NumSubIndices(1); };

  /**
     get the name of a variable 
     @param i the number of the variable
     @return the name of the variable
  */ 
  string VarName(int i) { return varlist->SubIndex(1, i); };

  /**
   * get the index of a variable.
   * @param name the name of the variable
   * @return the index of the variable
   */
  size_t VarIndex(string var_name);

  /** 
     Get the number of fields
     @return the number of fields
  */
  size_t NumFields() { return varlist->NumSubIndices(0); };

  /**
   * get the name of a field.
   * @param i the number of the field
   * @return the name of the field
   */
  string FieldName(size_t i) { return varlist->SubIndex(0, i); };

  /**
     Find out if a field is defined or not.
     @param field_name the name of the field to check
     @return true if this field exists
   */ 
  bool FieldDefined(string field_name);

  /**
     add new empty entries.
     @param n the number of new entries
  */
  void ExtendBy(size_t N) {
    for (size_t i = 0; i < N; i++) varlist->AddEntry();
  };

  /**
   * Reset the VarContainer (delete all variables).
   */
  void Reset();

  /**  Get a one dimensional variable.\\
       (see also the documentation for var1_t)
       @param i the field name of the variable
       @param j the variable name
       @return the one dimensional variable */ 
  var1_t Get1DVar(string i, string j);

  /**  Get a two dimensional variable.\\
       (see also the documentation for var2_t)
       @param i the field name of the variable
       @return the two dimensional variable */ 
  var2_t Get2DVar(string i);

#ifdef WITH_VAR3
  /**  Get a three dimensional variable.\\
       (see also the documentation for var3_t)
       @return the three dimensional variable */ 
  var3_t Get3DVar();
#endif
};

} // namespace

#endif
























