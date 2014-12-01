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

#include "VarContainer.hh"

#include <namespace_mouse.hh>
#include <StringTools.hh>

BEGIN_MOUSE

VarContainer::VarContainer(ConstructArg arg, bool init_list) 
  : MObject(arg)
{
  AddClassName();
  if (init_list) {
    InitList(100, 100);
    varlist = new TMDimList<real, string, 2, map_t>(this, 0.0, "");
  } else {
    varlist = NULL;
  };
};

VarContainer::VarContainer(ConstructArg arg, List *a_master) : MObject(arg)
{
  AddClassName();
  Init(a_master);
};

VarContainer::~VarContainer()
{
  delete varlist;
};

void VarContainer::Init(List *a_master) 
{
  Link(a_master);
  varlist = new TMDimList<real, string, 2, map_t>(a_master, 0.0, "");
};

void VarContainer::Reset()
{
  List *master = varlist->Master();
  delete varlist;
  varlist = new TMDimList<real, string, 2, map_t>(master, 0.0, "");
};
  
VarContainer::var1_t VarContainer::Get1DVar(string i, string j) 
{
  TMDimIndex<string> I(i);
  I = I + j;
  var1_t var;
  try { 
    var = var1_t(varlist, I); 
    return var;
  }
  catch (NotFound_error) { 
    string msg;
    StringTools::ToString(I,msg);
    msg = "unkown index '" + msg + "'";
    Error(msg); 
  };
  return var;
};

VarContainer::var2_t VarContainer::Get2DVar(string i) 
{
  TMDimIndex<string> I(i);
  var2_t var;
  try { 
    var = var2_t(varlist, I); 
    return var;
  }
  catch (NotFound_error) { 
    string msg;
    StringTools::ToString(I,msg);
    msg = "unkown index '" + msg + "'";
    Error(msg); 
  };
  return var;
};

#ifdef WITH_VAR3
VarContainer::var3_t VarContainer::Get3DVar() 
{
  TMDimIndex<string> I;
  var3_t var;
  try { 
    var = var3_t(varlist, I); 
    return var;
  }
  catch (NotFound_error) { 
    string msg;
    StringTools::ToString(I,msg);
    msg = "unkown index '" + msg + "'";
    Error(msg); 
  };
  return var;
};
#endif

bool VarContainer::FieldDefined(string field_name)
{
  bool ok = false;
  for (size_t i = 0; i < varlist->NumSubIndices(0); i++) {
    if (field_name == varlist->SubIndex(0, i)) ok = true;
  };
  return ok;
};

size_t VarContainer::VarIndex(string var_name)
{
  for (int i = 0; i < NumVars(); ++i) {
    if (VarName(i) == var_name) return i;
  };
  return NumVars();
};

END_MOUSE











