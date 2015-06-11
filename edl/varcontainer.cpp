// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2015 enGits GmbH                                         +
// +                                                                    +
// + enGitsDataLib is free software: you can redistribute it and/or     +
// + modify it under the terms of the GNU Lesser General Public License +
// + as published by the Free Software Foundation, either version 3 of  +
// + the License, or (at your option) any later version.                +
// +                                                                    +
// + enGitsDataLib is distributed in the hope that it will be useful,   +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of     +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      +
// + GNU Lesser General Public License for more details.                +
// +                                                                    +
// + You should have received a copy of the GNU Lesser General Public   +
// + License along with enGitsDataLib.                                  +
// + If not, see <http://www.gnu.org/licenses/>.                        +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "edl/varcontainer.h"
#include "edl/stringtools.h"

namespace EDL_NAMESPACE
{

VarContainer::VarContainer(bool init_list)
{
  if (init_list) {
    initList(100, 100);
    m_VarList = new TMDimList<real, string, 2, map_t>(this, 0.0, "");
  } else {
    m_VarList = NULL;
  }
}

VarContainer::VarContainer(List *a_master)
{
  init(a_master);
}

VarContainer::~VarContainer()
{
  delete m_VarList;
}

void VarContainer::init(List *a_master)
{
  link(a_master);
  m_VarList = new TMDimList<real, string, 2, map_t>(a_master, 0.0, "");
}

void VarContainer::reset()
{
  List *master = m_VarList->master();
  delete m_VarList;
  m_VarList = new TMDimList<real, string, 2, map_t>(master, 0.0, "");
}
  
VarContainer::var1_t VarContainer::get1DVar(string i, string j)
{
  TMDimIndex<string> I(i);
  I = I + j;
  var1_t var;
  try { 
    var = var1_t(m_VarList, I);
    return var;
  }
  catch (NotFound_error) { 
    string msg;
    StringTools::toString(I,msg);
    msg = "unkown index '" + msg + "'";
    throw EdlError(msg);
  }
  return var;
}

VarContainer::var2_t VarContainer::get2DVar(string i)
{
  TMDimIndex<string> I(i);
  var2_t var;
  try { 
    var = var2_t(m_VarList, I);
    return var;
  }
  catch (NotFound_error) { 
    string msg;
    StringTools::toString(I,msg);
    msg = "unkown index '" + msg + "'";
    throw EdlError(msg);
  }
  return var;
}

#ifdef WITH_VAR3
VarContainer::var3_t VarContainer::get3DVar()
{
  TMDimIndex<string> I;
  var3_t var;
  try { 
    var = var3_t(m_VarList, I);
    return var;
  }
  catch (NotFound_error) { 
    string msg;
    StringTools::toString(I,msg);
    msg = "unkown index '" + msg + "'";
    throw EdlError(msg);
  }
  return var;
}
#endif

bool VarContainer::fieldDefined(string field_name)
{
  bool ok = false;
  for (size_t i = 0; i < m_VarList->numSubIndices(0); i++) {
    if (field_name == m_VarList->subIndex(0, i)) ok = true;
  };
  return ok;
}

size_t VarContainer::varIndex(string var_name)
{
  for (int i = 0; i < numVars(); ++i) {
    if (varName(i) == var_name) return i;
  }
  return numVars();
}

} // namespace











