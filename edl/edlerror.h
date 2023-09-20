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
#ifndef EDLERROR_H
#define EDLERROR_H

#include "edl.h"

namespace EDL_NAMESPACE {
class EdlError;
}

#include <cstddef>
#include <iostream>

#include <string>
#include <QString>

namespace EDL_NAMESPACE
{

/** base class for exceptions.
 */
class EdlError
{
  /** the name of the exception. 
   */
  std::string m_Name;

  /** a message to print if the exception is not caught.
   */
  std::string m_Message;

protected:

  /** has to be called by every derived constructor.
   *  (similiar to AddClassName() for MObject)
   */
  virtual void setName();

public:

  /** constructor.
   */
  EdlError(std::string message);

  /** destructor.
   */
  virtual ~EdlError();

  /** get the name of this exception.
   *  @return the name of the exception
   */
  std::string name();

  /** get the message of the exception.
   *  @return the message
   */
  std::string message() { return m_Message; }

};

#define EDL_BUG \
{ \
  QString line; \
  line.setNum(__LINE__); \
  QString msg = "A problem occurred within enGitsDataLib"; \
  msg += QString("\n\nfile: ") + __FILE__ + "\nline:" + line + "\n\n"; \
  std::cerr << qPrintable(msg) << std::endl; \
  throw EdlError(qPrintable(msg)); \
};

} // namespace


#endif





