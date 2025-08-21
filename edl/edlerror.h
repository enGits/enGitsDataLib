// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef EDLERROR_H
#define EDLERROR_H

#include <iostream>
#include <string>
#include <sstream>

#include "edl.h"


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
  std::stringstream msg; \
  msg << "A problem occurred within enGitsDataLib"; \
  msg << "\nfile: " << __FILE__ << "\nline: " << __LINE__ << "\n"; \
  std::cerr << msg.str() << std::endl; \
  throw EdlError(msg.str()); \
};

} // namespace


#endif





