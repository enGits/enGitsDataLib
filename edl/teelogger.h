// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2026 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "edl/edl.h"

#include <string>
#include <thread>


namespace EDL_NAMESPACE
{

class TeeLogger
{

private:

  std::string m_LogPath;
  int         m_LogFd      = -1;
  int         m_ReadFd     = -1;
  int         m_OrigStdout = -1;
  int         m_OrigStderr = -1;
  std::thread m_Worker;

public:

  explicit TeeLogger(const std::string& path);

  ~TeeLogger();

  void startLogging();

  void stopLogging();

  void continueLogging();

private:

  bool beginCapture();

  void endCapture();

  void pump();

  void closeLog();
  
  void cleanup();

}; // class TeeLogger

} // namespace edl
