// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2026 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "teelogger.h"

#include <iostream>

#include <fcntl.h>
#include <unistd.h>


namespace EDL_NAMESPACE
{

TeeLogger::TeeLogger(const std::string& path)
{
  m_LogPath = path;

  m_OrigStdout = ::dup(STDOUT_FILENO);
  m_OrigStderr = ::dup(STDERR_FILENO);
}

TeeLogger::~TeeLogger()
{
  stopLogging();
  cleanup();
}

void TeeLogger::startLogging()
{
  if (!beginCapture()) {
    std::cerr << "Warning: unable to start edl::TeeLogger.\n";
    return;
  }
  closeLog();
  m_LogFd = ::open(m_LogPath.c_str(), O_CREAT | O_WRONLY | O_TRUNC, 0644);
}

void TeeLogger::stopLogging()
{
  endCapture();
  closeLog();
}

void TeeLogger::continueLogging()
{
  if (!beginCapture()) {
    std::cerr << "Warning: unable to start edl::TeeLogger.\n";
    return;
  }
  closeLog();
  m_LogFd = ::open(m_LogPath.c_str(), O_CREAT | O_WRONLY | O_APPEND, 0644);
}

bool TeeLogger::beginCapture()
{
  if (m_ReadFd >= 0) {
    return true;
  }
  if (m_OrigStdout < 0 || m_OrigStderr < 0) {
    return false;
  }

  int pipefd[2];
  if (::pipe(pipefd) != 0) {
    return false;
  }

  m_ReadFd = pipefd[0];
  int write_fd = pipefd[1];

  ::fflush(stdout);
  ::fflush(stderr);
  if (::dup2(write_fd, STDOUT_FILENO) < 0 || ::dup2(write_fd, STDERR_FILENO) < 0) {
    ::dup2(m_OrigStdout, STDOUT_FILENO);
    ::dup2(m_OrigStderr, STDERR_FILENO);
    ::close(m_ReadFd);
    ::close(write_fd);
    m_ReadFd = -1;
    return false;
  }
  ::close(write_fd);

  m_Worker = std::thread([this]() { pump(); });
  return true;
}

void TeeLogger::endCapture()
{
  if (m_ReadFd < 0) {
    return;
  }
  ::fflush(stdout);
  ::fflush(stderr);
  ::dup2(m_OrigStdout, STDOUT_FILENO);
  ::dup2(m_OrigStderr, STDERR_FILENO);

  if (m_Worker.joinable()) {
    m_Worker.join();
  }
  if (m_ReadFd >= 0) {
    ::close(m_ReadFd);
    m_ReadFd = -1;
  }
}

void TeeLogger::pump()
{
  char buffer[4096];
  while (true) {
    ssize_t n = ::read(m_ReadFd, buffer, sizeof(buffer));
    if (n <= 0) {
      break;
    }
    // write to terminal
    ::write(m_OrigStdout, buffer, static_cast<size_t>(n));
    // write to log file
    if (m_LogFd >= 0) {
      ::write(m_LogFd, buffer, static_cast<size_t>(n));
    }
  }
  if (m_ReadFd >= 0) {
    ::close(m_ReadFd);
    m_ReadFd = -1;
  }
}

void TeeLogger::closeLog()
{
  if (m_LogFd >= 0) {
    ::close(m_LogFd);
    m_LogFd = -1;
  }
}

void TeeLogger::cleanup()
{
  if (m_ReadFd >= 0) {
    ::close(m_ReadFd);
    m_ReadFd = -1;
  }
  if (m_OrigStdout >= 0) {
    ::close(m_OrigStdout);
    m_OrigStdout = -1;
  }
  if (m_OrigStderr >= 0) {
    ::close(m_OrigStderr);
    m_OrigStderr = -1;
  }
  closeLog();
}

} // namespace edl
