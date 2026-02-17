// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef QUICKTIMER_H
#define QUICKTIMER_H

#include <mutex>
#include <chrono>
#include <iostream>
#include <string>

namespace edl
{

class QuickTimer
{

  std::chrono::time_point<std::chrono::high_resolution_clock> m_Start;

  uint64_t    m_TotalNanoSecs = 0;
  bool        m_Running       = false;
  std::string m_Name          = "none";
  std::mutex  m_Mutex;

public:

  inline void setName(std::string name)
  {
    m_Name = name;
  }

  inline std::string name()
  {
    return m_Name;
  }

  inline void start()
  {
    std::lock_guard<std::mutex> lock(m_Mutex); // locks the mutex until it goes out of scope (unlock at end of start())
    if (!m_Running) {
      m_Start = std::chrono::high_resolution_clock::now();
      m_Running = true;
    }
  }

  inline void restart()
  {
    reset();
    start();
  }

  inline void stop()
  {
    std::lock_guard<std::mutex> lock(m_Mutex); // locks the mutex until it goes out of scope (unlock at end of stop())
    if (m_Running) {
      auto stop = std::chrono::high_resolution_clock::now();
      m_TotalNanoSecs += std::chrono::duration_cast<std::chrono::nanoseconds>(stop - m_Start).count();
      m_Running = false;
    }
  }

  inline void reset()
  {
    stop();
    m_TotalNanoSecs = 0;
  }

  inline uint64_t nanoseconds()
  {
    if (m_Running) {
      stop();
      start();
    }
    return m_TotalNanoSecs;
  }

  inline double microseconds()
  {
    return 1e-3*nanoseconds();
  }

  inline double milliseconds()
  {
    return 1e-6*nanoseconds();
  }

  inline double seconds()
  {
    return 1e-9*nanoseconds();
  }

  inline void print()
  {
    stop();
    if (m_Name != "none") {
      std::cout << "Timer \"" << m_Name << "\" : " << seconds() << " seconds" << std::endl;
    } else {
      std::cout << "Timer : " << seconds() << " seconds" << std::endl;
    }
    start();
  }

};

} // namespace EDL_NAMESPACE

#endif // QUICKTIMER_H