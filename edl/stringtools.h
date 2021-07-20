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
#ifndef StringTools_HH
#define StringTools_HH

#include "edl/edl.h"

#include <string>
#include <vector>
#include <list>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace EDL_NAMESPACE
{

namespace StringTools
{
  /**
   * convert a std::string into another type.
   * <b>Attention</b> this method only works for a type <i>T</i> for
   * which an input operator <i>operator>></i>exists.
   * @param s the std::string to convert
   * @param t a variable which will contain the converted value
   * @return true if the operation succeeded
   */
  template <class T> bool stringTo(std::string s, T &t);

  /**
   * convert another type into a string.
   * <b>Attention</b> this method only works for a type <i>T</i> for
   * which an output operator <i>operator<<</i>exists.
   * @param t the variable to convert
   * @param s the resulting string
   */
  template <class T> void toString(T t, std::string &s);

  /**
   * convert another type into a string.
   * <b>Attention</b> this method only works for a type <i>T</i> for
   * which an output operator <i>operator<<</i>exists.
   * @param t the variable to convert
   * @return the resulting string
   */
  template <class T> std::string toString(T t, int fill = 0);

  /**
   * replace a character with a different one.
   * @param s the source string
   * @param c_orig the character to replace
   * @param c_new the new character
   * return a std::string with the replaced characters
   */
  std::string replace(std::string s, char c_orig, char c_new);

  /**
   * append a character to the left side of a std::string until a given length has been reached.
   * @param s the source string
   * @param c the character to append
   * @param l the desired length
   * @return the filled string
   */
  std::string leftFill(std::string s, char c, size_t l);

  /**
   * append a character to the right side of a std::string until a given length has been reached.
   * @param s the source string
   * @param c the character to append
   * @param l the desired length
   * @return the filled string
   */
  std::string rightFill(std::string s, char c, size_t l);

  /**
   * read an entire line from an input stream
   * @param s the stream to read from
   * @return the line read
   */
  std::string readLine(std::istream &s);

  /**
   * extract a part of a string.
   * @param s the original string
   * @param i1 index of the first character of the intended sub-string
   * @param i2 index of the last character of the intended sub-string
   * @return the sub-string
   */
  std::string subString(std::string s, size_t i1, size_t i2);

  /**
   * extract a part at the end of a string.
   * @param s the original string
   * @param n number of characters of the intended sub-string
   * @return the sub-string
   */
  std::string right(std::string s, size_t n);

  /**
   * extract a part at the beginning of a string.
   * @param s the original string
   * @param n number of characters of the intended sub-string
   * @return the sub-string
   */
  std::string left(std::string s, size_t n);

  /**
   * split a std::string into several strings.
   * @param s the std::string to be split
   * @param sub the list which the sub strings will be appeneded to
   * @param c the dividing character
   * @return the number of sub strings
   */
  int split(std::string s, std::list<std::string> &sub, char c = ' ');

  /**
   * split a std::string into several strings.
   * @param s the std::string to be split
   * @param sub the vector which will contain the sub strings
   * @param c the dividing character
   * @return the number of sub strings
   */
  int split(std::string s, std::vector<std::string> &sub, char c = ' ');

  template <class T>
  inline bool stringTo(std::string s, T &t)
  {
    s += " -";
    std::istringstream stream(s);
    stream >> t;
    return stream.good();
  }

  template <>
  inline bool stringTo(std::string s, std::string &t)
  {
    t = s;
    return true;
  }

  template <class T>
  void toString(T t, std::string &s)
  {
    std::ostringstream stream;
    stream << t;
    s = stream.str();
  }
  
  template <class T>
  std::string toString(T t, int fill)
  {
    std::ostringstream stream;
    stream << t;
    std::string s = stream.str();
    if      (fill < 0) s = leftFill(s,' ',-fill);
    else if (fill > 0) s = rightFill(s,' ',fill);
    return s;
  }
  
  inline std::string replace(std::string s, char c_orig, char c_new)
  {
    std::string s_new = "";
    std::string::iterator i = s.begin();
    while (i != s.end()) {
      if (*i == c_orig) s_new += c_new;
      else s_new += *i;
      i++;
    };
    return s_new;
  }

  inline std::string remove(std::string s, char c)
  {
    std::string s_new = "";
    std::string::iterator i = s.begin();
    while (i != s.end()) {
      if (*i != c) s_new += *i;
      i++;
    };
    return s_new;
  }

  inline std::string readLine(std::istream &s)
  {
    std::string line = "";
    bool done = false;
    do {
      char c;
      s.get(c);
      if (c == '\n') {
        done = true;
      } else {
        line += c;
        done = (s.eof() || s.fail());
      };
    } while (!done);
    return line;
  }
  
  inline std::string subString(std::string s, size_t i1, size_t i2)
  {
    std::string sub = "";
    size_t i = i1;
    while ((i <= i2) && (i < s.size())) {
      sub += s[i];
      ++i;
    };
    return sub;
  }

  inline std::string leftFill(std::string s, char c, size_t l)
  { 
    while (s.size() < l) {
      s = c + s;
    };
    return s; 
  };

  inline std::string rightFill(std::string s, char c, size_t l)
  {
    while (s.size() < l) {
      s = s + c;
    };
    return s; 
  }

  inline std::string right(std::string s, size_t n)
  {
    return subString(s,s.size()-n,s.size()-1);
  }

  inline std::string left(std::string s, size_t n)
  {
    return subString(s,0,n-1);
  }

  inline std::string trim(std::string s)
  {
    std::string result = "";
    for (auto c : s) {
      if (!isspace(c) || result.size()) {
        result += c;
      }
    }
    while (result.size()) {
      if (isspace(result.back())) {
        result.pop_back();
      } else {
        break;
      }
    }
    return result;
  }

  inline std::vector<std::string> split(std::string s, std::string delimiter, bool trim_result=true)
  {
    size_t pos_start = 0;
    size_t pos_end;
    size_t delim_len = delimiter.length();
    //
    std::string              token;
    std::vector<std::string> result;
    //
    while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
      token = s.substr (pos_start, pos_end - pos_start);
      pos_start = pos_end + delim_len;
      if (trim_result) {
        result.push_back(trim(token));
      } else {
        result.push_back(token);
      }
    }
    //
    if (trim_result) {
      result.push_back(trim(s.substr(pos_start)));
    } else {
      result.push_back(s.substr(pos_start));
    }
    return result;
  }

  inline int countChar(std::string s, char count_char)
  {
    int count = 0;
    for (char c : s) {
      if (c == count_char) {
        ++count;
      }
    }
    return count;
  }

  inline std::vector<std::string> quotedSplit(std::string s, std::string delimiter, char quote_char='"', bool trim_result=true)
  {
    auto quote_str = std::string("") + quote_char;
    auto parts = split(s, delimiter, trim_result);
    std::vector<std::string> clean_parts;
    clean_parts.reserve(parts.size());
    int i = 0;
    while (i < parts.size()) {
      if (parts[i].empty()) {
        clean_parts.push_back(parts[i]);
        ++i;
      } else if (countChar(parts[i], quote_char) % 2 == 0) {
        clean_parts.push_back(parts[i]);
        ++i;
      } else {
        clean_parts.push_back(StringTools::remove(parts[i], '"'));
        while (i < parts.size() - 1) {
          ++i;
          clean_parts.back() += "," + StringTools::remove(parts[i], '"');
          if (StringTools::right(parts[i],1) == quote_str) {
            ++i;
            break;
          }
        }
      }
    }
    return clean_parts;
  }

  inline std::string toLower(std::string s)
  {
    for (auto& c :s) {
      c = std::tolower(c);
    }
    return s;
  }

  inline std::string toUpper(std::string s)
  {
    for (auto& c :s) {
      c = std::toupper(c);
    }
    return s;
  }

  inline std::string readQuotedLine(std::istream& s, char quote_char='"')
  {
    bool in_quotes = false;
    std::string line = "";
    while (!s.eof()) {
      char c = s.get();
      if (c == -1) {
        break;
      }
      if (c == quote_char) {
        in_quotes = !in_quotes;
      }
      if (c == '\n' && !in_quotes) {
        return line;
      } else {
        line += c;
      }
    }
    return line;
  }

} // namespace StringTools

} // namespace edl

#endif


