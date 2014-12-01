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

#ifndef TMDIMINDEX_H
#define TMDIMINDEX_H

namespace EDL_NAMESPACE
{
template <class TIndex> class TMDimIndex;
}

#include <iostream>
#include <string>

namespace EDL_NAMESPACE
{

using namespace std;

/**
 * a multidimensional symbolic index.
 */
template<class TIndex>
class TMDimIndex
{
protected:

  int dim;
  TIndex *I;


public:

  TMDimIndex(int a_dim);
  TMDimIndex();
  TMDimIndex(TIndex item);
  void resize(int a_dim);
  void operator=(TMDimIndex<TIndex> other);
  void operator=(TIndex item);
  bool operator==(const TMDimIndex<TIndex> other) const;
  TMDimIndex<TIndex> operator+(const TIndex item); 
  TIndex& operator[](int i) { return I[i]; }
  TIndex operator[](int i) const { return I[i]; }
  int Dim() const { return dim; }
  template<class aTIndex>
  friend ostream& operator<< (ostream& s, const TMDimIndex<aTIndex>& I);
};

template<class TIndex> ostream& operator<< (ostream& s, const TMDimIndex<TIndex>& I);


//
//.. constructor
//
template<class TIndex>
TMDimIndex<TIndex>::TMDimIndex(int a_dim)
{
  dim = a_dim;
  I = new TIndex[dim];
};

template<class TIndex>
TMDimIndex<TIndex>::TMDimIndex()
{
  dim = 0;
  I = new TIndex[dim];
};

template<class TIndex>
TMDimIndex<TIndex>::TMDimIndex(TIndex item)
{
  dim = 1;
  I = new TIndex[dim];
  I[0] = item;
};

//
//.. operator==
//
template<class TIndex>
bool TMDimIndex<TIndex>::operator==(const TMDimIndex<TIndex> other) const
{
  int i = 0;
  bool is_equal = true;
  while ((i < dim) && is_equal) {
    is_equal = (I[i] == other[i]);
    i++;
  };
  return is_equal;
};

//
//.. Resize
//
template<class TIndex>
void TMDimIndex<TIndex>::resize(int a_dim)
{
  delete [] I;
  dim = a_dim;
  I = new TIndex[dim];
};

//
//.. operator=
//
template<class TIndex>
void TMDimIndex<TIndex>::operator=(TMDimIndex<TIndex> other)
{
  resize(other.dim);
  for (int i = 0; i < dim; i++) I[i] = other[i];
};

//
//.. operator=
//
template<class TIndex>
void TMDimIndex<TIndex>::operator=(TIndex item)
{
  resize(1);
  I[0] = item;
};

//
//.. operator<<
//
template<class TIndex>
ostream& operator<<(ostream &s, const TMDimIndex<TIndex> &I)
{
  string sep = "";
  s << "(";
  for (int i = 0; i < I.Dim(); i++) {
    if (i == 1) sep = ",";
    s << sep << I[i];
  };
  s << ")";
  return s;
};

//
//.. operator+
//
template<class TIndex>
TMDimIndex<TIndex> TMDimIndex<TIndex>::operator+(const TIndex item)
{
  TMDimIndex<TIndex> new_I(dim+1);
  for (int i = 0; i < dim; i++) new_I[i] = I[i];
  new_I[dim] = item;
  return new_I;
};


} // namespace

#endif










