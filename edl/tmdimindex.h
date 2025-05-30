// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  friend std::ostream& operator<< (std::ostream& s, const TMDimIndex<aTIndex>& I);
};

template<class TIndex> std::ostream& operator<< (std::ostream& s, const TMDimIndex<TIndex>& I);


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
std::ostream& operator<<(std::ostream &s, const TMDimIndex<TIndex> &I)
{
  std::string sep = "";
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










