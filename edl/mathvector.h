// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef mathvector_h
#define mathvector_h

#include "edl.h"

namespace EDL_NAMESPACE
{
template <class V> struct MathVector;
template <class T, unsigned int DIM> class StaticVector;
}

#include <iostream>
#include <vector>
#include <cmath>
#include <typeinfo>

namespace EDL_NAMESPACE
{

typedef unsigned int uint_t;

#include "mathvector_structs.h"


// wrapper class for arrays to provide mathematical functionality
// =======================================================================================

template <class V>
struct MathVector : public V
{

  // types
  // -----
  typedef typename V::value_type scalar_t;
  typedef typename V::value_type value_type;


  // constructors
  // ------------

  MathVector() : V() {}
  MathVector(const MathVector<V>& vec);
  MathVector(const value_type *v);
  MathVector(const value_type v);
  MathVector(const value_type v1, const value_type v2, const value_type v3, const value_type v4);
  MathVector(const value_type v1, const value_type v2, const value_type v3);
  MathVector(const value_type v1, const value_type v2);
  
  template <class L, class O, class R>
  MathVector(const ParseNode<L,O,R> &expr);


  // operators
  // ---------
  void operator=  (const MathVector<V> &vec);
  void operator-= (const MathVector<V> &vec); 
  void operator+= (const MathVector<V> &vec);
  void operator*= (const scalar_t s);

  template <class C> void operator=  (const C &vec);

  // assignment to an expression
  template <class L, class O, class R> void operator= (const ParseNode<L,O,R> &expr);


  // other things
  // ------------
  MathVector<V> cross(const MathVector<V> &vec) const;
  scalar_t abs() const;
  scalar_t abs2() const;
  MathVector<V> normalise();
  MathVector<V> normalised();
  scalar_t* c_array() const;
  uint_t dim() { return this->size(); }
  void minimisePerCoord(const MathVector<V>& vec);
  void maximisePerCoord(const MathVector<V>& vec);
  

  // STL
  // ---
  class iterator;
  class const_iterator;

  class iterator
  {
    friend class const_iterator;
    MathVector<V> *vec;
    uint_t i;

  public:
    typedef std::random_access_iterator_tag iterator_category;
    typedef typename V::value_type value_type;
    typedef uint_t difference_type;
    typedef value_type* pointer;
    typedef value_type& reference;

    iterator(MathVector<V> *a_vec = NULL, uint_t an_i = 0) : vec(a_vec), i(an_i) {}
    bool operator==(const iterator &iter) const { return (iter.vec == vec) && (iter.i == i); }
    bool operator==(const const_iterator &iter) const;
    bool operator!=(const const_iterator &iter) const { return !operator==(iter); }
    iterator operator++() { ++i; return iterator(vec,i); }
    iterator operator++(int) { uint_t j = i; ++i; return iterator(vec,j); }
    iterator operator--() { --i; return iterator(vec,i); }
    iterator operator--(int) { uint_t j = i; --i; return iterator(vec,j); }
    void operator+=(uint_t n) { i += n; }
    void operator-=(uint_t n) { i -= n; }
    scalar_t& operator*() const { return (*vec)[i]; }
    scalar_t& operator[](uint_t n) const { return (*vec)[i+n]; }
    bool operator<(const iterator &iter) const { return i < iter.i; }
    uint_t operator-(const iterator &iter) const { return i - iter.i; }
  };

  class const_iterator
  {
    MathVector<V> *vec;
    uint_t i;

  public:
    typedef std::random_access_iterator_tag iterator_category;
    typedef typename V::value_type value_type;
    typedef uint_t difference_type;
    typedef value_type* pointer;
    typedef value_type& reference;

    const_iterator(MathVector<V> *a_vec = NULL, uint_t an_i = 0) : vec(a_vec), i(an_i) {}
    const_iterator(const iterator &iter) : vec(iter.vec), i(iter.i) {}
    void operator=(const iterator &iter) { vec = iter.vec; i = iter.i; }
    bool operator==(const const_iterator &iter) const { return (iter.vec == vec) && (iter.i == i); }
    bool operator==(const iterator &iter) const;
    bool operator!=(const iterator &iter) const { return !operator==(iter); }
    const_iterator operator++() { ++i; return iterator(vec,i); }
    const_iterator operator++(int) { uint_t j = i; ++i; return iterator(vec,j); }
    const_iterator operator--() { --i; return iterator(vec,i); }
    const_iterator operator--(int) { uint_t j = i; --i; return iterator(vec,j); }
    void operator+=(uint_t n) { i += n; }
    void operator-=(uint_t n) { i -= n; }
    scalar_t operator*() const { return (*vec)[i]; }
    scalar_t operator[](uint_t n) const { return (*vec)[i+n]; }
    bool operator<(const iterator &iter) const { return i < iter.i; }
    uint_t operator-(const iterator &iter) const { return i - iter.i; }
  };

  iterator       begin()       { return iterator(this); }
  iterator       end()         { return iterator(this,this->size()); }
  const_iterator begin() const { return iterator(this); }
  const_iterator end()   const { return iterator(this,this->size()); }
};




// static array class for small vectors
// ====================================

template <class T, uint_t DIM>
class StaticVector
{
public:

  typedef T value_type;

protected:

  T VALUE[DIM];

public:
  StaticVector() {}
  T& operator[](const uint_t &i) { return VALUE[i]; }
  T operator[](const uint_t &i) const { return VALUE[i]; }
  uint_t size() const { return DIM; }
  T* data() { return VALUE; }
};


typedef MathVector<StaticVector<double,2> > vec2_t;
typedef MathVector<StaticVector<double,3> > vec3_t;
typedef MathVector<StaticVector<double,4> > vec4_t;
typedef MathVector<StaticVector<double,5> > vec5_t;
typedef MathVector<StaticVector<double,6> > vec6_t;

#include "mathvector_operators.h"
#include "mathvector_methods.h"

} // namespace

#endif
