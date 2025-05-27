// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <class T>
struct mv_p
{
  static T apply(const T &a, const T &b) { return a+b; };
};

template <class T>
struct mv_m
{
  static T apply(const T &a, const T &b) { return a-b; };
};

template <class T>
struct mv_ml
{
  static T apply(const T &a, const T &b) { return a*b; };
};

template <class L, class O, class R>
struct ParseNode
{
  typedef typename R::value_type value_type;
  const L &l;
  const R &r;
  ParseNode(const L &a, const R &b) : l(a), r(b) {};
  value_type operator[](const uint_t &i) const { return O::apply(l[i], r[i]); };
  uint_t size() const { return r.size(); };
  value_type abs() const;
  value_type abs2() const;
};

template <class O, class R>
struct ParseNode<double, O, R>
{
  typedef typename R::value_type value_type;
  const double l;
  const R &r;
  ParseNode(const double a, const R &b) : l(a), r(b) {};
  value_type operator[](const uint_t &i) const { return O::apply(l, r[i]); };
  uint_t size() const { return r.size(); };
  value_type abs() const;
  value_type abs2() const;
};
