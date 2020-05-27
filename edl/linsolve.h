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
#ifndef linsolve_H
#define linsolve_H

#include "edl.h"

namespace EDL_NAMESPACE
{

struct LinSolveError
{
  double det;
  LinSolveError(double d) { det = d; }
};

// Rainers full matrix solver
template <class M, class V>
void linsolve(const M &Ain, const V &rsv, V &b)
{
  // Solve linear system Ax = rsv
  // copy vec to protect it. b will be the return value
  b = rsv;
  // copy yourself to protect matrix entries
  M A = Ain;
  
  int n = A.size();
  int k,i,j,p[n];
  double q,s,max,h,det;
  double ele_max = 0;
  
  // Find maximum element to get a relative value
  for (int i = 0; i < A.size(); ++i) {
    for (int j = 0; j < A[i].size(); ++j) {
      ele_max = std::max(ele_max,A[i][j]);
    }
  }
  
  // Get in matrix reduction
  det = 1;
  for (k = 0; k < n-1; k++) {
    max = 0.0;
    p[k] = 0;
    for(i = k; i < n; i++) {
      s=0.0;
      for(j = k; j < n; j++) {
        s = s + fabs(A[i][j]);
      }
      q = fabs(A[i][k])/s;
      if(q > max) {
        max=q;
        p[k]=i;
      }
    }
    if(!(p[k] == k)) {
      det = -det;
      for(j = 0; j < n; j++) {
        h = A[k][j];
        A[k][j] = A[p[k]][j];
        A[p[k]][j] = h;
      }
    }
    det = det*A[k][k];
    for(i = k+1; i < n; i++) {
      A[i][k] = A[i][k]/A[k][k];
      for(j = k+1; j < n; j++) {
        A[i][j] = A[i][j] - A[i][k]*A[k][j];
      }
    }
  }
  det = det*A[n-1][n-1];
  
  // Proceed with rest of system reduction
  for(k = 0; k < n-1; k++) {
    if(!(p[k]==k)) {
      h=b[k];
      b[k]=b[p[k]];
      b[p[k]]=h;
    }
  }
  for(i = 0; i < n; i++) {
    for(j = 0; j < i; j++) {
      b[i] = b[i] - A[i][j]*b[j];
    }
  }
  for(i = n-1;i >= 0; i--) {
    s = b[i];
    for(k = i+1; k < n; k++)
      s = s - A[i][k]*b[k];
    b[i] = s/A[i][i];
  }
  
  // Check Determinant and throw error, if needed
  if (fabs(det) < 1e-20) {
    throw LinSolveError(det);
  }
}

template <class M, class V>
void iterLinsolve(const M &A, const V &b, V &x, typename V::value_type tol=1e-4)
{
  typedef typename V::value_type real;
  M B, C, _A = A;
  typename V::value_type diag_max = 0;
  for (int i = 0; i < A.size(); ++i) {
    diag_max = std::max(diag_max, std::abs(A[i][i]));
  }
  for (int i = 0; i < A.size(); ++i) {
    for (int j = 0; j < A.size(); ++j) {
      if (i == j) {
        B[i][j] = absmax(real(0.1)*diag_max, A[i][j]);
        C[i][j] = A[i][j] - B[i][j];
      } else {
        B[i][j] = 0;
        C[i][j] = A[i][j];
      }
    }
  }
  B = B.inverse();
  C = B*C;
  typename V::value_type r = 0;
  x.normalise();
  int count = 0;

//  std::cout << _A << std::endl;
//  std::cout << B << std::endl;
//  std::cout << C << std::endl;

  do {
    V x_new = B*b - C*x;
    x_new.normalise();
    x_new = 0.5*x_new + 0.5*x;
    r = 0;
    for (int j = 0; j < 3; ++j) {
      r = std::max(r, std::abs(x_new[j] - x[j]));
    }
    x = x_new;
    ++count;
    if (count > 100000) {
      throw LinSolveError(_A.det());
    }
  } while (r > tol);
}

} // namespace

#endif
