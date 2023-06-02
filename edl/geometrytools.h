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
#ifndef geometrytools_H
#define geometrytools_H

#include "mathvector.h"
#include "smallsquarematrix.h"
#include "containertricks.h"

#include <QVector>
#include <QList>
#include <QPair>

#include "edlerror.h"
#include "edl.h"

namespace EDL_NAMESPACE
{

template <class T>
T rad2deg(T rad)
{
  return rad/M_PI*180;
}

template <class T>
T deg2rad(T deg)
{
  return deg/180*M_PI;
}

template <class VEC>
void rotate(VEC g1, VEC g2, VEC g3, VEC &b, typename VEC::value_type theta)
{
  SmallSquareMatrix<typename VEC::value_type,3> g_t;
  g_t[0] = g1;
  g_t[1] = g2;
  g_t[2] = g3;
  SmallSquareMatrix<typename VEC::value_type,3> g = g_t.transp();
  VEC bb = g_t*b;
  SmallSquareMatrix<typename VEC::value_type,3> rot = mat3_t::identity();
  rot[0][0] = cos(theta);
  rot[0][1] = -sin(theta);
  rot[1][0] = sin(theta);
  rot[1][1] = cos(theta);
  VEC bbb = rot*bb;
  b = g*bbb;
}

/** Rotates vector v around vector axis by an angle theta */
template <class VEC>
VEC rotate(VEC v, VEC axis, typename VEC::value_type theta)
{
  axis.normalise();

  // transposed base of rotate system
  SmallSquareMatrix<typename VEC::value_type,3> g_t;

  // compute projection of v in axis direction
  VEC v_axis = (axis*v)*axis;

  // compute first orthogonal vector (first base vector)
  g_t[0] = v-v_axis;

  //In case of points on the rotation axis, do nothing
  if(g_t[0].abs()==0) return v;

  g_t[0].normalise();

  // second base vector is the normalised axis
  g_t[1] = axis;

  // compute second orthogonal vector (third base vector)
  g_t[2] = g_t[0].cross(g_t[1]);

  // base of rotate system
  SmallSquareMatrix<typename VEC::value_type,3> g = g_t.transp();

  // matrix for rotation around g_t[1];
  SmallSquareMatrix<typename VEC::value_type,3> rot = SmallSquareMatrix<typename VEC::value_type,3>::identity();
  rot[0][0] =  cos(theta);
  rot[0][2] =  sin(theta);
  rot[2][0] = -sin(theta);
  rot[2][2] =  cos(theta);

  // transfer v to rotate system
  VEC v_r = g_t*v;

  // rotate the vector and transfer it back
  v_r = rot*v_r;
  v = g*v_r;

  return v;
}

/** Returns a normalized vector orthogonal to v */
template <class VEC>
VEC orthogonalVector(VEC v)
{
  // get absolute values
  real xx = v[0] < 0.0 ? -v[0] : v[0];
  real yy = v[1] < 0.0 ? -v[1] : v[1];
  real zz = v[2] < 0.0 ? -v[2] : v[2];
  // switch both biggest values and set the other one to zero
  VEC u;
  if (xx < yy) {
    u = xx < zz ? VEC(0,v[2],-v[1]) : VEC(v[1],-v[0],0);
  } else {
    u = yy < zz ? VEC(-v[2],0,v[0]) : VEC(v[1],-v[0],0);
  }
  u.normalise();
  return u;
}

/** Calculates the intersection between a line and a plane.
 * @param x_straight A point of the line.
 * @param v_straight Direction of the line.
 * @param x_plane A point of the plane
 * @param n_plane Normal of the plane
 */
template <class VEC>
typename VEC::value_type intersection(VEC x_straight, VEC v_straight, VEC x_plane, VEC n_plane)
{
  typename VEC::value_type k = (x_plane*n_plane - x_straight*n_plane)/(v_straight*n_plane);
  return k;
}

/** Calculates the intersection between a line and a plane.
 * @param x_straight A point of the line.
 * @param v_straight Direction of the line.
 * @param x_plane A point of the plane
 * @param u_plane A vector of the plane
 * @param v_plane Another vector of the plane not colinear with u_plane
 */
template <class VEC>
typename VEC::value_type intersection(VEC x_straight, VEC v_straight, VEC x_plane, VEC u_plane, VEC v_plane)
{
  VEC n = u_plane.cross(v_plane);
  return intersection(x_straight,v_straight,x_plane,n);
}

/**
 * Calculates the intersection of a segment [x1,x2] with a triangle (a,b,c)
 * Note: (xi,ri) will always be set to the intersection of the line (x1,x2) with the plane (a,b,c) even if the segment does not intersect the triangle!
 * @param a Input: Triangle point 1
 * @param b Input: Triangle point 2
 * @param c Input: Triangle point 3
 * @param x1 Input: Segment point 1
 * @param x2 Input: Segment point 2
 * @param xi Output: 3D Global coordinates of the intersection point
 * @param ri Output: 3D local triangle coordinates of the intersection point
 * @param tol Input: Relative tolerance: There can only be an intersection if:
 * 0-tol*(x1-x2).abs()<=k<=1+tol*(x1-x2).abs()
 * 0-tol<=ri[0]<=1+tol
 * 0-tol<=ri[1]<=1+tol
 * @return true if an intersection point was found, else false.
*/
template <class VEC>
bool intersectEdgeAndTriangle(const VEC& a, const VEC& b, const VEC& c,
                              const VEC& x1, const VEC& x2, VEC& xi, VEC& ri, typename VEC::value_type tol=1e-4)
{
  // triangle base
  VEC g1 = b - a;
  VEC g2 = c - a;
  VEC g3 = g1.cross(g2);
  g3.normalise();

  // direction of the edge
  VEC v = x2 - x1;

  // parallel?
  if (fabs(g3*v) < 1e-6) {
    return false;
  }

  // compute intersection between straight and triangular plane
  typename VEC::value_type k = intersection(x1, v, a, g3);
  xi = x1 + k*v;

  // transform xi to triangular base
  SmallSquareMatrix<typename VEC::value_type,3> G;
  G.setColumn(0, g1);
  G.setColumn(1, g2);
  G.setColumn(2, g3);

  SmallSquareMatrix<typename VEC::value_type,3> GI = G.inverse();
  ri = xi - a;
  ri = GI*ri;

  // intersection outside of edge range?
  if (k < 0 - tol) {
    return false;
  }
  if (k > 1 + tol) {
    return false;
  }

  // intersection outside of triangle?
  if (!isInsideTriangle(vec2_t(ri[0],ri[1]),tol)) {
    return false;
  }
  return true;
}

template <class VEC>
bool isInsideTriangle(VEC r, typename VEC::value_type tol=1e-4)
{
  if (r[0] < 0-tol || 1+tol < r[0] || r[1] < 0-tol || 1+tol < r[1] || r[0]+r[1] > 1+tol) {
    return false;
  }
  return true;
}

template <class VEC>
bool isInsideCartesianBox(VEC x, VEC x1, VEC x2)
{
  for (int i = 0; i < 3; ++i) {
    if (x[i] < x1[i] || x[i] >= x2[i]) {
      return false;
    }
  }
  return true;
}

/** Calculates the intersection point M between the lines (r1,u1) and (r2,u2).
 * @param k1 Returned by reference. Verifies M = r1+k1*u1
 * @param k2 Returned by reference. Verifies M = r2+k2*u2
 * @param r1 point of line 1
 * @param r2 point of line 2
 * @param u1 direction of line 1
 * @param u2 direction of line 2
 * @return true if an intersection has been found, else false.
 */
template <class VEC>
bool intersection (typename VEC::value_type &k1, typename VEC::value_type &k2, VEC r1, VEC u1, VEC r2, VEC u2)
{
  typename VEC::value_type ave_length = .5*(sqrt(u1[0]*u1[0]+u1[1]*u1[1]) + sqrt(u2[0]*u2[0]+u2[1]*u2[1]));
  typename VEC::value_type DET = (u1[0]*u2[1]-u1[1]*u2[0]);
  if (fabs(DET) > 1e-6*ave_length) {
    k1 = -(u2[0]*r2[1]-u2[0]*r1[1]-r2[0]*u2[1]+r1[0]*u2[1])/DET;
    k2 = -(-u1[1]*r2[0]+u1[0]*r2[1]-u1[0]*r1[1]+u1[1]*r1[0])/DET;
    return true;
  } else {
    return false;
  }
}

template <class VEC>
void sliceTriangle(const std::vector<VEC> &Tin, VEC x, VEC n, std::vector<std::vector<VEC> > &Tout)
{
  VEC a = Tin[0];
  VEC b = Tin[1];
  VEC c = Tin[2];
  typename VEC::value_type kab = intersection(a,b-a,x,n);
  typename VEC::value_type kbc = intersection(b,c-b,x,n);
  typename VEC::value_type kca = intersection(c,a-c,x,n);
  bool ab_cut = ((kab >= 0) && (kab <= 1));
  bool bc_cut = ((kbc >= 0) && (kbc <= 1));
  bool ca_cut = ((kca >= 0) && (kca <= 1));
  if (ab_cut && bc_cut && ca_cut) {
    //std::cerr << "invalid triangle (SliceTriangle) A" << std::endl;
    //exit(EXIT_FAILURE);
    if      ((kab <= kbc) && (kab <= kca)) ab_cut = false;
    else if ((kbc <= kab) && (kbc <= kca)) bc_cut = false;
    else                                   ca_cut = false;
  }
  if (ab_cut && bc_cut) {
    VEC ab = a + kab*(b-a);
    VEC bc = b + kbc*(c-b);
    Tout.resize(3, std::vector<VEC>(3));
    clinit(Tout[0]) = a,ab,bc;
    clinit(Tout[1]) = ab,b,bc;
    clinit(Tout[2]) = bc,c,a;
  } else if (bc_cut && ca_cut) {
    VEC bc = b + kbc*(c-b);
    VEC ca = c + kca*(a-c);
    Tout.resize(3, std::vector<VEC>(3));
    clinit(Tout[0]) = a,bc,ca;
    clinit(Tout[1]) = a,b,bc;
    clinit(Tout[2]) = bc,c,ca;
  } else if (ca_cut && ab_cut) {
    VEC ca = c + kca*(a-c);
    VEC ab = a + kab*(b-a);
    Tout.resize(3, std::vector<VEC>(3));
    clinit(Tout[0]) = a,ab,ca;
    clinit(Tout[1]) = ab,b,ca;
    clinit(Tout[2]) = b,c,ca;
  } else {
    Tout.resize(1, std::vector<VEC>(3));
    clinit(Tout[0]) = a,b,c;
  }
}

/** Returns the volume of a tetrahedron.
 * V= v1*(v2^v3) with vi=xi-x0
 * If neg = false and V<0, it will return V=-1e99, else it returns V.
 * @param x0 point 0 of the tetrahedron
 * @param x1 point 1 of the tetrahedron
 * @param x2 point 2 of the tetrahedron
 * @param x3 point 3 of the tetrahedron
 * @param neg If neg = false and V<0, it will return V=-1e99, else it returns V.
 * @return volume of the tetrahedron
 */
template <class VEC>
typename VEC::value_type tetraVol(const VEC& x0, const VEC& x1, const VEC& x2, const VEC& x3, bool neg=false)
{
  static typename VEC::value_type f16 = 1.0/6.0;
  VEC v1(x1[0]-x0[0], x1[1]-x0[1], x1[2]-x0[2]);
  VEC v2(x2[0]-x0[0], x2[1]-x0[1], x2[2]-x0[2]);
  VEC v3(x3[0]-x0[0], x3[1]-x0[1], x3[2]-x0[2]);
  typename VEC::value_type V = v1[0]*(v2[1]*v3[2]-v2[2]*v3[1]) + v1[1]*(v2[2]*v3[0]-v2[0]*v3[2]) + v1[2]*(v2[0]*v3[1]-v2[1]*v3[0]);
  V *= f16;
  if (!neg && (V < 0)) {
    V = -1e99;
  }
  return V; //fabs(V);
}

template <class VEC>
typename VEC::value_type pyraVol(VEC x0, VEC x1, VEC x2, VEC x3, VEC x4, bool neg=false)
{
  typename VEC::value_type V1 = tetraVol(x0, x1, x3, x4, neg) + tetraVol(x1, x2, x3, x4, neg);
  typename VEC::value_type V2 = tetraVol(x0, x1, x2, x4, neg) + tetraVol(x2, x3, x0, x4, neg);
  return std::min(V1,V2);
}

template <class VEC>
typename VEC::value_type prismVol(VEC x0, VEC x1, VEC x2, VEC x3, VEC x4, VEC x5, bool neg=false)
{
  typename VEC::value_type V = 0;
  VEC p = 1.0/6.0*(x0+x1+x2+x3+x4+x5);
  V += tetraVol(x0, x2, x1, p, neg);
  V += tetraVol(x3, x4, x5, p, neg);
  V += pyraVol (x0, x1, x4, x3, p, neg);
  V += pyraVol (x1, x2, x5, x4, p, neg);
  V += pyraVol (x0, x3, x5, x2, p, neg);
  return V;
}

template <class VEC>
typename VEC::value_type hexaVol(VEC x0, VEC x1, VEC x2, VEC x3, VEC x4, VEC x5, VEC x6, VEC x7, bool neg=false)
{
  typename VEC::value_type V = 0;
  VEC p = 1.0/8.0*(x0+x1+x2+x3+x4+x5+x6+x7);
  V += pyraVol(x0, x1, x3, x2, p, neg);
  V += pyraVol(x0, x4, x5, x1, p, neg);
  V += pyraVol(x4, x6, x7, x5, p, neg);
  V += pyraVol(x2, x3, x7, x6, p, neg);
  V += pyraVol(x1, x5, x7, x3, p, neg);
  V += pyraVol(x0, x2, x6, x4, p, neg);
  return V;
}

template <class VEC>
typename VEC::value_type triArea(VEC x0, VEC x1, VEC x2)
{
  VEC a = x1-x0;
  VEC b = x2-x0;
  typename VEC::value_type A = 0.5*((a.cross(b)).abs());
  return A;
}

template <class VEC>
typename VEC::value_type quadArea(VEC x0, VEC x1, VEC x2, VEC x3)
{
  typename VEC::value_type A = 0;
  VEC p = .25*(x0+x1+x2+x3);
  A += triArea(x0,x1,p);
  A += triArea(x1,x2,p);
  A += triArea(x2,x3,p);
  A += triArea(x3,x0,p);
  return A;
}


template <class VEC>
VEC triNormal(VEC x0, VEC x1, VEC x2)
{
  VEC a = x1-x0;
  VEC b = x2-x0;
  VEC n = 0.5*(a.cross(b));
  return n;
}

template <class VEC>
VEC quadNormal(VEC x0, VEC x1, VEC x2, VEC x3)
{
  VEC n;
  clinit(n) = 0,0,0;
  VEC p = .25*(x0+x1+x2+x3);
  n += triNormal(x0,x1,p);
  n += triNormal(x1,x2,p);
  n += triNormal(x2,x3,p);
  n += triNormal(x3,x0,p);
  return n;
}

template <class VEC>
VEC turnRight(const VEC &v)
{
  VEC u;
  u[0] =  v[1];
  u[1] = -v[0];
  return u;
}

template <class VEC>
VEC turnLeft(const VEC &v)
{
  VEC u;
  u[0] = -v[1];
  u[1] =  v[0];
  return u;
}

//polygon must be numbered clockwise
template <class VEC>
bool isConvex(VEC a,VEC b,VEC c,VEC d)
{
  VEC u[4];
  u[0]=b-a;
  u[1]=c-b;
  u[2]=d-c;
  u[3]=a-d;
  
  for(int i=0;i<4;i++) {
    VEC n=u[i].cross(u[(i+1)%4]);
    if(n[2]>0) return(false);
  }
  return(true);
}

template <class T>
bool isConvex(MathVector<StaticVector<T,2>> a_2D, MathVector<StaticVector<T,2>> b_2D, MathVector<StaticVector<T,2>> c_2D, MathVector<StaticVector<T,2>> d_2D)
{
  MathVector<StaticVector<T,3>> a_3D(a_2D[0],a_2D[1]);
  MathVector<StaticVector<T,3>> b_3D(b_2D[0],b_2D[1]);
  MathVector<StaticVector<T,3>> c_3D(c_2D[0],c_2D[1]);
  MathVector<StaticVector<T,3>> d_3D(d_2D[0],d_2D[1]);
  return(isConvex(a_3D,b_3D,c_3D,d_3D));
}

/// return the angle with relation to another 3-vector
template <class VEC>
typename VEC::value_type angle(const VEC & u, const VEC & v)
{
  // return the angle w.r.t. another 3-vector
  real ptot2 = u.abs2()*v.abs2();
  if(ptot2 <= 0) {
      return 0.0;
  } else {
    real arg = (u*v)/sqrt(ptot2);
    if(arg >  1.0) arg =  1.0;
    if(arg < -1.0) arg = -1.0;
    return acos(arg);
  }
}

/** Compute the circumscribed circle of a triangle in 3D coordinates.
  * @param a first node of the triangle
  * @param b second node of the triangle
  * @param c third node of the triangle
  * @param x on return this will be the centre of the circumscribed circle
  * @param radius on return this will be the radius of the circumscribed circle
  */
template <class VEC>
void computeCircumscribedCircle(VEC a, VEC b, VEC c, VEC &x, typename VEC::value_type &radius)
{
  typename VEC::value_type la = (b-c).abs();
  typename VEC::value_type lb = (a-c).abs();
  typename VEC::value_type lc = (a-b).abs();
  typename VEC::value_type bca = sqr(la)*(sqr(lb) + sqr(lc) - sqr(la));
  typename VEC::value_type bcb = sqr(lb)*(sqr(la) + sqr(lc) - sqr(lb));
  typename VEC::value_type bcc = sqr(lc)*(sqr(la) + sqr(lb) - sqr(lc));
  typename VEC::value_type sum = bca + bcb + bcc;
  if (fabs(sum) < 1e-6) {
    x = (1.0/3.0)*(a + b + c);
    radius = 1e99;
  } else {
    x = bca*a + bcb*b + bcc*c;
    x *= 1.0/sum;
    radius = (x-a).abs();
  }
}

template <class T>
MathVector<StaticVector<T,3>> getBarycentricCoordinates(T x, T y)
{
#if __cplusplus >= 201103L
  if(std::isnan(x) || std::isinf(x) || std::isnan(y) || std::isinf(y))
#else
  if(isnan(x) || isinf(x) || isnan(y) || isinf(y))
#endif
  {
    std::cout << "x="<<x;
    std::cout << "y="<<y;
    EDL_BUG;
  }

  T x_1=0;
  T y_1=0;
  T x_2=1;
  T y_2=0;
  T x_3=0;
  T y_3=1;

  SmallSquareMatrix<T,2> M;
  M[0][0]=x_1-x_3; M[0][1]=x_2-x_3;
  M[1][0]=y_1-y_3; M[1][1]=y_2-y_3;

  if(M.det()==0) {
    std::cout << "T.det()="<<M.det();
    std::cout << M[0][0]<<M[0][1];
    std::cout << M[1][0]<<M[1][1];
    std::cout << "T[0][0]*T[1][1]-T[1][0]*T[0][1]="<<M[0][0]*M[1][1]-M[1][0]*M[0][1];
    EDL_BUG;
  }

  T lambda_1 = ((y_2-y_3)*(x-x_3)-(x_2-x_3)*(y-y_3))/(M.det());
  T lambda_2 = (-(y_1-y_3)*(x-x_3)+(x_1-x_3)*(y-y_3))/(M.det());
  T lambda_3 = 1-lambda_1-lambda_2;

  MathVector<StaticVector<T,3>> bary_coords(lambda_1,lambda_2,lambda_3);
  return bary_coords;
}

template <class VEC>
VEC intersectionOnPlane(VEC v, VEC A, VEC nA, VEC B, VEC nB)
{
  VEC u = B-A;
  v.normalise();
  v = u.abs()*v;

  typedef MathVector<StaticVector<typename VEC::value_type,2>> vec2_t;

  vec2_t p_A(0,0);
  vec2_t p_B(1,0);
  vec2_t p_nA = projectVectorOnPlane(nA,u,v);
  vec2_t p_nB = projectVectorOnPlane(nB,u,v);

  vec2_t p_tA = turnRight(p_nA);
  vec2_t p_tB = turnRight(p_nB);

  typename VEC::value_type k1, k2;
  vec2_t p_K;
  if (!intersection(k1, k2, p_A, p_tA, p_B, p_tB)) {
    p_K = 0.5*(p_A + p_B);
  } else {
    p_K = p_A + k1*p_tA;
  }

  if (p_K[0] <0 ) {
    p_K[0] = 0;
  }
  if (p_K[0] >1 ) {
    p_K[0] = 1;
  }
  VEC K = A + p_K[0]*u + p_K[1]*v;
  return K;
}

/** Projects vector V onto plane (O,i,j)
 * @param V The vector to project
 * @param i A vector of the plane
 * @param j A vector of the plane
 * @return Returns a 2D vector (x,y) so that V = x*i +y*j
 */
template <class VEC>
MathVector<StaticVector<typename VEC::value_type,2>> projectVectorOnPlane(VEC V, VEC i, VEC j)
{
  if (i.abs2()==0) EDL_BUG;
  if (j.abs2()==0) EDL_BUG;
  typename VEC::value_type x = V*i/i.abs2();
  typename VEC::value_type y = V*j/j.abs2();
  return MathVector<StaticVector<typename VEC::value_type,2>>(x,y);
}

template <class VEC>
VEC projectPointOnEdge(const VEC& M,const VEC& A, const VEC& u)
{
  if (u.abs2()==0) EDL_BUG;
  typename VEC::value_type k = ((M-A)*u)/u.abs2();
  return A + k*u;
}

template <class VEC>
VEC projectPointOnPlane(const VEC& M, const VEC& A, const VEC& N)
{
  typename VEC::value_type k = ((M-A)*N)/N.abs2();
  return( M - k*N );
}

template <class VEC>
void cart2spherical(VEC x, typename VEC::value_type& alpha, typename VEC::value_type& beta, typename VEC::value_type& r)
{
  r = x.abs();
  static const VEC ex(1,0,0);
  VEC xy(x[0],x[1],0);
  if (x[1] >= 0) {
    alpha = angle(ex, xy);
  } else {
    alpha = 2*M_PI - angle(xy, ex);
  }
  if (xy.abs2() > 0) {
    if (x[2] >= 0) {
      beta = angle(xy, x);
    } else {
      beta = -angle(xy, x);
    }
  } else {
    beta = 0.5*M_PI;
  }
}

template <class VEC>
void cart2spherical(VEC x, const VEC& x0, typename VEC::value_type& alpha, typename VEC::value_type& beta, typename VEC::value_type& r)
{
  return cart2spherical(x - x0, alpha, beta, r);
}

template <class T>
MathVector<StaticVector<T,3>> spherical2cart(T alpha, T beta, T r)
{
  return r*MathVector<StaticVector<T,3>>(cos(alpha)*cos(beta), sin(alpha)*cos(beta), sin(beta));
}

template <class VEC>
VEC spherical2cart(VEC x0, typename VEC::value_type alpha, typename VEC::value_type beta, typename VEC::value_type r)
{
  return x0 + spherical2cart(alpha, beta, r);
}

template <class C, class VEC>
void planeFit(const C &pts, VEC &x0, VEC &n, bool closed_loop = false)
{
  x0 = VEC(0,0,0);
  int N = pts.size();
  if (!closed_loop) {
    ++N;
  }
  QVector<VEC> x(N);
  for (int i = 0; i < pts.size(); ++i) {
    x[i] = pts[i];
    x0 += x[i];
  }
  if (!closed_loop) {
    x.last() = x.first();
  }
  x0 *= 1.0/pts.size();
  for (int i = 0; i < x.size(); ++i) {
    x[i] -= x0;
  }
  n = VEC(0,0,0);
  for (int i = 0; i < x.size() - 1; ++i) {
    n += x[i].cross(x[i+1]);
  }
  real a11 = 0, a12 = 0, a13 = 0;
  real a21 = 0, a22 = 0, a23 = 0;
  real a31 = 0, a32 = 0, a33 = 0;
  for (int i = 0; i < x.size() - 1; ++i) {
    a11 += x[i][0]*x[i][0];
    a12 += x[i][0]*x[i][1];
    a13 += x[i][0]*x[i][2];

    a21 += x[i][0]*x[i][1];
    a22 += x[i][1]*x[i][1];
    a23 += x[i][1]*x[i][2];

    a31 += x[i][0]*x[i][2];
    a32 += x[i][1]*x[i][2];
    a33 += x[i][2]*x[i][2];
  }
  //n[2] = n[2];
  //n[1] = -a33*n[2]/(a22*a31 - a21*a32);
  //n[0] = -(a12*n[1] + a13*n[2])/a11;
  n.normalise();
}

template <class C>
typename C::value_type polyNormal(const C &pts, bool closed_loop=false)
{
  typedef typename C::value_type vec_t;
  vec_t x0(0,0,0);
  int N = pts.size();
  if (!closed_loop) {
    ++N;
  }
  QVector<vec_t> x(N);
  for (int i = 0; i < pts.size(); ++i) {
    x[i] = pts[i];
    x0 += x[i];
  }
  if (!closed_loop) {
    x.last() = x.first();
  }
  x0 *= 1.0/pts.size();
  for (int i = 0; i < x.size(); ++i) {
    x[i] -= x0;
  }
  vec_t n(0,0,0);
  for (int i = 0; i < x.size() - 1; ++i) {
    n += x[i].cross(x[i+1]);
  }
  return n;
}

template <class VEC>
QList<VEC> orderNodesAroundCentre(QList<VEC> x3, VEC x_centre, VEC normal)
{
  VEC g1 = orthogonalVector(normal).normalised();
  VEC g2 = (normal.cross(g1)).normalised();

  QList<MathVector<StaticVector<typename VEC::value_type,2>>> x2;
  for (int i = 0; i < x3.size(); ++i) {
    VEC xc3 = x3[i] - x_centre;
    x2 << projectVectorOnPlane(xc3, g1, g2).normalised();
  }

  QList<QPair<typename VEC::value_type, int> > indices;
  for (int i = 0; i < x3.size(); ++i) {
    typename VEC::value_type x     = x2[i][0];
    typename VEC::value_type y     = x2[i][1];
    typename VEC::value_type angle = asin(y);
    //
    if (x < 0) {
      angle = typename VEC::value_type(M_PI) - angle;
    }
    //
    indices << QPair<typename VEC::value_type, int>(angle, i);
  }
  std::sort(indices.begin(), indices.end());
  QList<VEC> x_sorted;
  for (int i = 0; i < indices.size(); ++i) {
    x_sorted << x3[indices[i].second];
  }
  //
  return x_sorted;
}



/*
  checks if point p is in the in the side of the plane formed by (x0,x1,x2) as x3.
  Returns true if p and x3 are in the same side or lay in the plane, and false otherwise
 */
template <class VEC>
bool checkPointPlaneSide(const VEC& x0, const VEC& x1, const VEC& x2, const VEC& x3, const VEC& p)
{
  VEC normal, aux;
  aux = (x1 - x0);
  normal = aux.cross(x2 - x0);
  typename VEC::value_type dot_x3 = normal*(x3 - x0);
  typename VEC::value_type dot_p = normal*(p - x0);
  return ((dot_x3*dot_p)>=0);
}

/*
  returns true if the point p is inside the tetrahedron formed by x0,x1,x2 and x3.
 */
template <class VEC>
bool isPointInTetra(const VEC& x0, const VEC& x1, const VEC& x2, const VEC& x3, const VEC& p)
{
  return  checkPointPlaneSide(x0, x1, x2, x3, p) &&
          checkPointPlaneSide(x1, x2, x3, x0, p) &&
          checkPointPlaneSide(x2, x3, x0, x1, p) &&
          checkPointPlaneSide(x3, x0, x1, x2, p);   
}

template <class VEC>
void findBoundingBox(std::vector<VEC> & points, VEC & xyzmin, VEC &xyzmax)
{
  size_t n = points[0].size();
  for (size_t i =0; i<n; ++i) {
    xyzmin[i] = FLT_MAX;
    xyzmax[i] = - FLT_MAX;
  }
  for (auto p: points) {
    for (size_t i=0; i<n; i++) {
      if (xyzmin[i] > p[i]) {
        xyzmin[i] = p[i];
      }
      if (xyzmax[i] < p[i]) {
        xyzmax[i] = p[i];
      }
    }
  }
}

}

TEST_CASE("test function findBoundingBox")
{
  using namespace edl;
  typedef MathVector<StaticVector<float,3> > vec3_t;
  std::vector<vec3_t> points;
  points.push_back(vec3_t{1,0,0});
  points.push_back(vec3_t{0,5.2,0});
  points.push_back(vec3_t{2,1,14});
  points.push_back(vec3_t{-9,3,17});
  points.push_back(vec3_t{3.3,-12.4,3});
  points.push_back(vec3_t{1.8,0,-15});
  vec3_t xyzmin,xyzmax;
  edl::findBoundingBox<vec3_t>(points,xyzmin,xyzmax);
  CHECK(edl::almostEqual(xyzmin[0],float(-9.0))==true);
  CHECK(edl::almostEqual(xyzmin[1],float(-12.4))==true);
  CHECK(edl::almostEqual(xyzmin[2],float(-15))==true);
  CHECK(edl::almostEqual(xyzmax[0],float(3.3))==true);
  CHECK(edl::almostEqual(xyzmax[1],float(5.2))==true);
  CHECK(edl::almostEqual(xyzmax[2],float(17))==true);
}


#endif
