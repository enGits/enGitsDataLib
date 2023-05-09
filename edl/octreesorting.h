
#ifndef OCTREE_H
#define OCTREE_H


template <typename VEC>
struct OctreeBucket 
{  
  std::vector<VEC> m_Points;
  VEC m_xyzmin, m_xyzmax;
  OctreeBucket(VEC xyzmin, VEC xyzmax) 
  {
    m_xyzmin = xyzmin;
    m_xyzmax = xyzmax;
  }
};


template <typename CONT>
class OctreeSorting 
{
  typedef typename CONT::value_type vector_t;
  typedef typename vector_t::value_type scalar_t;

private:
  std::vector<OctreeBucket<vector_t>*> m_buckets;
    void divideBucket(OctreeBucket<vector_t>* b, int n, OctreeSorting::scalar_t tol, int n_3=10, int n_2 = 5);

public:
    void dividePointsIntoBuckets(int n, CONT points, OctreeSorting::scalar_t tol, int num_points);
    OctreeSorting(){};
    std::vector<OctreeBucket<vector_t>*> getBuckets(){return m_buckets;};
};


template <typename CONT>
void OctreeSorting<CONT>::dividePointsIntoBuckets(int n, CONT points, OctreeSorting::scalar_t tol, int num_points)
{
  std::vector< OctreeBucket<vector_t>* > buckets;
  auto xyzmin = points[0];
  auto xyzmax = points[0];
  //creating bounding box around points
  for (int i_points = 0; i_points < num_points; i_points++) 
  {
    auto xyz = points[i_points];
    for (int i=0; i<3; i++) {
      if (xyz[i] < xyzmin[i]) {
        xyzmin[i] = xyz[i];
      } else if (xyz[i] > xyzmax[i]) {
        xyzmax[i] = xyz[i];
      }
    }
  }
  //to make sure we include points in the outer edges
  xyzmin[0] = xyzmin[0]-tol; xyzmin[1] = xyzmin[1]-tol; xyzmin[2] = xyzmin[2]-tol;
  xyzmax[0] = xyzmax[0]+tol; xyzmax[1] = xyzmax[1]+tol; xyzmax[2] = xyzmax[2]+tol;
  //creating octree and copying data
  OctreeBucket<vector_t>*  b;
  b = new OctreeBucket<vector_t>(xyzmin, xyzmax);
  b->m_Points.resize(num_points);
  for (size_t i = 0; i < num_points; i++) {
    b->m_Points[i] = points[i];
  }
  // //starting division
 divideBucket(b,n,tol);
}

template <typename CONT>
void OctreeSorting< CONT>::OctreeSorting::divideBucket(OctreeBucket<vector_t>* b, int n, OctreeSorting::scalar_t tol, int n_3, int n_2) 
{
  //if the bucket has less cell than the maximun, store it in the mesh
  if (b->m_Points.size() <= n) {
    if (b->m_Points.size() > 0) {   
      m_buckets.push_back(b);
    } 
  }
// //   //otherwise divide it and check the new smaller buckets
  else {
    std::vector<OctreeBucket<vector_t>*> new_buckets;
    int n_div;
    //we decide how many divisions to make
    if (b->m_Points.size() > n_3*n) {
      n_div = 3;
    } else if (b->m_Points.size() > n_2*n) {
      n_div = 2;
    } else {
      n_div = 1;
    }
// //     //check in which direction to divide
    vector_t dxdydz0 = b->m_xyzmax-b->m_xyzmin;
    vector_t dxdydz = dxdydz0;
    for (int i_div = 0; i_div < n_div; ++i_div) {
      if (dxdydz[0] >= std::max(dxdydz[1],dxdydz[2])) {
        dxdydz[0] = dxdydz[0]/2.0;
      } else if (dxdydz[1] >= std::max(dxdydz[0],dxdydz[2])) {
        dxdydz[1] = dxdydz[1]/2.0;
      } else {
        dxdydz[2] = dxdydz[2]/2.0;
      }
    }
    // //create new buckets
    int nx = round(dxdydz0[0]/dxdydz[0]);
    int ny = round(dxdydz0[1]/dxdydz[1]);
    int nz = round(dxdydz0[2]/dxdydz[2]);
    scalar_t x0 = b->m_xyzmin[0], y0 = b->m_xyzmin[1], z0 = b->m_xyzmin[2]; 
    scalar_t x=x0, y=y0, z=z0;
    for (int ix = 0; ix<nx; ix++) {
      y=y0;
      for ( int iy=0; iy<ny; iy++) {
        z=z0;
        for (int iz=0; iz<nz; iz++) {
          vector_t xyz {x,y,z};
          OctreeBucket<vector_t>*  b_new;
          b_new = new OctreeBucket<vector_t>(xyz, xyz+dxdydz);      
          new_buckets.push_back(b_new);
          z = z+dxdydz[2];
        }
        y = y+dxdydz[1];
      }
        x = x + dxdydz[0];
    }

//   //sort cells in new buckets
    for (auto point : b->m_Points) {
      for (auto b_new : new_buckets) {
        if (isInsideCartesianBox(point,b_new->m_xyzmin, b_new->m_xyzmax))
        {
          b_new->m_Points.push_back(point);
        }
      }
    }
    // try to divide new buckets
    for (auto b_new: new_buckets) {
      divideBucket(b_new, n,tol);
    }
   }
  }


#endif