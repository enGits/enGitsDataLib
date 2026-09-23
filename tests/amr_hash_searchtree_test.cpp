// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2026 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <edl/searchtree.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>

static void require(bool ok, const char* message)
{
  if (!ok) throw std::runtime_error(message);
}

// Explicitly test a power-of-two reduction on BOTH platforms: GNU's default
// prime bucket policy can conceal this portability bug. No timing thresholds.
template<class Integer> void testHash()
{
  using Key = edl::AMRIndex<Integer>;
  std::vector<size_t> buckets(32768, 0);
  std::unordered_map<Key, int> items;
  for (int i=0; i<32; ++i) for (int j=0; j<32; ++j) for (int k=0; k<32; ++k) {
    Key key(i,j,k,5);
    ++buckets[std::hash<Key>{}(key) & (buckets.size()-1)];
    items.emplace(key, i+j+k);
  }
  size_t occupied=0, longest=0;
  for (auto n : buckets) { occupied += n!=0; longest=std::max(longest,n); }
  std::cout << "hash occupied=" << occupied << " longest=" << longest << '\n';
  require(occupied>16000 && longest<32, "AMR hash has poor low-bit distribution");
  for (int i=0; i<32; ++i) for (int j=0; j<32; ++j) for (int k=0; k<32; ++k)
    require(items.at(Key(i,j,k,5)) == i+j+k, "AMR lookup changed");
  for (int level=0; level<10; ++level) {
    Key a(1,2,3,level), b(1,2,3,level);
    require(a==b && std::hash<Key>{}(a)==std::hash<Key>{}(b), "equal keys must hash equally");
  }
  if constexpr (std::is_signed_v<Integer>) {
    const Key negative(-2, 3, 4, 5);
    items.emplace(negative, 42);
    require(items.at(Key(-2, 3, 4, 5)) == 42, "negative AMR index lookup changed");
  }
}

template<class Scalar> void testSearch()
{
  using Vec = edl::MathVector<edl::StaticVector<Scalar,3>>;
  // Match DrNUM's int index type, not just EDL's uint16_t point-tree alias.
  using Tree = edl::TSearchTree<int,int,Vec,edl::SearchTreeVectorCheck<Vec>>;
  std::vector<Vec> points;
  for(int i=0;i<8;++i) for(int j=0;j<8;++j) for(int k=0;k<8;++k)
    points.emplace_back(Scalar(i)*Scalar(.1),Scalar(j)*Scalar(.1),Scalar(k)*Scalar(.1));
  Tree tree;
  tree.setMaxBucketSize(8);
  tree.setMaxSearchDist(Scalar(.15));
  tree.setItems(points);
  for(size_t n=0;n<points.size();++n) {
    const Vec query = points[n] + Vec(Scalar(.013),Scalar(.017),Scalar(.019));
    int expected=-1;
    Scalar best=Scalar(.15);
    for(size_t i=0;i<points.size();++i) {
      Scalar distance=(query-points[i]).abs();
      if(distance<best) {best=distance;expected=static_cast<int>(i);}
    }
    require(tree.nearestItemIndex(query)==expected, "tree query differs from brute force");
  }
  require(tree.nearestItemIndex(Vec(10,10,10))==-1, "out-of-range search must miss");
}

int main()
{
  try {
    const auto start=std::chrono::steady_clock::now();
    testSearch<float>();
    testSearch<double>();
    testHash<int>();
    testHash<uint16_t>();
    std::cout << "EDL search/hash regression passed in "
              << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
              << " seconds\n";
    return 0;
  } catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
