#include "bsptree.hpp"
#include "stl.hpp"

#include <boost/qvm/vec.hpp>
#include <boost/qvm/vec_operations.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/qvm/map_vec_mat.hpp>
#include <boost/qvm/swizzle.hpp>
#include <boost/qvm/vec_mat_operations.hpp>
#include <boost/qvm/mat.hpp>

#include <vector>

// this example demonstrates some other features of the library

using namespace boost::qvm;

struct Vertex
{
  vec<float, 3> pos;
};

namespace bsp
{
  template <> struct bsp_vertex_traits<Vertex>
  {
    typedef vec<float, 3> position_type;
    static const position_type & getPosition(const Vertex & v)
    {
      return v.pos;
    }
    static Vertex getInterpolatedVertex(const Vertex & a, const Vertex & b, float i)
    {
      return { a.pos*(1-i) + b.pos*i };
    }

    // because we use the transform functionality, we need to add this
    // function to the trait
    static void transform(Vertex & v, const boost::qvm::mat<float, 4, 4> & m)
    {
      boost::qvm::vec<float, 4> pp = boost::qvm::XYZ1(v.pos);

      pp = m * pp;
      pp /= boost::qvm::W(pp);
      v.pos = boost::qvm::XYZ(pp);
     }
  };
}

int main()
{
  // Example 3

  std::vector<Vertex> v3 {
    {0, 0, 1},{1, 0, 1},{1, 1, 1},    {0, 0, 1}, {1, 1, 1}, {0, 1, 1},
    {0, 0, 0},{1, 1, 0},{1, 0, 0},    {0, 0, 0}, {0, 1, 0}, {1, 1, 0},

    {0, 1, 0},{1, 1, 1},{1, 1, 0},    {0, 1, 0}, {0, 1, 1}, {1, 1, 1},
    {0, 0, 0},{1, 0, 0},{1, 0, 1},    {0, 0, 0}, {1, 0, 1}, {0, 0, 1},

    {1, 0, 0},{1, 1, 0},{1, 1, 1},    {1, 0, 0}, {1, 1, 1}, {1, 0, 1},
    {0, 0, 0},{0, 1, 1},{0, 1, 0},    {0, 0, 0}, {0, 0, 1}, {0, 1, 1},
  };


  std::vector<Vertex> v1 {
    {0, 0, 1},{1, 0, 1},{1, 2, 1},    {0, 0, 1}, {1, 2, 1}, {0, 2, 1},
    {0, 0, 0},{1, 2, 0},{1, 0, 0},    {0, 0, 0}, {0, 2, 0}, {1, 2, 0},

    {0, 2, 0},{1, 2, 1},{1, 2, 0},    {0, 2, 0}, {0, 2, 1}, {1, 2, 1},
    {0, 0, 0},{1, 0, 0},{1, 0, 1},    {0, 0, 0}, {1, 0, 1}, {0, 0, 1},

    {1, 0, 0},{1, 2, 0},{1, 2, 1},    {1, 0, 0}, {1, 2, 1}, {1, 0, 1},
    {0, 0, 0},{0, 2, 1},{0, 2, 0},    {0, 0, 0}, {0, 0, 1}, {0, 2, 1},
  };

  bsp::BspTree<std::vector<Vertex>, std::vector<uint16_t>> bsp1(std::move(v1));
  bsp::BspTree<std::vector<Vertex>, std::vector<uint16_t>> bsp3(std::move(v3));

  auto a3 = bsp3.sort(vec<float, 3>{-5, 5, 5});

  if (bsp3.isInside(vec<float, 3>{-5, 5, 5})) printf("oops 1\n");
  if (!bsp3.isInside(vec<float, 3>{0.5, 0.5, 0.5})) printf("oops 3\n");

  bsp3.transform(translation_mat(vec<double,3>{-0.5, -0.5, -0.5}));
  bsp3.transform(rot_mat<4>(vec<double,3>{-0.5, -0.5, -0.5}, 1.0));
//  bsp3.transform(translation_mat(vec<double,3>{1.5, 1.5, 1.5}));

  if (bsp3.isInside(vec<float, 3>{1, 1, 0.5})) printf("oops 4\n");
  if (!bsp3.isInside(vec<float, 3>{1, 1, 1.5})) printf("oops 5\n");

  auto u = bsp1.unify(bsp1, bsp3);
//  auto u = bsp1.intersect(bsp1, bsp3);
//  auto u = bsp1.subtract(bsp1, bsp3);
//  auto u = bsp1.subtract(bsp3, bsp1);

  printSTL(u.getVertices(), u.sort(vec<float, 3>{-5, 5, 5}));
}
