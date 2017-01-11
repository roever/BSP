#include "bsptree.hpp"

#include <boost/qvm/vec_access.hpp>
#include <boost/qvm/vec_traits_array.hpp>

#include <vector>

// similar to example 1 but this time we directly use qvm::vert as type to avoid
// registering the return type... probably the best solution if you can do it
struct Vertex
{
  boost::qvm::vec<float, 3> pos;

  Vertex(float x1, float y1, float z1) : pos {x1, y1, z1} {}
  Vertex(const boost::qvm::vec<float, 3> & p) : pos(p) {}
  Vertex() : pos {0, 0, 0} {}


};

namespace bsp {

  template <> struct bsp_vertex_traits<Vertex>
  {
    typedef const boost::qvm::vec<float, 3> & position_type;
    static position_type getPosition(const Vertex & v)
    {
      return v.pos;
    }
    static Vertex getInterpolatedVertex(const Vertex & a, const Vertex & b, float i)
    {
      using boost::qvm::operator+;
      using boost::qvm::operator*;
      return Vertex(a.pos*(1-i) + b.pos*i);
    }
  };
}

int main()
{
  std::vector<Vertex> v2 {
    {0, 0, 0},{1, 0, 0},{1, 1, 0},    {0, 0, 0}, {1, 1, 0}, {0, 1, 0},
    {0, 0, 1},{1, 1, 1},{1, 0, 1},    {0, 0, 1}, {0, 1, 1}, {1, 1, 1},

    {0, 0, 0},{1, 0, 1},{1, 0, 0},    {0, 0, 0}, {0, 0, 1}, {1, 0, 1},
    {0, 1, 0},{1, 1, 0},{1, 1, 1},    {0, 1, 0}, {1, 1, 1}, {0, 1, 1},

    {0, 0, 0},{0, 1, 0},{0, 1, 1},    {0, 0, 0}, {0, 1, 1}, {0, 0, 1},
    {1, 0, 0},{1, 1, 1},{1, 1, 0},    {1, 0, 0}, {1, 0, 1}, {1, 1, 1},

    {0.5, 0.5, 0.5},{1.5, 0.5, 0.5},{1.5, 1.5, 0.5},    {0.5, 0.5, 0.5}, {1.5, 1.5, 0.5}, {0.5, 1.5, 0.5},
    {0.5, 0.5, 1.5},{1.5, 1.5, 1.5},{1.5, 0.5, 1.5},    {0.5, 0.5, 1.5}, {0.5, 1.5, 1.5}, {1.5, 1.5, 1.5},

    {0.5, 0.5, 0.5},{1.5, 0.5, 1.5},{1.5, 0.5, 0.5},    {0.5, 0.5, 0.5}, {0.5, 0.5, 1.5}, {1.5, 0.5, 1.5},
    {0.5, 1.5, 0.5},{1.5, 1.5, 0.5},{1.5, 1.5, 1.5},    {0.5, 1.5, 0.5}, {1.5, 1.5, 1.5}, {0.5, 1.5, 1.5},

    {0.5, 0.5, 0.5},{0.5, 1.5, 0.5},{0.5, 1.5, 1.5},    {0.5, 0.5, 0.5}, {0.5, 1.5, 1.5}, {0.5, 0.5, 1.5},
    {1.5, 0.5, 0.5},{1.5, 1.5, 1.5},{1.5, 1.5, 0.5},    {1.5, 0.5, 0.5}, {1.5, 0.5, 1.5}, {1.5, 1.5, 1.5},
  };

  std::vector<uint8_t> indices(v2.size());
  for (size_t i = 0; i < indices.size(); i++) indices[i] = i;

  // as we know that we only have very few vertices, we use uint8_t indices this time... this limits us
  // to 256 vertices though but saves some memory and would require less upload to the GPU
  bsp::BspTree<std::vector<Vertex>, std::vector<uint8_t>> bsp2(std::move(v2), indices);

  auto a2 = bsp2.sort({-5, 5, 5});
}
