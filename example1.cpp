#include "bsptree.hpp"

#include <boost/qvm/vec.hpp>
#include <boost/qvm/vec_operations.hpp>

#include <vector>

// This is the shortest example I can come up with, it creates a BSP tree from a set of
// triangles and sorts them according to a vantage point

using namespace boost::qvm;

// This is the vertex structure we use in this example. You are free to use whatever
// structure you want, with whatever content. But for demonstration purposes we keep
// is to the essentials. The only required member for the BSP tree is the position of
// the vertices, so we stay with it.
// Also this library needs the position to be compatible with boost::qvm. We don't want
// to bribe qvm into using our own type, so we stick with the qvm types for this example
// see a later example for a possibility of your own struct
struct Vertex
{
  vec<float, 3> pos;
};

// You always have to create this specialization for your Vertex type. This allows the
// library to access the position and to create intermediate vertices, if needed. There
// are 3 members in this struct
namespace bsp
{
  template <> struct vertex_traits<Vertex>
  {
    // the position type is the type you want the library to use to handle positions.
    // optimally this is identical to your own position type in the Vertex structure,
    // but it doesn't need to be
    typedef vec<float, 3> position_type;

    // getPosition is the function used by the library to read the position field from
    // a vertex, you can return a copy or a const reference... it is up to you
    static const position_type & getPosition(const Vertex & v)
    {
      return v.pos;
    }

    // this function creates an interpolated vertex, it is intermediate between the
    // two given vertices a and b
    // the example implementation only has to care for the position, if your Vertex
    // structure contains more fields you have to interpolate them as well, or not
    // (e.g. when the colour is flat you don't)
    static Vertex getInterpolatedVertex(const Vertex & a, const Vertex & b, float i)
    {
      return { a.pos*(1-i) + b.pos*i };
    }
  };
}

int main()
{
  // simple example, mesh of a cube, the library is prepared to use vector and deque as container
  // types for the vertices and indices, but if you want to use another container you can,
  // look at the later examples to see how
  std::vector<Vertex> v2 {
    {0, 0, 0},{1, 0, 0},{1, 1, 0},    {0, 0, 0}, {1, 1, 0}, {0, 1, 0},
    {0, 0, 1},{1, 1, 1},{1, 0, 1},    {0, 0, 1}, {0, 1, 1}, {1, 1, 1},

    {0, 0, 0},{1, 0, 1},{1, 0, 0},    {0, 0, 0}, {0, 0, 1}, {1, 0, 1},
    {0, 1, 0},{1, 1, 0},{1, 1, 1},    {0, 1, 0}, {1, 1, 1}, {0, 1, 1},

    {0, 0, 0},{0, 1, 0},{0, 1, 1},    {0, 0, 0}, {0, 1, 1}, {0, 0, 1},
    {1, 0, 0},{1, 1, 1},{1, 1, 0},    {1, 0, 0}, {1, 0, 1}, {1, 1, 1},
  };

  // create the tree
  bsp::BspTree<std::vector<Vertex>, std::vector<uint8_t>> bsp2(std::move(v2));

  // sort the triangles
  auto a2 = bsp2.sort(vec<float, 3>{-5, 5, 5});
}
