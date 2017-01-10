#include "bsptree.hpp"

#include <boost/qvm/vec_access.hpp>
#include <boost/qvm/vec_traits_array.hpp>

// a simple example that uses c-type float arrays in the vertex structure to
// store the position
//
// as we need to return a proper type to the bsp tree that also encapsulates the size
// we define us a c-type float3. This type points to an array with 3 floats, but as
// c doesn't allow proper array pointers... well the knowledge is implicit
using float3=float *;

// so starting from there we need to teach qvm that float3 is a vector of 3 floats because
// the bsp tree uses qvm for its vector calculations.. this is pretty straightforward more or
// less lifted from qvm docutmentation
namespace boost
{
  namespace qvm
  {
    template <>
      struct vec_traits<float3>
      {
        static int const dim=3;
        typedef float scalar_type;

        template <int I> static inline scalar_type & write_element( float3 & v ) { return v[I]; }
        template <int I> static inline scalar_type read_element( float3 const & v ) { return v[I]; }

        static inline scalar_type & write_element_idx( int i, float3 & v ) { return v[i]; }
        static inline scalar_type read_element_idx( int i, float3 const & v ) { return v[i]; }
      };
  }
}

// our vertex structure... containing the position in p and a normal in n
struct Vertex
{
  float p[3];

  // constructors mainly for our use below
  Vertex(float x1, float y1, float z1) : p{x1, y1, z1} {}
  Vertex() : p {0, 0, 0} {}

  // this is required by the bsp-tree to create intermediate vertices
  // when an triangle needs to be split
  // when i is 0 the function must create a vertex identical to a, if i is 1 it must
  // be identical to b and inbetween it must be linearly interpolated between the two
  Vertex(const Vertex a, const Vertex b, float i)
  {
    for (int j = 0; j < 3; j++)
    {
      p[j] = (1-i)*a.p[j] + i*b.p[j];
    }
  }

};

// we need to tell the bsp tree library how to read out the position from a vertex
// struct and what type the value will have. The value needs to be qvm compatible
namespace bsp {

  template<> struct bsp_traits<std::vector<Vertex>>
  {
    typedef const float3 position_type;
    static position_type getPosition(const std::vector<Vertex> & v, size_t i) { return (const float3)v[i].p; }
    static size_t addInterpolatedVertex(std::vector<Vertex> & dest, size_t a, size_t b, float i)
    {
      size_t res = dest.size();
      dest.emplace_back(dest[a], dest[b], i);
      return res;
    }
  };
}

using namespace boost::qvm;

int main()
{

  // initialize with 2 intersecting cubes
  std::vector<Vertex> v {
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

  // create the tree
  bsp::BspTree<std::vector<Vertex>, uint16_t> bsp(std::move(v));

  // sort when looking from the given position
  auto a = bsp.sort({-5, 5, 5});

  // output the resulting mesh as a stl file... we also use qvm here for simpler normal calculation
  printf("solid \n");

  for (size_t i = 0; i < a.size(); i += 3)
  {
    auto v1 = bsp::bsp_traits<std::vector<Vertex>>::getPosition(bsp.getVertices(), a[i  ]);
    auto v2 = bsp::bsp_traits<std::vector<Vertex>>::getPosition(bsp.getVertices(), a[i+1]);
    auto v3 = bsp::bsp_traits<std::vector<Vertex>>::getPosition(bsp.getVertices(), a[i+2]);
    auto n = cross(vref(v2)-v1, vref(v3)-v1);

    printf("facet normal %f %f %f\n", X(n), Y(n), Z(n));
    printf("  outer loop\n");
    printf("    vertex %f %f %f\n", X(v1), Y(v1), Z(v1));
    printf("    vertex %f %f %f\n", X(v2), Y(v2), Z(v2));
    printf("    vertex %f %f %f\n", X(v3), Y(v3), Z(v3));
    printf("  endloop\n");
    printf("endfacet\n");
  }

  printf("endsolid\n");
}
