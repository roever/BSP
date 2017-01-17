#include "bsptree.hpp"
#include "stl.hpp"

#include <boost/qvm/vec_access.hpp>
#include <boost/qvm/vec_traits_array.hpp>

#include <vector>

using namespace boost::qvm;

// this example demonstrates how to use none-qvm position types and a non-standard
// container

// let's start with our vertex structure, we will use a float array for this example
struct Vertex
{
  float p[3];
  // more fields ....
};

// and now teach qvm to use this array... for details please look at the documentation
// of boost::qvm

// define a type that can be used with our vertex position,
using float3=float *;

// teach qvm how to use that type
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

// and now teach the bsp library how to use our vertex type
//
namespace bsp {

  template<> struct bsp_vertex_traits<Vertex>
  {
    typedef const float3 position_type;
    static position_type getPosition(const Vertex & v) { return (const float3)v.p; }
    static Vertex getInterpolatedVertex(const Vertex & a, const Vertex & b, float i)
    {
      return Vertex { a.p[0]*(1-i) + b.p[0]*i, a.p[1]*(1-i) + b.p[1]*i, a.p[2]*(1-i) + b.p[2]*i };
    }
  };
}

// and now create our own container. For the sake of the example let's just wrap
// a large array
struct VertexContainer {
  Vertex data[100];
  int used = 0;
};

struct IndexContainer {
  uint16_t data[100];
  int used = 0;
};

// and teach the bsp library how to use these container, depending on which container
// (the one for the vertices or the one for the indices) and the functions of the BSP
// tree class you are actually using you may not need all of the functions below...
namespace bsp {

  template<> struct bsp_container_traits<VertexContainer>
  {
    typedef Vertex value_type;
    static auto get(const VertexContainer & v, size_t i) { return v.data[i]; }
    template <class F>
    static size_t appendInterpolate(VertexContainer & v, const value_type & a, const value_type & b, F f)
    {
      size_t res = v.used;
      v.data[v.used] = bsp_vertex_traits<value_type>::getInterpolatedVertex(a, b, f);
      v.used++;
      return res;
    }

    static void append(VertexContainer & v, const value_type & val)
    {
      v.data[v.used] = val;
      v.used++;
    }
    static void append(VertexContainer & v, const VertexContainer & v2)
    {
      for (int i = 0; i < v2.used; i++)
      {
        v.data[v.used] = v2.data[i];
        v.used++;
      }
    }
    static size_t getSize(const VertexContainer & v)
    {
      return v.used;
    }
    static void resize(VertexContainer & v, size_t i)
    {
      v.used = i;
    }
  };

  template<> struct bsp_container_traits<IndexContainer>
  {
    typedef uint16_t value_type;
    static auto get(const IndexContainer & v, size_t i) { return v.data[i]; }
    static void append(IndexContainer & v, const value_type & val)
    {
      v.data[v.used] = val;
      v.used++;
    }
    static void append(IndexContainer & v, const IndexContainer & v2)
    {
      for (int i = 0; i < v2.used; i++)
      {
        v.data[v.used] = v2.data[i];
        v.used++;
      }
    }
    static size_t getSize(const IndexContainer & v)
    {
      return v.used;
    }
    static void resize(IndexContainer & v, size_t i)
    {
      v.used = i;
    }
  };
}

using namespace boost::qvm;

int main()
{

  // initialize with 2 intersecting cubes
  VertexContainer v
  {
    {
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
    },

    72
  };

  bsp::BspTree<VertexContainer, IndexContainer> bsp(std::move(v));

  // sort when looking from the given position
  printSTL(bsp.getVertices(), bsp.sort(vec<float, 3>{-5, 5, 5}));
}
