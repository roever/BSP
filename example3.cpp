#include "bsptree.hpp"

#include <boost/qvm/vec_access.hpp>
#include <boost/qvm/vec_traits_array.hpp>
#include <boost/qvm/vec_mat_operations.hpp>
#include <boost/qvm/mat.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/qvm/swizzle.hpp>
#include <boost/qvm/map_vec_mat.hpp>

#include <deque>

// a third example, this time using std::arrays with 3 elements as the type to store
// the position

// this registers all std::arrays with qvm... so be careful with operator overloading... and right now is seems
// to only work with clang... not with gcc
namespace boost
{
    namespace qvm
    {
        template <class S, int D>
        struct vec_traits<std::array<S, D>>
        {
            static int const dim=D;
            typedef S scalar_type;

            template <int I> static inline scalar_type & write_element( std::array<scalar_type, dim> & v ) { return v[I]; }
            template <int I> static inline scalar_type read_element( const std::array<scalar_type, dim> & v ) { return v[I]; }

            static inline scalar_type & write_element_idx( int i, std::array<scalar_type, dim> & v ) { return v[i]; } //optional
            static inline scalar_type read_element_idx( int i, const std::array<scalar_type, dim> & v ) { return v[i]; } //optional
        };
    }
}

// our vertex class is a template and we can choose the scalar type with which we want to store the position (float, double, ...)
template <class S>
struct Vertex
{
  std::array<S, 3> pos, norm;

  Vertex(S x1, S y1, S z1) : pos {x1, y1, z1}, norm {0, 0, 0} {}
  Vertex() : pos {0, 0, 0}, norm {0, 0, 0} {}

};

namespace bsp {

  template<class S> struct bsp_vertex_traits<Vertex<S>>
  {
    typedef const std::array<S, 3> & position_type;

    static position_type getPosition(const Vertex<S> & v)
    {
      return v.pos;
    }
    static Vertex<S> getInterpolatedVertex(const Vertex<S> & a, const Vertex<S> & b, S i)
    {
      using boost::qvm::operator+;
      using boost::qvm::operator*;

      auto pos  = a.pos*(1-i) + b.pos*i;
      return Vertex<S>(boost::qvm::X(pos), boost::qvm::Y(pos), boost::qvm::Z(pos));
    }

    static void transform(Vertex<S> & v, const boost::qvm::mat<S, 4, 4> & m)
    {
      boost::qvm::vec<S, 4> pp = boost::qvm::XYZ1(v.pos);

      pp = m * pp;
      pp /= boost::qvm::W(pp);
      v.pos = boost::qvm::XYZ(pp);
    }
  };
}

template <class B>
void printSTL(const B & b)
{
  using namespace boost::qvm;

  auto a = b.sort(vec<float, 3>{-5, 5, 5});

  // output the resulting mesh as a stl file... we also use qvm here for simpler normal calculation
  printf("solid \n");

  for (size_t i = 0; i < a.size(); i += 3)
  {
    const auto v1 = bsp::bsp_vertex_traits<typename B::vertex_type>::getPosition(b.getVertices()[a[i  ]]);
    const auto v2 = bsp::bsp_vertex_traits<typename B::vertex_type>::getPosition(b.getVertices()[a[i+1]]);
    const auto v3 = bsp::bsp_vertex_traits<typename B::vertex_type>::getPosition(b.getVertices()[a[i+2]]);

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

int main()
{
  // Example 3

  std::deque<Vertex<double>> v3 {
    {0, 0, 1},{1, 0, 1},{1, 1, 1},    {0, 0, 1}, {1, 1, 1}, {0, 1, 1},
    {0, 0, 0},{1, 1, 0},{1, 0, 0},    {0, 0, 0}, {0, 1, 0}, {1, 1, 0},

    {0, 1, 0},{1, 1, 1},{1, 1, 0},    {0, 1, 0}, {0, 1, 1}, {1, 1, 1},
    {0, 0, 0},{1, 0, 0},{1, 0, 1},    {0, 0, 0}, {1, 0, 1}, {0, 0, 1},

    {1, 0, 0},{1, 1, 0},{1, 1, 1},    {1, 0, 0}, {1, 1, 1}, {1, 0, 1},
    {0, 0, 0},{0, 1, 1},{0, 1, 0},    {0, 0, 0}, {0, 0, 1}, {0, 1, 1},
  };


  std::deque<Vertex<double>> v1 {
    {0, 0, 1},{1, 0, 1},{1, 2, 1},    {0, 0, 1}, {1, 2, 1}, {0, 2, 1},
    {0, 0, 0},{1, 2, 0},{1, 0, 0},    {0, 0, 0}, {0, 2, 0}, {1, 2, 0},

    {0, 2, 0},{1, 2, 1},{1, 2, 0},    {0, 2, 0}, {0, 2, 1}, {1, 2, 1},
    {0, 0, 0},{1, 0, 0},{1, 0, 1},    {0, 0, 0}, {1, 0, 1}, {0, 0, 1},

    {1, 0, 0},{1, 2, 0},{1, 2, 1},    {1, 0, 0}, {1, 2, 1}, {1, 0, 1},
    {0, 0, 0},{0, 2, 1},{0, 2, 0},    {0, 0, 0}, {0, 0, 1}, {0, 2, 1},
  };

  std::deque<uint16_t> indices(v3.size());
  for (size_t i = 0; i < indices.size(); i++) indices[i] = i;

  bsp::BspTree<std::deque<Vertex<double>>, std::deque<uint16_t>, 4, double> bsp1(std::deque<Vertex<double>>(v1), indices);
  bsp::BspTree<std::deque<Vertex<double>>, std::deque<uint16_t>, 4, double> bsp3(std::deque<Vertex<double>>(v3), indices);

  auto a3 = bsp3.sort(boost::qvm::vec<float, 3>{-5, 5, 5});

  if (bsp3.isInside(boost::qvm::vec<float, 3>{-5, 5, 5})) printf("oops 1\n");
  if (!bsp3.isInside(boost::qvm::vec<float, 3>{0.5, 0.5, 0.5})) printf("oops 3\n");

  bsp3.transform(boost::qvm::translation_mat(boost::qvm::vec<double,3>{-0.5, -0.5, -0.5}));
  bsp3.transform(boost::qvm::rot_mat<4>(boost::qvm::vec<double,3>{-0.5, -0.5, -0.5}, 1.0));
//  bsp3.transform(boost::qvm::translation_mat(boost::qvm::vec<double,3>{1.5, 1.5, 1.5}));

  if (bsp3.isInside(boost::qvm::vec<float, 3>{1, 1, 0.5})) printf("oops 4\n");
  if (!bsp3.isInside(boost::qvm::vec<float, 3>{1, 1, 1.5})) printf("oops 5\n");

  auto u = bsp1.unify(bsp1, bsp3);
//  auto u = bsp1.intersect(bsp1, bsp3);
//  auto u = bsp1.subtract(bsp1, bsp3);
//  auto u = bsp1.subtract(bsp3, bsp1);

  printSTL(u);
}
