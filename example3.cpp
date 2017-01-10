#include "bsptree.hpp"

#include <boost/qvm/vec_access.hpp>
#include <boost/qvm/vec_traits_array.hpp>

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

  template<class S> struct bsp_traits<std::deque<Vertex<S>>>
  {
    typedef const std::array<S, 3> & position_type;

    static position_type getPosition(const std::deque<Vertex<S>> & v, size_t i)
    {
      return v[i].pos;
    }
    static size_t addInterpolatedVertex(std::deque<Vertex<S>> & dest, size_t a, size_t b, double i)
    {
      using boost::qvm::operator+;
      using boost::qvm::operator*;

      auto pos  = dest[a].pos*(1-i) + dest[b].pos*i;
      size_t res = dest.size();
      dest.emplace_back(boost::qvm::X(pos), boost::qvm::Y(pos), boost::qvm::Z(pos));
      return res;
    }
  };
}

int main()
{
  // Example 3

  std::deque<Vertex<double>> v3 {
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


  bsp::BspTree<std::deque<Vertex<double>>, uint16_t, 4, double> bsp3(std::move(v3));

  auto a3 = bsp3.sort({-5, 5, 5});

}
