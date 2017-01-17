#include <boost/qvm/vec_access.hpp>
#include <boost/qvm/vec_traits_array.hpp>

// output the resulting mesh as a stl file... we also use qvm here for simpler normal calculation
template <class C, class I>
void printSTL(const C & vertices, const I & indices)
{
  using namespace boost::qvm;

  printf("solid \n");

  for (size_t i = 0; i+2 < bsp::container_traits<I>::getSize(indices); i += 3)
  {
    const auto v1 = bsp::vertex_traits<typename bsp::container_traits<C>::value_type>::getPosition(bsp::container_traits<C>::get(vertices, bsp::container_traits<I>::get(indices, i  )));
    const auto v2 = bsp::vertex_traits<typename bsp::container_traits<C>::value_type>::getPosition(bsp::container_traits<C>::get(vertices, bsp::container_traits<I>::get(indices, i+1)));
    const auto v3 = bsp::vertex_traits<typename bsp::container_traits<C>::value_type>::getPosition(bsp::container_traits<C>::get(vertices, bsp::container_traits<I>::get(indices, i+2)));

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

