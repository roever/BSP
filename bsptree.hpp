#pragma once

#include <cstdint>
#include <tuple>
#include <memory>
#include <limits>
#include <cmath>

#include <boost/qvm/vec.hpp>
#include <boost/qvm/vec_operations.hpp>

namespace bsp {

namespace qvm=boost::qvm;

/// specialize this template for your vertex type
/// you must provide the following two functions:
/// getPosition(const V & v), which returns a qvm handleable vector of the position
/// getInterpolatedVertex(v1, v2, f) which returns an interpolated vertex between
///      the two given vertices
template <class V> struct bsp_vertex_traits;

/// specialize this template for the container that
/// you want to use, the required fields are shown in this
/// default version that works for normal std containers with random access (e.g. vector und deque)
template<class V> struct bsp_container_traits
{
  typedef typename V::value_type value_type;
  static auto get(const V & v, size_t i) { return v[i]; }
  template <class F>
  static size_t appendInterpolate(V & v, size_t a, size_t b, F f)
  {
    size_t res = v.size();
    v.emplace_back(bsp_vertex_traits<value_type>::getInterpolatedVertex(v[a], v[b], f));
    return res;
  }
  static void append(V & v, const value_type & val)
  {
    v.push_back(val);
  }
  static void append(V & v, const V & v2)
  {
    v.insert(v.end(), v2.begin(), v2.end());
  }
};

/// A class for a bsp-Tree. The tree is meant for OpenGL usage. You input container of vertices and
/// a container of indices into the vertice container and you get a bsp tree that can order your
/// polygons from back to front.
/// Limitations:
///   - triangles only
///   - coplanar triangles must not overlap, drawing order will be undefined in this case
/// The class tries to minimize the number of triangles that need to be cut... only then does it
/// try to balance the tree because traversing the tree will always take the same amount of time
/// because the number of triangles in there is always the same
//
/// \tpar C the vertex container type, you need to specialize bsp_container_traits for this, if it
///       is not a vector or deque, you also need to specialize bsp_vertex_traits for the
///       contained vertices
/// \tpar I the index container type, the contained integer types need to be big enough to index
///       all vertices
/// \tpar E exponent of the epsilon value used to find out if a vertex is on the plane of a triangle.
///       Specify according to the scale of your object, E=0 mans 1, E=1 0.1, E=2 0.01 and so on, if
///       the vertex is less than that amount away from a plane it will be considered on the plane
/// \tpar F the floating point type used for internal position representation
template <class C, class I, int E = 4, class F = float>
class BspTree
{
  public:

    // the vertex type and the type for indexing the vertex container
    using vertex_type=typename bsp_container_traits<C>::value_type;
    using index_type=typename bsp_container_traits<I>::value_type;

  private:

    // Some internal helper functions

    // a function that is used to contract the numbers between -1 and 1 into one, used
    // for the categorisation of a triangles relation to the cutting plane
    static constexpr int splitType(int a, int b, int c)
    {
      return (a+1)*16 + (b+1)*4 + (c+1);
    }

    // calculate distance of a point from a plane
    template <class T>
    F distance(const std::tuple<qvm::vec<F, 3>, F> & plane, const T & t) const
    {
      return qvm::dot(std::get<0>(plane), t) - std::get<1>(plane);
    }

    // the epsilon value for diciding if a point is on a plane
    const F epsilon = std::pow(0.1, E);

    // calculate the sign of a number,
    int sign(F i) const
    {
      if (i >  epsilon) return 1;
      if (i < -epsilon) return -1;
      return 0;
    }

    // calculate relative position of a point along a line which is a
    // units from one and b units from the other
    F relation(F a, F b) const
    {
      return std::abs(a) / ( std::abs(a) + std::abs(b) );
    }

    // type for the node of the bsp-tree
    typedef struct Node {
      std::tuple<qvm::vec<F, 3>, F> plane; // the plane that intersects the space
      I triangles;                         // the triangles that are on this plane
      std::unique_ptr<struct Node> behind; // all that is behind the plane (relative to normal of plane)
      std::unique_ptr<struct Node> infront;// all that is in front of the plane
    } Node;

    // pointer to root of bsp-tree
    std::unique_ptr<Node> root_;

    // the vertices of all triangles within the tree
    C vertices_;

    // calculate the plane in hessian normal form for the triangle with the indices given in the triple p
    std::tuple<qvm::vec<F, 3>, F> calculatePlane(int a, int b, int c)
    {
      auto p1 = bsp_vertex_traits<vertex_type>::getPosition(bsp_container_traits<C>::get(vertices_, a));
      auto p2 = bsp_vertex_traits<vertex_type>::getPosition(bsp_container_traits<C>::get(vertices_, b));
      auto p3 = bsp_vertex_traits<vertex_type>::getPosition(bsp_container_traits<C>::get(vertices_, c));

      using boost::qvm::operator-;
      qvm::vec<F, 3> norm = qvm::normalized(qvm::cross(qvm::vref(p2)-p1, qvm::vref(p3)-p1));
      F p = dot(norm, p1);

      return std::make_tuple(norm, p);
    }

    // append indices for a triangle to the index container
    void append(I & v, index_type v1, index_type v2, index_type v3)
    {
      bsp_container_traits<I>::append(v, v1);
      bsp_container_traits<I>::append(v, v2);
      bsp_container_traits<I>::append(v, v3);
    }

    // separate the triangles within indices into the 3 lists of triangles that are behind, infront and on the
    // plane of the triangle given in pivot
    // when needed triangles are split and the smaller triangles are added to the proper lists
    // return the plane
    std::tuple<qvm::vec<F, 3>, F> separateTriangles(size_t pivot, const I & indices, I & behind, I & infront, I & onPlane)
    {
      // get the plane of the pivot triangle
      auto plane = calculatePlane(indices[pivot], indices[pivot+1], indices[pivot+2]);

      // go over all triangles and separate them
      for (size_t i = 0; i < indices.size(); i+=3)
      {
        // calculate distance of the 3 vertices from the choosen partitioning plane
        std::array<F, 3> dist
        {
          distance(plane, bsp_vertex_traits<vertex_type>::getPosition(bsp_container_traits<C>::get(vertices_, indices[i  ]))),
          distance(plane, bsp_vertex_traits<vertex_type>::getPosition(bsp_container_traits<C>::get(vertices_, indices[i+1]))),
          distance(plane, bsp_vertex_traits<vertex_type>::getPosition(bsp_container_traits<C>::get(vertices_, indices[i+2])))
        };

        // check on which side of the plane the 3 points are
        std::array<int, 3> side { sign(dist[0]), sign(dist[1]), sign(dist[2]) };

        // if necessary create intermediate points for triangle
        // edges that cross the plane
        // the new points will be on the plane and will be new
        // vertices for new triangles
        // we only need to calculate the intermediate points for an
        // edge, when one vertex of the edge is on one side of the plan
        // and the other one on the other side, so when the product of
        // the signs is negative
        std::array<index_type, 3> A;

        if (side[0] * side[1] == -1)
        {
          A[0] = bsp_container_traits<C>::appendInterpolate(vertices_, indices[i  ], indices[i+1], relation(dist[0], dist[1]));
        }
        if (side[1] * side[2] == -1)
        {
          A[1] = bsp_container_traits<C>::appendInterpolate(vertices_, indices[i+1], indices[i+2], relation(dist[1], dist[2]));
        }
        if (side[2] * side[0] == -1)
        {
          A[2] = bsp_container_traits<C>::appendInterpolate(vertices_, indices[i+2], indices[i  ], relation(dist[2], dist[0]));
        }

        // go over all possible positions of the 3 vertices relative to the plane
        switch (splitType(side[0], side[1], side[2]))
        {
          // all point on one side of the plane (or on the plane)
          // in this case we simply add the complete trigangle
          // to the proper halve of the subtree
          case splitType(-1, -1, -1):
          case splitType(-1, -1,  0):
          case splitType(-1,  0, -1):
          case splitType(-1,  0,  0):
          case splitType( 0, -1, -1):
          case splitType( 0, -1,  0):
          case splitType( 0,  0, -1):
            append(behind, indices[i  ], indices[i+1], indices[i+2]);
            break;

          case splitType( 0,  0,  1):
          case splitType( 0,  1,  0):
          case splitType( 0,  1,  1):
          case splitType( 1,  0,  0):
          case splitType( 1,  0,  1):
          case splitType( 1,  1,  0):
          case splitType( 1,  1,  1):
            append(infront, indices[i  ], indices[i+1], indices[i+2]);
            break;

          // triangle on the dividing plane
          case splitType( 0,  0,  0):
            append(onPlane, indices[i  ], indices[i+1], indices[i+2]);
            break;

          // and now all the ways that the triangle can be cut by the plane
          case splitType( 1, -1,  0):
            append(behind,  indices[i+1], indices[i+2], A[0]);
            append(infront, indices[i+2], indices[i+0], A[0]);
            break;

          case splitType(-1,  0,  1):
            append(behind,  indices[i+0], indices[i+1], A[2]);
            append(infront, indices[i+1], indices[i+2], A[2]);
            break;

          case splitType( 0,  1, -1):
            append(behind,  indices[i+2], indices[i+0], A[1]);
            append(infront, indices[i+0], indices[i+1], A[1]);
            break;

          case splitType(-1,  1,  0):
            append(behind,  indices[i+2], indices[i+0], A[0]);
            append(infront, indices[i+1], indices[i+2], A[0]);
            break;

          case splitType( 1,  0, -1):
            append(behind,  indices[i+1], indices[i+2], A[2]);
            append(infront, indices[i+0], indices[i+1], A[2]);
            break;

          case splitType( 0, -1,  1):
            append(behind,  indices[i+0], indices[i+1], A[1]);
            append(infront, indices[i+2], indices[i+0], A[1]);
            break;

          case splitType( 1, -1, -1):
            append(infront, indices[i+0], A[0],         A[2]);
            append(behind,  indices[i+1], A[2],         A[0]);
            append(behind,  indices[i+1], indices[i+2], A[2]);
            break;

          case splitType(-1,  1, -1):
            append(infront, indices[i+1], A[1],         A[0]);
            append(behind,  indices[i+2], A[0],         A[1]);
            append(behind,  indices[i+2], indices[i+0], A[0]);
            break;

          case splitType(-1, -1,  1):
            append(infront, indices[i+2], A[2],         A[1]);
            append(behind,  indices[i+0], A[1],         A[2]);
            append(behind,  indices[i+0], indices[i+1], A[1]);
            break;

          case splitType(-1,  1,  1):
            append(behind,  indices[i+0], A[0],         A[2]);
            append(infront, indices[i+1], A[2],         A[0]);
            append(infront, indices[i+1], indices[i+2], A[2]);
            break;

          case splitType( 1, -1,  1):
            append(behind,  indices[i+1], A[1],         A[0]);
            append(infront, indices[i+0], A[0],         A[1]);
            append(infront, indices[i+2], indices[i+0], A[1]);
            break;

          case splitType( 1,  1, -1):
            append(behind,  indices[i+2], A[2],         A[1]);
            append(infront, indices[i+0], A[1],         A[2]);
            append(infront, indices[i+0], indices[i+1], A[1]);
            break;

        }
      }

      return plane;
    }

    // check what would happen if the plane of pivot is used as a cutting plane for the triangles in indices
    // returns the number of triangles that would end up on the plane of pivot, behind it or in front of it
    std::tuple<int, int, int> evaluatePivot(int pivot, const I & indices)
    {
      int cut = 0;
      int behind = 0;
      int infront = 0;

      // count how many tringles would need to be cut, would lie behind and in front of the plane
      auto plane = calculatePlane(indices[pivot], indices[pivot+1], indices[pivot+2]);

      // this is a simplification of the algorithm above to just count the numbers of triangles
      for (size_t i = 0; i < indices.size(); i+=3)
      {
        std::array<int, 3> side
        {
          sign(distance(plane, bsp_vertex_traits<vertex_type>::getPosition(bsp_container_traits<C>::get(vertices_, i  )))),
          sign(distance(plane, bsp_vertex_traits<vertex_type>::getPosition(bsp_container_traits<C>::get(vertices_, i+1)))),
          sign(distance(plane, bsp_vertex_traits<vertex_type>::getPosition(bsp_container_traits<C>::get(vertices_, i+2))))
        };

        switch (splitType(side[0], side[1], side[2]))
        {
          case splitType(-1, -1, -1):
          case splitType(-1, -1,  0):
          case splitType(-1,  0, -1):
          case splitType(-1,  0,  0):
          case splitType( 0, -1, -1):
          case splitType( 0, -1,  0):
          case splitType( 0,  0, -1):
            behind++;
            break;

          case splitType( 0,  0,  1):
          case splitType( 0,  1,  0):
          case splitType( 0,  1,  1):
          case splitType( 1,  0,  0):
          case splitType( 1,  0,  1):
          case splitType( 1,  1,  0):
          case splitType( 1,  1,  1):
            infront++;
            break;

          case splitType( 0,  0,  0):
            cut++;
            break;

          case splitType(-1, -1,  1):
          case splitType(-1,  1, -1):
          case splitType( 1, -1, -1):
            behind += 2;
            infront++;
            break;

          case splitType(-1,  0,  1):
          case splitType( 1,  0, -1):
          case splitType(-1,  1,  0):
          case splitType( 1, -1,  0):
          case splitType( 0, -1,  1):
          case splitType( 0,  1, -1):
            behind++;
            infront++;
            break;

          case splitType(-1,  1,  1):
          case splitType( 1, -1,  1):
          case splitType( 1,  1, -1):
            behind ++;
            infront+=2;
            break;
        }
      }

      return std::make_tuple(cut, behind, infront);
    }

    // create the bsp tree for the triangles given in the indices vector, it returns the
    // pointer to the root of the tree
    // the function cooses a cutting plane and recursively recursively calls itself with
    // the lists of triangles that are behind and in front of the choosen plane
    std::unique_ptr<Node> makeTree(const I & indices)
    {
      if (indices.size() > 0)
      {
        // we assume triangles, so the number of vertices in the index list
        // must be divisible by 3
        assert(indices.size() % 3 == 0);

        // find a good pivot element
        std::tuple<int, int, int> pivot = evaluatePivot(0, indices);
        size_t best = 0;

        for (size_t i = 3; i < indices.size(); i+=3)
        {
          auto newPivot = evaluatePivot(i, indices);

          // new pivot is better, if
          // the total number of trignales is lower (less division)
          // or equal and more triangles on the plane (less to distribute into subtrees)
          // or equal and the triangles more equally distributed between left and right

          int np = std::get<0>(newPivot);
          int nb = std::get<1>(newPivot);
          int nf = std::get<2>(newPivot);
          int ns = np+nb+nf;

          int bp = std::get<0>(pivot);
          int bb = std::get<1>(pivot);
          int bf = std::get<2>(pivot);
          int bs = bp+bb+bf;

          if ( (ns < bs) ||
               ((ns == bs) && ((np > bp) ||
                               ((np == bp) && (abs(nb-nf) < abs(bb-bf)))
                              )
               )
             )
          {
            pivot = newPivot;
            best = i;
          }
        }

        // create the node for this part of the tree
        auto node = std::make_unique<Node>();

        // container for the triangle indices for the triangles in front and behind the plane
        I behind, infront;

        // sort the triangles into the 3 containers
        node->plane = separateTriangles(best, indices, behind, infront, node->triangles);
        node->behind = makeTree(behind);
        node->infront = makeTree(infront);

        return node;
      }
      else
      {
        // this tree is empty, return nullptr
        return std::unique_ptr<Node>();
      }
    }

    // sort the triangles in the tree into the out container so that triangles far from p are
    // in front of the output vector
    template <class P>
    void sortBackToFront(const P & p, const Node * n, I & out) const
    {
      if (!n) return;

      if (distance(n->plane, p) < 0)
      {
        sortBackToFront(p, n->infront.get(), out);
        bsp_container_traits<I>::append(out, n->triangles);
        sortBackToFront(p, n->behind.get(), out);
      }
      else
      {
        sortBackToFront(p, n->behind.get(), out);
        bsp_container_traits<I>::append(out, n->triangles);
        sortBackToFront(p, n->infront.get(), out);
      }
    }

    template <class P>
    bool isInside(const P & p, const Node * n)
    {
      F dist = distance(n->plane, p);

      if (dist < 0)
      {
          if (n->behind)
              return isInside(p, n->behind.get());
          else
              return true;
      }
      else
      {
          if (n->infront)
              return isInside(p, n->infront.get());
          else
              return false;
      }
    }

    void reCalculatePlanes(Node * n)
    {
        // recalculate plane for this node with the first triangle on this plane
        n->plane = calculatePlane(
                bsp_container_traits<I>::get(n->triangles, 0),
                bsp_container_traits<I>::get(n->triangles, 1),
                bsp_container_traits<I>::get(n->triangles, 2)
                );

        // do the two subtrees
        if (n->infront) return reCalculatePlanes(n->infront.get());
        if (n->behind) return reCalculatePlanes(n->behind.get());
    }

  public:

    /// construct the tree, vertices are taken over, indices not
    /// the tree is constructed in such a way that the least number
    /// of triangles need to be split, if there is a way to build this
    /// tree without splitting, it will be found
    /// \param vertices, container with vertices, will be taken over and
    ///        new vertices appended, when necessary
    /// \param indices, container with indices into the vertices, each group of
    ///        3 corresponds to one triangle
    BspTree(C && vertices, const I & indices) : vertices_(std::move(vertices))
    {
      root_ = makeTree(indices);
    }

    /// get the vertex container
    const C & getVertices() const { return vertices_; }

    /// get a container of indices for triangles so that the triangles are sorted
    /// from back to front when viewed from the given position
    /// \param p the point from where to look
    /// \return container of indices into the vertex container
    template <class P>
    I sort(const P & p) const
    {
      I out;

      // TODO do we want to check, if the container elements are big enough?
      sortBackToFront(p, root_.get(), out);

      return out;
    }

    /// check if a point is inside the polyhedron defined by the bsp-tree or outside of
    /// it. This only works, when you defined all your triangles in counter clockwise
    /// fashion when seen from the outside and also none of the polytopes that you add
    /// to the tree may intersect
    /// \tpar P type of the position vector, qvm must be able to handle it
    /// \param p position you want to check
    /// \return true, when inside of one of the polytopes
    template <class P>
    bool isInside(const P & p)
    {
        return isInside(p, root_.get());
    }

    /// transform the polygon that this tree describes by the given matrix
    /// \tpar M matrix type, you have to be able to handle whit matrix
    ///       in the bsp_vertex_traits::transform function, otherwise you can use
    ///       whatever you want
    template <class M>
    void transform(const M & m)
    {
        // transform all vertices
        for (auto & v : vertices_)
        {
            bsp_vertex_traits<vertex_type>::transform(v, m);
        }

        // traverse tree and re-calculate all planes
        reCalculatePlanes(root_.get());
    }
};


}

