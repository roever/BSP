#pragma once

#include <cstdint>
#include <tuple>
#include <memory>
#include <limits>
#include <cmath>
#include <vector>

#include <boost/qvm/vec.hpp>
#include <boost/qvm/vec_operations.hpp>

namespace bsp {

namespace qvm=boost::qvm;

/// specialize these 2 templates for your type V that you use within the bsp tree
/// VertexReturnType represents the type of the getPosition function
/// try to use a reference and it can be const
template <class V> struct bsp_traits;


/// A class for a bsp-Tree. The tree is meant for OpenGL usage. You input a vertex array and the tree will sort
/// them in order from back to front so that you can draw transparency properly
/// The vertex datatype is up to you. Intended usage pattern:
///  1) create a std::vector of vertices with all the vertices you want to sort
///  2) create the bsp tree from that vector, the class will take over the content of the vector created in 1)
///  3) step 2 might need to add some interpolated vertices to the vector, so get that final vector
///     with all it's data and transfer to the GPU
///  4) sort the vertices back to front as seen from multiple different angles and draw using the drawElements family
/// The class tries to minimize the number of triangles that need to be cut... only then does it try to balance the tree
/// because traversing the tree will always take the same amount of time because the number of triangles in there is
/// always the same
/// The vertex array you input may use vertex sharing, but the class will not go to the trouble of finding
/// duplicates withint he vertices it adds
//
/// \tpar V the vertex type, the type must be possible to put into a vector, it must have an interpolation constructor
///      of type V(const V & v1, const V & v2, F i) that created an interpolating vertex and there must be a getPosition
///      function (template) that allows to get the position of the vertex
/// \tpar I type of the indices to use into the vertex buffer. use the type that can contain the expected number of trignales
///      e.g uint16_t, when less than 65536 triangles can be expected
/// \tpar E exponent of the epsilon value used to find out if a vertex is on the plane of a triangle. Specify according to
///      the scale of your object, E=0 mans 1, E=1 0.1, E=2 0.01 and so on, if the vertex is less than tat amount away from
///      a plane it will be considered on the plane
/// \tpar F the floating point type used for internal position representation
template <class V, class I, int E = 4, class F = float>
class BspTree
{
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
    F distance(const std::tuple<qvm::vec<F, 3>, F> & plane, const T & t)
    {
      return qvm::dot(std::get<0>(plane), t) + std::get<1>(plane);
    }

    const F epsilon = std::pow(0.1, E);

    // calculate the sign of a number,
    int sign(F i)
    {
      if (i >  epsilon) return 1;
      if (i < -epsilon) return -1;
      return 0;
    }

    // calculate relative position of a point along a line which is a
    // units from one and b units from the other
    F relation(F a, F b)
    {
      return std::abs(a) / ( std::abs(a) + std::abs(b) );
    }

    // type for the node of the bsp-tree
    typedef struct Node {
      std::tuple<qvm::vec<F, 3>, F> plane; // the plane that intersects the space
      std::vector<I> triangles; // the triangles that are on this plane
      std::unique_ptr<struct Node> behind; // all that is behind the plane (relative to normal of plane)
      std::unique_ptr<struct Node> infront; // all that is in front of the plane
    } Node;

    // pointer to root
    std::unique_ptr<Node> root_;

    // the vertices within the tree
    V vertices_;

    // calculate the plane in hessian normal form for the triangle with the indices given in the triple p
    std::tuple<qvm::vec<F, 3>, F> calculatePlane(int a, int b, int c)
    {
      auto p1 = bsp_traits<V>::getPosition(vertices_, a);
      auto p2 = bsp_traits<V>::getPosition(vertices_, b);
      auto p3 = bsp_traits<V>::getPosition(vertices_, c);

      using boost::qvm::operator-;
      qvm::vec<F, 3> norm = qvm::normalized(qvm::cross(qvm::vref(p2)-p1, qvm::vref(p3)-p1));
      F p = -dot(norm, p1);

      return std::make_tuple(norm, p);
    }

    void append(std::vector<I> & v, I v1, I v2, I v3)
    {
      v.push_back(v1);
      v.push_back(v2);
      v.push_back(v3);
    }

    // separate the triangles within indices into the 3 lists of triangles that are behind, infront and on the
    // plane of the triangle given in pivot
    // when needed triangles are split and the part triangles are added to the proper lists
    void separateTriangles(
        int pivot,
        const std::vector<I> indices,
        std::vector<I> & behind,
        std::vector<I> & infront,
        std::vector<I> & onPlane)
    {
      auto plane = calculatePlane(indices[pivot], indices[pivot+1], indices[pivot+2]);

      for (size_t i = 0; i < indices.size(); i+=3)
      {
        std::array<F, 3> dist
        {
          distance(plane, bsp_traits<V>::getPosition(vertices_, indices[i  ])),
          distance(plane, bsp_traits<V>::getPosition(vertices_, indices[i+1])),
          distance(plane, bsp_traits<V>::getPosition(vertices_, indices[i+2]))
        };

        std::array<int, 3> side { sign(dist[0]), sign(dist[1]), sign(dist[2]) };

        std::array<I, 3> A;

        if (side[0] * side[1] == -1)
        {
          A[0] = bsp_traits<V>::addInterpolatedVertex(vertices_, indices[i  ], indices[i+1], relation(dist[0], dist[1]));
        }
        if (side[1] * side[2] == -1)
        {
          A[1] = bsp_traits<V>::addInterpolatedVertex(vertices_, indices[i+1], indices[i+2], relation(dist[1], dist[2]));
        }
        if (side[2] * side[0] == -1)
        {
          A[2] = bsp_traits<V>::addInterpolatedVertex(vertices_, indices[i+2], indices[i  ], relation(dist[2], dist[0]));
        }

        switch (splitType(side[0], side[1], side[2]))
        {
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

          case splitType( 0,  0,  0):
            append(onPlane, indices[i  ], indices[i+1], indices[i+2]);
            break;

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
    }

    // check what would happen if the plane of pivot is used as a cutting plane for the triangles in indices
    // returns the number of triangles that would end up on the plane of pivot, behind it or in front of it
    std::tuple<int, int, int> evaluatePivot(int pivot, const std::vector<I> & indices)
    {
      int cut = 0;
      int behind = 0;
      int infront = 0;

      // count how many tringles would need to be cut, would lie behind and in front of the plane
      auto plane = calculatePlane(indices[pivot], indices[pivot+1], indices[pivot+2]);

      for (size_t i = 0; i < indices.size(); i+=3)
      {
        std::array<F, 3> dist
        {
          distance(plane, bsp_traits<V>::getPosition(vertices_, i  )),
          distance(plane, bsp_traits<V>::getPosition(vertices_, i+1)),
          distance(plane, bsp_traits<V>::getPosition(vertices_, i+2))
        };

        std::array<int, 3> side { sign(dist[0]), sign(dist[1]), sign(dist[2]) };

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
    std::unique_ptr<Node> makeTree(const std::vector<I> & indices)
    {
      if (indices.size() > 0)
      {
        assert(indices.size() % 3 == 0);

        std::tuple<int, int, int> pivot = evaluatePivot(0, indices);
        size_t best = 0;

        // find a good pivot element
        for (size_t i = 3; i < indices.size(); i+=3)
        {
          auto newPivot = evaluatePivot(i, indices);

          // new pivot is better, if
          // the total number of trignales is lower (less division)
          // or equal and more triangles on the plane (less distribute into subtrees)
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

        auto node = std::make_unique<Node>();

        std::vector<I> behind, infront;

        separateTriangles(best, indices, behind, infront, node->triangles);

        node->plane = calculatePlane(indices[best], indices[best+1], indices[best+2]);
        node->behind = makeTree(behind);
        node->infront = makeTree(infront);

        return node;
      }
      else
      {
        return std::unique_ptr<Node>();
      }
    }

    // sort the triangles in the tree into the out vector so that triangles far from p are
    // in front of the output vector
    void sortBackToFront(const qvm::vec<F, 3> & p, const Node * n, std::vector<I> & out)
    {
      if (!n) return;

      if (distance(n->plane, p) < 0)
      {
        sortBackToFront(p, n->infront.get(), out);
        out.insert(out.end(), n->triangles.begin(), n->triangles.end());
        sortBackToFront(p, n->behind.get(), out);
      }
      else
      {
        sortBackToFront(p, n->behind.get(), out);
        out.insert(out.end(), n->triangles.begin(), n->triangles.end());
        sortBackToFront(p, n->infront.get(), out);
      }
    }

  public:

    // construct the tree, vertices are taken over
    // the tree is constructed in such a way that the least number
    // of triangles need to be split... this is dones greedy though
    // only after that criterium it tries to balance the triangles
    // as balancing is good for creation but not important for the more
    // sort functionality
    // vertices in input array are taken as groups of 3 forming one triangle
    BspTree(V && vertices) : vertices_(std::move(vertices))
    {
      size_t n = vertices_.size()-2;
      I i = 0;

      std::vector<I> indices;

      while (i < n)
      {
        indices.push_back(i++);
        indices.push_back(i++);
        indices.push_back(i++);
      }

      root_ = makeTree(indices);
    }

    /// same as the other constructor but this time you need to provide the indices
    /// as well, make sure, that your indices don't overstep the
    /// size of the vertives vector and also make sure that the number
    /// of elements in indices is divisible by 3. Each group of
    /// 3 elements in indices is the 3 vertices of one triangle
    BspTree(V && vertices, const std::vector<I> & indices) : vertices_(std::move(vertices))
    {
      root_ = makeTree(indices);
    }

    /// because it might be necessary to split the polygons
    /// in the list given, this list of vertices might be longer
    /// than the one initially given
    const V & getVertices() const { return vertices_; }

    /// get the index list of vertices to draw in the right order
    /// from back to fron relative ot the position given in p
    std::vector<I> sort(const qvm::vec<F, 3> & p)
    {
      std::vector<I> out;

      if (std::numeric_limits<I>::max() < vertices_.size())
      {
        // todo throw expression
      }
      else
      {
        sortBackToFront(p, root_.get(), out);
      }

      return out;
    }
};


}

