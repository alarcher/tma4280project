#ifndef TMA_MESH_H_
#define TMA_MESH_H_

#include <tma/types.h>

#include <iostream>
#include <iomanip>

namespace tma
{

/*
 * This is a simple topology class to store cell-vertices connectivities in
 * an array and return the number of entities.
 * To make it easier, a template argument is used: it defines the cell type.
 */
template <class T>
struct topology
{
  // The constructor requires the number of cells and vertices
  topology(uidx ncells, uidx nverts) :
    ncells_(ncells),
    nverts_(nverts),
    cv_(ncells ? new uidx[ncells * T::num(0)]() : NULL),// Should be non-zero
    offset_(0)
  {
  }

  ~topology()
  {
    delete [] cv_;
  }

  // Return the topological dimension
  uint dim() const { return T::dim(); }

  // Return the number of cells in the mesh
  uidx ncells() const { return ncells_; }

  // Return the number of vertices in the mesh
  uidx nverts() const { return nverts_; }

  // Accessor for the vertices of the i-th cell
  uidx       * operator()(uint i)       { return cv_ + i * T::num(0); }
  uidx const * operator()(uint i) const { return cv_ + i * T::num(0); }

  // Display the connectivities
  void dump() const
  {
    uidx const * cv = cv_;
    for (uidx c = 0; c < ncells_; ++c)
    {
      std::cout << std::setw(4) << c << ":";
      for (uint v = 0; v < T::num(0); ++v, ++cv)
      {
        std::cout << std::setw(4) << *cv;
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }

private:

  uidx const ncells_;
  uidx const nverts_;
  uidx * const cv_;
  uidx offset_; // Use in parallel

};



/*
 * This is a simple geometry class to store vertex coordinates in an array.
 */
template <uint D>
struct geometry
{
  geometry(uidx nverts) :
    nverts_(nverts),
    vx_(nverts ? new real[nverts*D]() : NULL)
  {
  }

  ~geometry()
  {
    delete [] vx_;
  }

  // Accessor for the coordinates of the i-th vertex
  real       * operator()(uidx i)       { return vx_ + i * D; }
  real const * operator()(uidx i) const { return vx_ + i * D; }

  // Display the coordinates
  void dump() const
  {
    real const * vx = vx_;
    for (uidx v = 0; v < nverts_; ++v)
    {
      std::cout << std::setw(4) << v << ":";
      for (uint x = 0; x < D; ++x, ++vx)
      {
        std::cout << std::setw(8) << *vx;
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }

private:

  uidx const nverts_;
  real * const vx_;

};


template <class T, uint D>
struct mesh
{

  /*
   * Now let us define a mesh class combining both ingredients to build a mesh
   * for a given cell type in a given Euclidean space.
   * We use template arguments to define the topology and the geometry and we
   * only need to specify the number of cells and vertices.
   */

  mesh(uint ncells, uint nverts) :
    T_(ncells, nverts),
    G_(nverts)
  {
  }

  uidx ncells() const { return T_.ncells(); }
  uidx nverts() const { return T_.nverts(); }

  // Topology
  topology<T>& topo() { return T_; }

  // Return array of vertex indices for the given cell
  uidx       * K(uidx i)       { return T_(i); }
  uidx const * K(uidx i) const { return T_(i); }
  

  // Geometry
  geometry<D>& geom() { return G_; }

  // Return array of coordinates for the given vertex
  real       * x(uidx i)       { return G_(i); }
  real const * x(uidx i) const { return G_(i); }

  // Display
  void dump() const
  {
    message("Topology");
    T_.dump();
    message("Geometry");
    G_.dump();
  }

private:

  topology<T> T_;
  geometry<D> G_;


};

} /* namespace tma */

#endif /* TMA_MESH_H_ */
