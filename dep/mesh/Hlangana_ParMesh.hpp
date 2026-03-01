// @HEADER
// ****************************************************************************
//                Hlangana: Copyright S. Mabuza Enterprises LLC
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hlangana_ParMesh_HPP__
#define __Hlangana_ParMesh_HPP__

#include <mpi.h>
#include <parmetis.h>
#include <Teuchos_RCP.hpp>
#include <Tpetra_Core.hpp>
#include <Hydrofem_Mesh.hpp>

#include "Hlangana.hpp"

namespace hlangana
{

class ParMesh
{

  ParMesh() = delete;

public:

  using SPoint = hydrofem::SPoint;
  using Element = hydrofem::Element;
  using Edge = hydrofem::Edge;
  using SES = hydrofem::SearchEdgeStruct;

  //! Ctor from a serial mesh with similar interface
  ParMesh(const hydrofem::Mesh &mesh,
          const Teuchos::MpiComm<int> &tcomm);

  //! Dtor
  ~ParMesh();

  //! Some get methods for computing externally used mesh info
  //@{

  //! computes the center point of an edge
  SPoint getEdgeCenterPoint(const int edge_ind) const;

  //! gets the center of mass in an element
  SPoint getElementCenterOfMass(const int elem_ind) const;

  //! gets the element vertices
  const std::vector<SPoint> getElementVertices(const int elem_ind) const;

  //! get the edge vertices
  const std::vector<SPoint> getEdgeVertices(const int edge_ind) const;

  //! get the center of an edge
  SPoint edgeEvalCenterPoint(const int edge_ind) const;

  //! computes the outward normal to an element given the edge index
  SPoint getNormalOfEdge(const int elem_ind, const int local_edge_ind) const;

  //! returns the area of the element
  double getElementArea(const int elem_ind) const;

  //! gives the length of the edges
  double getEdgeLength(const int edge_ind) const;

  //! get the number of points
  int numOfPoints() const;

  //! get the number of elements
  int numOfElements() const;

  //@}

private:
  //! The main data of the mesh
  //@{

  //! get the points of the whole mesh
  const std::vector<SPoint> &points() const;

  //! get the elements of the mesh on the proc including ghost
  const std::vector<Element *> &elems() const;

  //! get the edges of the mesh including ghosts
  const std::vector<Edge *> &edges() const;

  //! local to global index mapping.
  std::vector<int> loc2glob;

  //! including ghosts
  std::vector<int> glob2loc;

  //! all points of the grid (serial)
  std::vector<SPoint> m_points;

  //! all points of the grid on this process (indexes excluding ghosts)
  std::vector<int> m_points_ind_loc;

  //! all points of the grid on this process (indexes including ghosts)
  std::vector<int> m_points_ind;

  //! all the elements of the grid (including ghosts)
  std::vector<Element *> m_elems;

  //! all the edges of the grid on this processor (including ghosts)
  std::vector<Edge *> m_edges;

  //! neighbors sharing an edge with an element (not sharing nodes)
  std::vector<std::vector<int>> m_neighbors;

  //! node to element mapping (exlcuding ghost nodes)
  std::vector<std::map<int, int>> m_stencil;

  //! ghost elements
  std::vector<int> ghosts;

  //! number of ghost cells per local cell
  std::vector<int> nghosts;

  //! number of local neighbours+1 (the cell itsself) per local cell
  std::vector<int> nlocals;

  //! mapping from global elements to local
  std::vector<int> glob2locs;

  //! mapping from global elements to processes
  std::vector<int> glob2proc;

  //! buffer?
  std::vector<int> buf;

  //! number of local elements (not including ghosts)
  int _n_elems_local;

  //! number of local elements (including ghosts)
  int _n_elems_local_incl;

  //! number of elements in the full serial mesh
  int _n_elems;

  //! own displacement
  int displown;

  //@}

  //! Some set methods for forming the mesh data structures
  //@{

  //! computes the area of an element
  double computeElementArea(const int elem_ind);

  //! computes the edge length
  double computeEdgeLength(const int edge_ind);

  //! sets the length for all edges
  void setLengthEachEdge();

  //! computes the area for all elems
  void setAreaEachElement();

  //! Adds an edge
  int addEdge(const int node1, const int node2);

  //! sets the edges of elems in process
  void setEdges();

  //! gets index of edge (-1 if not present)
  int getIndexOfTheEdge(const int node1, const int node2);

  //@}

protected:
  //! The local edge builder
  std::vector<std::vector<SES>> m_search_edge;
};

} // namespace hlangana

#endif /** __Hlangana_ParMesh_HPP__ */