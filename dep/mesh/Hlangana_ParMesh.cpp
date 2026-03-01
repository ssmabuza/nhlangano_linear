// @HEADER
// ****************************************************************************
//                Hlangana: Copyright S. Mabuza Enterprises LLC
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hlangana_ParMesh2D.hpp"

namespace hlangana
{

  ParMesh::ParMesh(const hydrofem::Mesh &mesh,
                   const Teuchos::MpiComm<int> &comm)
  {
    // initilization
    //@{
    idx_t *elmwgt = 0;
    int wgtflag = 0;
    int numflag = 0;
    int ncon = 1;
    int ncommonnodes = 2;
    int options = 0;
    int edgecut;
    //@}

    // get MPI information
    //@{
    int size, rank;
    MPI_Comm rawComm = Teuchos::getRawMpiComm(comm);
    rank = comm.getRank();
    size = comm.getSize();
    //@}

    // distribution of elements to processes
    //@{
    idx_t elmdist[size + 1];
    int nelemperp = std::ceil(mesh.numOfElements() / Teuchos::as<double>(size));
    int count = 0;
    int toomuch = nelemperp * size - mesh.numOfElements();

    for (int proc = 0; proc < size; proc++)
    {
      elmdist[proc] = count;
      if (proc == size - toomuch)
      {
        nelemperp--;
        toomuch = 0;
      }

      count += nelemperp;
      if (count >= int(mesh.numOfElements()))
        count = mesh.numOfElements();
    }

    elmdist[size] = mesh.numOfElements();
    int nparts = size;
    //@}

    // create local elements
    //@{
    int n_elems_init_loc = elmdist[rank + 1] - elmdist[rank];
    int tot_num_points = 0;
    for (int i = 0; i < n_elems_init_loc; i++)
      tot_num_points += mesh.getElement(elmdist[rank] + i).getNodes().size();

    idx_t eptr[n_elems_init_loc + 1];
    idx_t eind[tot_num_points];
    for (int i = 0; i < n_elems_init_loc; i++)
    {
      int loc_num_points = mesh.getElement(elmdist[rank] + i).getNodes().size();
      eptr[i] = i * loc_num_points;
      for (int point_ind(0); point_ind < loc_num_points; ++point_ind)
        eind[i * loc_num_points + point_ind] = mesh.getElement(elmdist[rank] + i).getNode(point_ind);
    }
    eptr[n_elems_init_loc] = tot_num_points;

    float tpwgts[ncon * nparts];
    for (int proc = 0; proc < size; proc++)
    {
      tpwgts[proc] = 1.0 / (double)size;
    }

    float ubvec[ncon];
    ubvec[0] = 1.05;

    idx_t part[n_elems_init_loc];
    //@}

    ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, elmwgt,
                             &wgtflag, &numflag, &ncon, &ncommonnodes,
                             &nparts, tpwgts, ubvec, &options,
                             &edgecut, part, &rawComm);

    // set number of global elements in mesh
    _n_elems = mesh.numOfElements();

    // determine number of local nodes
    int n_elems_proc[size];
    for (int proc = 0; proc < size; proc++)
    {
      n_elems_proc[proc] = 0;
    }

    for (int i = 0; i < n_elems_init_loc; i++)
    {
      n_elems_proc[part[i]]++;
    }

    int n_elems_fromproc[size];
    MPI_Alltoall(n_elems_proc, 1, MPI_INT, n_elems_fromproc, 1, MPI_INT, rawComm);

    _n_elems_local = 0;
    for (int proc = 0; proc < size; proc++)
    {
      _n_elems_local += n_elems_fromproc[proc];
    }

    // send the indices to the nodes
    std::vector<bool> handled;
    handled.reserve(n_elems_init_loc);

    for (int i = 0; i < n_elems_init_loc; i++)
    {
      handled.emplace_back(false);
    }

    buf.reserve(n_elems_init_loc);
    for (int i = 0; i < n_elems_init_loc; i++)
    {
      if (!handled[i])
      {
        int n = 0;
        int p = part[i];
        for (int j = i; j < n_elems_init_loc; j++)
        {
          if (part[j] == p)
          {
            buf.push_back(elmdist[rank] + j);
            handled[j] = true;
            n++;
          }
        }

        if (p != rank)
        {
          MPI_Request request;
          int rc;
          rc = MPI_Isend(&buf[buf.size() - n], n, MPI_INT, p, 0, rawComm, &request);
          assert(rc == MPI_SUCCESS);
        }
        else
        {
          for (std::size_t j = buf.size() - n; j < buf.size(); j++)
          {
            loc2glob.push_back(buf[j]);
          }
        }
      }
    }

    while (loc2glob.size() != std::size_t(_n_elems_local))
    {
      int buf[_n_elems_local];
      MPI_Status status;
      MPI_Recv(buf, _n_elems_local, MPI_INT, MPI_ANY_SOURCE, 0, rawComm, &status);
      int count;
      MPI_Get_count(&status, MPI_INT, &count);
      for (int i = 0; i < count; i++)
      {
        loc2glob.push_back(buf[i]);
      }
    }
    std::sort(loc2glob.begin(), loc2glob.end());

    // global to local mapping; should be avoided
    int ntlocals[size];
    MPI_Allgather(&_n_elems_local, 1, MPI_INT, ntlocals, 1, MPI_INT, rawComm);

    int locs2glob[_n_elems];
    int displs[size];
    displs[0] = 0;
    for (int proc = 1; proc < size; proc++)
    {
      displs[proc] = displs[proc - 1] + ntlocals[proc - 1];
    }

    displown = displs[rank];
    MPI_Allgatherv(&loc2glob[0], loc2glob.size(), MPI_INT, locs2glob, ntlocals, displs, MPI_INT, rawComm);

    glob2locs.resize(_n_elems);
    glob2loc.resize(_n_elems);
    glob2proc.resize(_n_elems);

    int idx = 0;
    for (int proc = 0; proc < size; proc++)
    {
      for (int i = 0; i < ntlocals[proc]; ++i)
      {
        glob2locs[locs2glob[idx]] = i;
        glob2proc[locs2glob[idx]] = proc;
        if (proc == rank)
        {
          glob2loc[locs2glob[idx]] = i;
        }
        idx++;
      }
    }

    std::vector<bool> needghost;
    needghost.reserve(_n_elems);
    for (int i = 0; i < _n_elems; i++)
    {
      needghost.emplace_back(false);
    }

    nghosts.resize(_n_elems_local);
    for (int i = 0; i < _n_elems_local; i++)
    {
      nghosts[i] = 0;
    }

    // FIXME: change the number of ghosts to include vertex connected elements
    for (int i = 0; i < _n_elems_local; i++)
    {
      int iglob = loc2glob[i];
      for (std::size_t k(0); k < mesh.getElement(iglob).getNodes().size(); ++k)
      {
        int pointInd = mesh.getElement(iglob).getNode(k);
        // check if the elements are already in the current proc
        // if yes, skip the if no, add then and increase number of ghosted elements and set glob2loc indexes
        auto &node_stencil_ = mesh->m_stencil[pointInd];
        for (auto it = node_stencil_.begin(); it != node_stencil_.end(); ++it)
        {
          int gg_ind = it->first;
          if (glob2proc[gg_ind] != rank)
          {
            if (needghost[gg_ind] == false)
            {
              needghost[gg_ind] = true;
              nghosts[i]++;
            }
          }
        }
      }
    }

    ghosts.reserve(_n_elems);
    for (int i = 0; i < _n_elems; i++)
    {
      if (needghost[i])
      {
        ghosts.push_back(i);
      }
    }

    for (std::size_t i = 0; i < ghosts.size(); i++)
      glob2loc[ghosts[i]] = _n_elems_local + i;

    // get the global vertices of the mesh
    m_points.resize(mesh->numOfPoints());
    for (std::size_t i = 0; i < mesh->m_points.size(); i++)
    {
      m_points[i] = mesh->m_points[i];
    }

    _n_elems = _n_elems_local + ghosts.size();
    m_elems.resize(_n_elems, nullptr);
    m_neighbors.resize(_n_elems_local);

    for (int i = 0; i < _n_elems; i++)
    {
      int icell;
      if (i < _n_elems_local)
      {
        icell = loc2glob[i];
      }
      else
      {
        icell = ghosts[i - _n_elems_local];
      }

      if (mesh->m_elems[icell]->GetElementType() == typeTriangle)
        m_elems[i] = new Triangle(mesh->m_elems[icell]->m_nodes[0], mesh->m_elems[icell]->m_nodes[1], mesh->m_elems[icell]->m_nodes[2]);
      else
        m_elems[i] = new Quadrilateral(mesh->m_elems[icell]->m_nodes[0], mesh->m_elems[icell]->m_nodes[1], mesh->m_elems[icell]->m_nodes[2], mesh->m_elems[icell]->m_nodes[3]);

      m_elems[i]->area = mesh->m_elems[icell]->area;
      m_elems[i]->orientation = mesh->m_elems[icell]->orientation;
    }

    // get the local vertices of the mesh
    m_points_ind_loc.clear();
    for (int i = 0; i < _n_elems_local; i++)
    {
      for (std::size_t lp_ind = 0; lp_ind < m_elems[i]->m_nodes.size(); ++lp_ind)
      {
        int gindex = m_elems[i]->m_nodes[lp_ind];
        if (std::find(m_points_ind_loc.begin(), m_points_ind_loc.end(), gindex) == m_points_ind_loc.end())
          m_points_ind_loc.push_back(gindex);
      }
    }

    // get all vertices including ghost nodes (to create bigger Epetra map)
    for (int i = 0; i < _n_elems; i++)
    {
      for (std::size_t lp_ind = 0; lp_ind < m_elems[i]->m_nodes.size(); ++lp_ind)
      {
        int gindex = m_elems[i]->m_nodes[lp_ind];
        if (std::find(m_points_ind.begin(), m_points_ind.end(), gindex) == m_points_ind.end())
          m_points_ind.push_back(gindex);
      }
    }

    /////////////////////////////////////////////////////////////////

    // set the node to element stencil
    m_stencil.resize(m_points_ind_loc.size());
    for (std::size_t i = 0; i < m_points_ind_loc.size(); ++i)
    {
      for (auto it = mesh->m_stencil[m_points_ind_loc[i]].begin(); it != mesh->m_stencil[m_points_ind_loc[i]].end(); ++it)
      {
        int elem_ind = it->first;
        int loc_node_ind = it->second;
        m_stencil[i][glob2loc.at(elem_ind)] = loc_node_ind;
      }
    }

    /////////////////////////////////////////////////////////////////

    for (int i = 0; i < _n_elems_local; i++)
    {
      int icell = loc2glob[i];

      if (mesh->m_elems[icell]->GetElementType() == typeQuadrilateral)
        m_neighbors[i].resize(4);
      else
        m_neighbors[i].resize(3);

      if (mesh->m_neighbors[icell][0] >= 0)
      {
        m_neighbors[i][0] = glob2loc[mesh->m_neighbors[icell][0]];
      }
      else
      {
        m_neighbors[i][0] = mesh->m_neighbors[icell][0];
      }

      if (mesh->m_neighbors[icell][1] >= 0)
      {
        m_neighbors[i][1] = glob2loc[mesh->m_neighbors[icell][1]];
      }
      else
      {
        m_neighbors[i][1] = mesh->m_neighbors[icell][1];
      }

      if (mesh->m_neighbors[icell][2] >= 0)
      {
        m_neighbors[i][2] = glob2loc[mesh->m_neighbors[icell][2]];
      }
      else
      {
        m_neighbors[i][2] = mesh->m_neighbors[icell][2];
      }

      if (mesh->m_elems[icell]->GetElementType() == typeQuadrilateral)
      {
        if (mesh->m_neighbors[icell][3] >= 0)
        {
          m_neighbors[i][3] = glob2loc[mesh->m_neighbors[icell][3]];
        }
        else
        {
          m_neighbors[i][3] = mesh->m_neighbors[icell][3];
        }
      }
    }

    //////////////////////////////////////////////////////////////////

    glob2proc.clear();
    glob2locs.clear();

    //////////////////////////////////////////////////////////////////

    setEdges();
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  ParMesh2D::~ParMesh2D()
  {
    for (std::size_t elem_ind(0); elem_ind < m_elems.size(); ++elem_ind)
    {
      delete m_elems[elem_ind];
    }
    m_elems.clear();

    for (std::size_t edge_ind(0); edge_ind < m_edges.size(); ++edge_ind)
    {
      delete m_edges[edge_ind];
    }
    m_edges.clear();
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  const std::vector<Point2D> Parmesh2D::getElementVertices(const int elem_ind) const
  {
    if (m_elems.at(elem_ind)->GetElementType() == typeTriangle)
      return std::vector<Point2D>({m_points.at(m_elems.at(elem_ind)->m_nodes[0]),
                                   m_points.at(m_elems.at(elem_ind)->m_nodes[1]),
                                   m_points.at(m_elems.at(elem_ind)->m_nodes[2])});
    else if (m_elems.at(elem_ind)->GetElementType() == typeQuadrilateral)
      return std::vector<Point2D>({m_points.at(m_elems.at(elem_ind)->m_nodes[0]),
                                   m_points.at(m_elems.at(elem_ind)->m_nodes[1]),
                                   m_points.at(m_elems.at(elem_ind)->m_nodes[2]),
                                   m_points.at(m_elems.at(elem_ind)->m_nodes[3])});
    else
      throw("Error in Parmesh2D::getElementVertices!! Element type not valid.");
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  const std::vector<Point2D> Parmesh2D::getEdgeVertices(const int edge_ind) const
  {
    return std::vector<Point2D>({m_points.at(m_edges.at(edge_ind)->m_nodes[0]), m_points.at(m_edges.at(edge_ind)->m_nodes[1])});
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  void Parmesh2D::setEdges()
  {
    for (unsigned int ElemInd(0); ElemInd < m_elems.size(); ++ElemInd)
    {
      if (m_elems[ElemInd]->GetElementType() == typeTriangle)
      {
        m_elems[ElemInd]->m_edges[0] = addEdge(m_elems[ElemInd]->m_nodes[1], m_elems[ElemInd]->m_nodes[2]);
        m_elems[ElemInd]->m_edges[1] = addEdge(m_elems[ElemInd]->m_nodes[2], m_elems[ElemInd]->m_nodes[0]);
        m_elems[ElemInd]->m_edges[2] = addEdge(m_elems[ElemInd]->m_nodes[0], m_elems[ElemInd]->m_nodes[1]);
      }
      else if (m_elems[ElemInd]->GetElementType() == typeQuadrilateral)
      {
        m_elems[ElemInd]->m_edges[0] = addEdge(m_elems[ElemInd]->m_nodes[0], m_elems[ElemInd]->m_nodes[1]);
        m_elems[ElemInd]->m_edges[1] = addEdge(m_elems[ElemInd]->m_nodes[1], m_elems[ElemInd]->m_nodes[2]);
        m_elems[ElemInd]->m_edges[2] = addEdge(m_elems[ElemInd]->m_nodes[2], m_elems[ElemInd]->m_nodes[3]);
        m_elems[ElemInd]->m_edges[3] = addEdge(m_elems[ElemInd]->m_nodes[3], m_elems[ElemInd]->m_nodes[0]);
      }
      else
      {
        std::cout << "In Parmesh2D::SetEdges : unsupported element type!!!\n";
        exit(1);
      }
    }

    for (unsigned int ElemInd = 0; ElemInd < m_elems.size(); ++ElemInd)
    {
      for (unsigned int EdgeInd = 0; EdgeInd < m_elems[ElemInd]->m_edges.size(); ++EdgeInd)
      {
        int edge_ind = m_elems[ElemInd]->m_edges[EdgeInd];
        if (m_edges[edge_ind]->m_elems[0] < 0)
          m_edges[edge_ind]->m_elems[0] = ElemInd;
        else
          m_edges[edge_ind]->m_elems[1] = ElemInd;
      }
    }

    for (unsigned int EdgeInd = 0; EdgeInd < m_edges.size(); ++EdgeInd)
      if (m_edges[EdgeInd]->m_elems[1] < 0)
        m_edges[EdgeInd]->m_is_boundary = true;
    setLengthEachEdge();
  }
  /// end void Parmesh2D::setEdges

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  int ParGrid2D::addEdge(const int node1, const int node2)
  {
    int edge_ind = getIndexOfTheEdge(node1, node2);
    if (edge_ind != -1)
      return edge_ind;
    m_edges.push_back(new Edge2D(node1, node2));
    int node = node1;
    if (node > node2)
      node = node2;
    if ((int)m_search_edge.size() < node + 1)
      m_search_edge.resize(node + 1);
    m_search_edge[node].push_back(SearchEdgeStruct(node1, node2, (int)m_edges.size() - 1));
    return (int)m_edges.size() - 1;
  }
  /// end int ParGrid2D::addEdge

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  Point2D ParGrid2D::edgeEvalCenterPoint(const int edge_ind) const
  {
    return 0.5 * (m_points.at(m_edges.at(edge_ind)->m_nodes[0]) + m_points.at(m_edges.at(edge_ind)->m_nodes[1]));
  }
  /// end Point2D ParGrid2D::edgeEvalCenterPoint

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  int ParGrid2D::getIndexOfTheEdge(const int node1, const int node2)
  {
    if (m_search_edge.empty())
    {
      m_search_edge.resize(m_points.size());
      return -1;
    }
    int node = (node1 < node2) ? node1 : node2;
    if (int(m_search_edge.size()) < node + 1)
      m_search_edge.resize(node + 1);
    for (unsigned int i = 0; i < m_search_edge[node].size(); ++i)
      if (m_search_edge[node][i].IsTheSameEdge(node1, node2))
        return m_search_edge[node][i].m_edge_ind;
    return -1;
  }
  /// end int ParGrid2D::GetIndexOfTheEdge

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  Point2D ParGrid2D::getElementCenterOfMass(const int elem_ind) const
  {
    Point2D center(0.0, 0.0);
    for (unsigned int i = 0; i < m_elems[elem_ind]->m_nodes.size(); ++i)
      center += m_points[m_elems[elem_ind]->m_nodes[i]];
    center /= m_elems[elem_ind]->m_nodes.size();
    return center;
  }
  /// end Point2D ParGrid2D::ElementEvalCenterPoint

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  Point2D ParGrid2D::getNormalOfEdge(const int elem_ind, const int local_edge_ind) const
  {
    Point2D elem_center = getElementCenterOfMass(elem_ind);
    int edge_ind = m_elems[elem_ind]->m_edges[local_edge_ind];
    Point2D edge_center = getEdgeCenterPoint(edge_ind);
    Point2D test = edge_center - elem_center;
    Point2D tangent = m_points[m_edges[edge_ind]->m_nodes[0]] - m_points[m_edges[edge_ind]->m_nodes[1]];
    Point2D normal = Point2D(tangent.y, -tangent.x);
    normal.Normalize();
    if (normal * test < 0.0)
      normal.Reverse();
    return normal;
  }
  /// end void ParGrid2D::GetVectorOfTwoNormalFluxDOFPerEdge

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  double ParGrid2D::computeElementArea(const int elem_ind)
  {
    ObjType thisCellType = m_elems[elem_ind]->GetElementType();
    double area = -666.66;
    switch (thisCellType)
    {
    case typeTriangle:
    {
      area = EvalTriangleArea(m_points[m_elems[elem_ind]->m_nodes[0]],
                              m_points[m_elems[elem_ind]->m_nodes[1]],
                              m_points[m_elems[elem_ind]->m_nodes[2]]);
    }
    break;

    case typeQuadrilateral:
    {
      area = EvalTriangleArea(m_points[m_elems[elem_ind]->m_nodes[0]],
                              m_points[m_elems[elem_ind]->m_nodes[1]],
                              m_points[m_elems[elem_ind]->m_nodes[2]]) +
             EvalTriangleArea(m_points[m_elems[elem_ind]->m_nodes[2]],
                              m_points[m_elems[elem_ind]->m_nodes[3]],
                              m_points[m_elems[elem_ind]->m_nodes[0]]);
    }
    break;
    default:
    {
    }
    }
    m_elems[elem_ind]->area = area;
    return area;
  }
  /// end double ParGrid2D::ElementEvalArea

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  double ParGrid2D::computeEdgeLength(const int edge_ind)
  {
    double length = EvalDistance(m_points[m_edges[edge_ind]->m_nodes[0]], m_points[m_edges[edge_ind]->m_nodes[1]]);
    m_edges[edge_ind]->length = length;
    return length;
  }
  /// end double ParGrid2D::EdgeEvalLength

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  void ParGrid2D::setLengthEachEdge()
  {
    for (int edgeInd = 0; edgeInd < (int)m_edges.size(); ++edgeInd)
      computeEdgeLength(edgeInd);
  }
  /// end void ParGrid2D::SetLengthEachEdge

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  void ParGrid2D::setAreaEachElement()
  {
    for (int cellInd = 0; cellInd < (int)m_elems.size(); ++cellInd)
      getEdgeCenterPoint(cellInd);
  }
  /// end void ParGrid2D::SetAreaEachElement

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  const std::vector<Edge2D *> &ParGrid2D::edges() const
  {
    return m_edges;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  const std::vector<Element2D *> &ParGrid2D::elems() const
  {
    return m_elems;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  Point2D ParGrid2D::getEdgeCenterPoint(const int edge_ind) const
  {
    return edgeEvalCenterPoint(edge_ind);
    //  return 0.5*(m_points.at(m_edges.at(edge_ind)->m_nodes[0]) + m_points.at(m_edges.at(edge_ind)->m_nodes[1]));
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  double ParGrid2D::getEdgeLength(const int edge_ind) const
  {
    return getEdgeLength(edge_ind);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  double ParGrid2D::getElementArea(const int elem_ind) const
  {
    return getElementArea(elem_ind);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  int ParGrid2D::numOfElements() const
  {
    return _n_elems_local_incl;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  int ParGrid2D::numOfPoints() const
  {
    return int(m_points_ind_loc.size());
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  const std::vector<Point2D> &ParGrid2D::points() const
  {
    return m_points;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace hlangana
