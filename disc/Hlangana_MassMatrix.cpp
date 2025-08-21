

#include "Hlangana_MassMatrix.hpp"

#include <Hydrofem_FEBasis.hpp>
#include <Hydrofem_FEQuadrature.hpp>


namespace hlangana {

void EvaluateMassMatrix(const Teuchos::RCP<const ParMesh>& parMesh,
                        const Teuchos::RCP<Tpetra::CrsMatrix<double>>& massMatrix)
{
  for (int elem_ind = 0; elem_ind < int(m_grid->_n_elems_local); ++elem_ind)
  {
    auto locToglobDOFs = m_grid->m_elems[elem_ind]->m_nodes;
    std::vector<Point2D> element(locToglobDOFs.size());
    for (int i = 0; i < int(locToglobDOFs.size()); ++i)
      element[i] = m_grid->m_points[locToglobDOFs[i]];
    auto phi = P1Function2D::getBasisFunctions(element);
    for (int i = 0; i < int(phi.size()); ++i)
    {
      for (int j = 0; j < int(phi.size()); ++j)
      {
        auto fun = [&i,&j,&phi](double x,double y)->double
        {
          return phi[i](x,y)*phi[j](x,y);
        };
        
        int row = locToglobDOFs[i];
        int cols (locToglobDOFs[j]);
        double vals (FEQuadrature_Dunavant::IntegrateOverTriangle<5>(element,fun));
        m_mat->insertGlobalValues(row,1,&vals,&cols);
//         int retVal = m_mat->sumIntoGlobalValues(row,1,&vals,&cols);
//         assert(retVal == 1);
      }
    }
  }
  
  // fill and complete matrix
  if (m_mat->getMap()->isOneToOne() && (m_mat->getMap()->getComm()->getSize() > 1) )
  {
    std::cout << "The map is 1-1 so, no fill complete needed" << std::endl;
  } else {
    std::cout << "The map is not 1-1, which means we have to call fillComplete" << std::endl;
    m_mat->fillComplete();
  }
  
  // at the end, set computed flag to true
  m_is_computed = true;
}

}
