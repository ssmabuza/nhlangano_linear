
#ifndef __Hlangana_MassMatrix_HPP__
#define __Hlangana_MassMatrix_HPP__

#include <Teuchos_RCP.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "Hlangana_ParMesh.hpp"

namespace hlangana {

void EvaluateMassMatrix(const Teuchos::RCP<const ParMesh>& parMesh,
                        const Teuchos::RCP<Tpetra::CrsMatrix<double>>& massMatrix);

}

#endif /** __Hlangana_MassMatrix_HPP__ */
