
#ifndef __Hlangana_LinearSolverBelos_HPP__
#define __Hlangana_LinearSolverBelos_HPP__


#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

namespace hlangana {

/**
 * \brief A routine for solving a linear system using Belos with Ifpack2 Preconditioners
 * 
 * 
 */
void SolveLinearSystemBelos(const Teuchos::RCP<const Tpetra::CrsMatrix<double>>& A, 
                            const Teuchos::RCP<const Tpetra::Vector<double>>& rhs, 
                            const Teuchos::RCP<Tpetra::Vector<double>>& sol,
                            const Teuchos::RCP<Teuchos::ParameterList>& prec_pl,
                            const Teuchos::RCP<Teuchos::ParameterList>& solver_pl);

}

#endif /** __Hlangana_LinearSolverBelos_HPP__ */