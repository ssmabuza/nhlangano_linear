

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosSolverManager.hpp>
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Preconditioner.hpp>
#include "BelosTpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <BelosSolverFactory.hpp>
#include "BelosTpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_Details_PackTraits.hpp>


#include "LinearSolverBelos.hpp"



void SolveLinearSystemBelos(const Teuchos::RCP< Tpetra::CrsMatrix< double > >& A, 
                            const Teuchos::RCP< Tpetra::Vector< double > >& rhs, 
                            Teuchos::RCP< Tpetra::Vector< double > >& sol, 
                            Teuchos::RCP< Teuchos::ParameterList >& prec_pl, 
                            Teuchos::RCP< Teuchos::ParameterList >& solver_pl)
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////// Belos solver setup //////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    Belos::SolverFactory<double,Tpetra::MultiVector<double>,Tpetra::Operator<double>> solverFactory;
    Teuchos::RCP<const Tpetra::Operator<double>> tmp_A = Teuchos::rcp_dynamic_cast<const Tpetra::Operator<double>,Tpetra::CrsMatrix<double>>(A);
    Teuchos::RCP<Belos::LinearProblem<double,Tpetra::MultiVector<double>,Tpetra::Operator<double>>> problem 
               = Teuchos::rcp(new Belos::LinearProblem<double,Tpetra::MultiVector<double>,Tpetra::Operator<double>>(tmp_A,sol,rhs));
    
    Teuchos::RCP<Ifpack2::Preconditioner<double>> ifpack_prec;
    Ifpack2::Factory prec_factory;
    ifpack_prec = prec_factory.template create("CHEBYSHEV", Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<double>,Tpetra::CrsMatrix<double>>(A));
    ifpack_prec->setParameters( *prec_pl );
    ifpack_prec->initialize();
    ifpack_prec->compute();    
    
    Teuchos::RCP<Tpetra::Operator<double>> M = Teuchos::rcp_dynamic_cast<Tpetra::Operator<double>,Ifpack2::Preconditioner<double>>(ifpack_prec);
    problem->setRightPrec( M );
    
    // set the problem:
    bool set = problem->setProblem();
    if (set == false)
    {
      std::cout << "Error Belos::LinearProblem is not correctly set" << std::endl;
      exit(1);
    }
    
    Teuchos::RCP<Belos::SolverManager<double,Tpetra::MultiVector<double>,Tpetra::Operator<double>>> solver = solverFactory.create("Block CG", solver_pl);
    solver->setProblem(problem);
    Belos::ReturnType ret = solver->solve();
    
    if (ret == Belos::Converged)
    {
      std::cout << "Belos converged" << std::endl;
    } else {
      std::cout << "Belos did not converge" << std::endl;
      exit(1);
    }

}

