// @HEADER
// ****************************************************************************
//                Hlangana: Copyright S. Mabuza Enterprises LLC
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include <mpi.h>
#include <parmetis.h>

#include "Hlangana.hpp"

#include "BuildMap.hpp"
#include "ParGrid2D.hpp"
#include "Assembler.hpp"
#include "MassMatrix.hpp"
#include "ExportArray.hpp"

#include "LinearSolvers.hpp"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
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

#include <Hydrofem_Mesh.hpp>
#include <Hydrofem_Mesh2D.hpp>

int main(int argc, char **argv)
{

  // initialize MPI
  Teuchos::GlobalMPISession mpisess(&argc, &argv, &std::cout);
  Teuchos::RCP<Teuchos::MpiComm<int>> comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  const int myRank = comm->getRank();
  bool success = true;
  bool verbose = true;

  try
  {
    // The problem is defined on a 2D grid, global size is nx * nx.
    Teuchos::RCP<hydrofem::Mesh> mesh = Teuchos::rcp(new hydrofem::Mesh2D());
    //    mesh->GenerateTriangulatedUniformRectangularMesh(4, 4, 5.0, 5.0, typeUnionJack);
    if (myRank == 0)
      std::cout << "Serial mesh with " << grid->m_elems.size() << " elements." << std::endl;

    // build the parallel mesh
    Teuchos::RCP<ParGrid2D> par_grid = Teuchos::rcp(new ParGrid2D(grid, comm));
    //     grid.reset();

    // build all the maps needed for the problem
    Teuchos::RCP<BuildMap> mapSet = Teuchos::rcp(new BuildMap(comm, par_grid));
    {
      // building the maps
      mapSet->buildAllMaps();
    }

    if (myRank == 0)
      std::cout << "Forming global mass matrix" << std::endl;

    // build raw mass matrix on map with shared nodes
    Teuchos::RCP<Tpetra::CrsMatrix<double>> raw_mass_mat = Teuchos::rcp(new Tpetra::CrsMatrix<double>(mapSet->getConstStartingMap(), mapSet->getNumberOfOwnedPoints()));
    Teuchos::RCP<MassMatrix> massMat = Teuchos::rcp(new MassMatrix(par_grid, raw_mass_mat));

    // assemble the mass matrix in parallel
    if (myRank == 0)
      std::cout << "Assembling global mass matrix" << std::endl;
    massMat->evaluate();
    if (myRank == 0)
      std::cout << "Global mass matrix formed and assembled" << std::endl;

    comm->barrier();

    // print some info about the matrix
    massMat->getvalue()->print(std::cout);

    // build mass matrix on a 1-1 map
    Teuchos::RCP<Tpetra::CrsMatrix<double>> clean_mass_mat = Teuchos::rcp(new Tpetra::CrsMatrix<double>(mapSet->getConstOneToOneMap(), mapSet->getNumberOfOwnedPoints()));

    // create the exporter from the original map to 1-1 map
    Teuchos::RCP<ExportArray<Tpetra::CrsMatrix<double>>> exporter = Teuchos::rcp(new ExportArray<Tpetra::CrsMatrix<double>>(mapSet->getConstStartingMap(), mapSet->getConstOneToOneMap()));

    // export the matrix info
    exporter->exportData(massMat->getvalue(), clean_mass_mat);

    // clean the old info from the matrix assembly
    raw_mass_mat.reset();
    massMat.reset();

    // close the new 1-1 matrix
    clean_mass_mat->fillComplete();

    comm->barrier();
    // print the info about this new matrix
    clean_mass_mat->print(std::cout);

    // set up the sol and RHS
    Teuchos::RCP<Tpetra::Vector<>> X, B;
    // building X with the column map of A
    //     X = Teuchos::rcp(new Tpetra::Vector<>(clean_mass_mat->getColMap()));

    // building X using the 1-1 map
    X = Teuchos::rcp(new Tpetra::Vector<double>(mapSet->getConstOneToOneMap()));

    B = Teuchos::rcp(new Tpetra::Vector<>(clean_mass_mat->getDomainMap()));
    X->putScalar(0.0);
    B->putScalar(1.0);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////// Belos solver setup //////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Belos::SolverFactory<double, Tpetra::MultiVector<double>, Tpetra::Operator<double>> solverFactory;
    Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();
    solverParams->set("Num Blocks", 1);
    solverParams->set("Maximum Iterations", 1000);
    solverParams->set("Convergence Tolerance", 1.0e-08);

    Teuchos::RCP<const Tpetra::Operator<double>> A = Teuchos::rcp_dynamic_cast<const Tpetra::Operator<double>, Tpetra::CrsMatrix<double>>(clean_mass_mat);
    Teuchos::RCP<Belos::LinearProblem<double, Tpetra::MultiVector<double>, Tpetra::Operator<double>>> problem = Teuchos::rcp(new Belos::LinearProblem<double, Tpetra::MultiVector<double>, Tpetra::Operator<double>>(A, X, B));

    Teuchos::RCP<Ifpack2::Preconditioner<double>> ifpack_prec;
    Teuchos::RCP<Teuchos::ParameterList> precParams = Teuchos::parameterList();
    precParams->set("chebyshev: degree", 1);
    precParams->set("chebyshev: min eigenvalue", 0.5);
    precParams->set("chebyshev: max eigenvalue", 2.0);

    Ifpack2::Factory prec_factory;
    ifpack_prec = prec_factory.template create("CHEBYSHEV", Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<double>, Tpetra::CrsMatrix<double>>(clean_mass_mat));
    ifpack_prec->setParameters(*precParams);
    ifpack_prec->initialize();
    ifpack_prec->compute();

    Teuchos::RCP<Tpetra::Operator<double>> M = Teuchos::rcp_dynamic_cast<Tpetra::Operator<double>, Ifpack2::Preconditioner<double>>(ifpack_prec);
    problem->setRightPrec(M);

    // set the problem:
    bool set = problem->setProblem();
    if (set == false)
    {
      std::cout << "Error Belos::LinearProblem is not correctly set" << std::endl;
      exit(1);
    }

    Teuchos::RCP<Belos::SolverManager<double, Tpetra::MultiVector<double>, Tpetra::Operator<double>>> solver = solverFactory.create("Block CG", solverParams);
    solver->setProblem(problem);
    Belos::ReturnType ret = solver->solve();

    if (ret == Belos::Converged)
    {
      std::cout << "Belos converged" << std::endl;
    }
    else
    {
      std::cout << "Belos did not converge" << std::endl;
    }

    //     std::cout << "The final solution completed" << *X << std::endl;

    auto data_ = X->getData();

    //     std::cout << data_ << std::endl;

    //     mapSet->getOneToOneMap()->getGlobalElement();

    //     auto val111 = data_[0];

    int i = 0;
    for (auto it = data_.begin(); it != data_.end(); ++it)
    {
      //       X->getMap()->getLocalElement();
      std::cout << "On Proc : " << myRank << " local entry: " << i << " global entry: " << X->getMap()->getGlobalElement(i) << " Val = " << /**it*/ data_[i] << std::endl;
      ++i;
    }

    Teuchos::ArrayRCP<double> test_data;
    //     test_data.assign(i,0.0);
    test_data.resize(i, 0.0);
    test_data[0] = 100;

    //     Teuchos::ArrayRCP<double> data2_ = X->template getLocalView();

    if (myRank == 0)
    {
      std::shared_ptr<SparseMatrix> Mc = std::make_shared<SparseMatrix>(par_grid->m_points.size(), par_grid->m_points.size());
      DenseVector dX(Mc->GetNumberOfColumns(), 0.0);
      DenseVector dB(Mc->GetNumberOfColumns(), 1.0);

      Assembler::ComputeMassMatrix(*Mc.get(), *grid);

      SolveLinearSystemPCG(Mc, dB, dX);

      dX.showme();
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success)

  comm->barrier();
  if (myRank == 0)
  {
    std::cout << "**************** The end of the test reached, all tests have passed! ********************" << std::endl;
  }
  return 0;
}
