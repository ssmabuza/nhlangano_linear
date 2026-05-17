// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#include "Flujo_DriverFactory.hpp"

#include "Flujo_Driver.hpp"

#include <Teuchos_TestForException.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

#include <Panzer_ClosureModel_Factory_TemplateManager.hpp>

#include "Flujo_EquationSetFactory.hpp"
#include "Flujo_ClosureModel_Factory_TemplateBuilder.hpp"
#include "Flujo_BCStrategy_Factory.hpp"

namespace flujo {

Teuchos::RCP<Driver> DriverFactory::build(Teuchos::RCP<Teuchos::ParameterList> input_params,
                                          Teuchos::RCP<const Teuchos::Comm<int>> comm) {
  TEUCHOS_TEST_FOR_EXCEPTION(input_params == Teuchos::null, std::logic_error,
                             "Flujo::DriverFactory::build: input_params is null.");
  TEUCHOS_TEST_FOR_EXCEPTION(comm == Teuchos::null, std::logic_error,
                             "Flujo::DriverFactory::build: comm is null.");

  Teuchos::RCP<const Teuchos::MpiComm<int>> mpi_comm =
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm, true);

  Teuchos::RCP<const panzer::EquationSetFactory> eqset_factory =
      Teuchos::rcp(new EquationSetFactory);

  ClosureModelFactory_TemplateBuilder cm_builder;
  Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits>> cm_factory =
      Teuchos::rcp(new panzer::ClosureModelFactory_TemplateManager<panzer::Traits>);
  cm_factory->buildObjects(cm_builder);

  Teuchos::RCP<const panzer::BCStrategyFactory> bc_factory = Teuchos::rcp(new BCStrategyFactory);

  return Teuchos::rcp(new Driver(mpi_comm, eqset_factory, cm_factory, bc_factory));
}

}  // namespace flujo
