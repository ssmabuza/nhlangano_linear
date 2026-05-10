// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_Driver_HPP__
#define __Flujo_Driver_HPP__

#include <stdexcept>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
// #include <Teuchos_DefaultOStream.hpp>

#include <Panzer_EquationSet_Factory.hpp>
#include <Panzer_ClosureModel_Factory_TemplateManager.hpp>
#include <Panzer_BCStrategy_Factory.hpp>

namespace flujo {

/**
 * @brief A base class for all driver classes.
 * 
 */
class Driver {
public:

  /**
   * @brief Constructor for the driver.
   */
  Driver(const Teuchos::RCP<const Teuchos::MpiComm<int> >& comm,
         const Teuchos::RCP<const panzer::EquationSetFactory>& eqset_factory,
         const Teuchos::RCP<const panzer::ClosureModelFactory_TemplateManager<panzer::Traits> >& cm_factory,
         const Teuchos::RCP<const panzer::BCStrategyFactory>& bc_factory
  ) : comm_(comm), eqset_factory_(eqset_factory), cm_factory_(cm_factory), bc_factory_(bc_factory) {}

  /**
   * @brief Set up the driver.
   */
  virtual void setup(const Teuchos::RCP<Teuchos::ParameterList>& input_params) {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                               "ERROR: Driver::setup not implemented!");
  }

  /**
   * @brief Solve the driver.
   */
  virtual void solve() const {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                               "ERROR: Driver::solve not implemented!");
  }

  /**
   * @brief Destructor for the driver.
   */
  virtual ~Driver() {}

private:

  Teuchos::RCP<const Teuchos::MpiComm<int>> comm_;
  Teuchos::RCP<const panzer::EquationSetFactory> eqset_factory_;
  Teuchos::RCP<const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>> cm_factory_;
  Teuchos::RCP<const panzer::BCStrategyFactory> bc_factory_;

};

} // end namespace flujo

#endif /** __Flujo_Driver_HPP__ */