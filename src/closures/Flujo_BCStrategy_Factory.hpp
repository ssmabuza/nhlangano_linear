// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_BCStrategy_Factory_HPP__
#define __Flujo_BCStrategy_Factory_HPP__

#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_Factory_Defines.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_Traits.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Flujo_BCStrategy_Dirichlet_Constant.hpp"

namespace flujo {

PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(flujo::BCStrategy_Dirichlet_Constant,
                                           BCStrategy_Dirichlet_Constant)

class BCStrategyFactory : public panzer::BCStrategyFactory {
public:
  Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits>> buildBCStrategy(
      const panzer::BC& bc,
      const Teuchos::RCP<panzer::GlobalData>& global_data) const override {
    auto bc_tm = Teuchos::rcp(new panzer::BCStrategy_TemplateManager<panzer::Traits>);
    bool found = false;

    PANZER_BUILD_BCSTRATEGY_OBJECTS("Constant", BCStrategy_Dirichlet_Constant)

    TEUCHOS_TEST_FOR_EXCEPTION(
        !found, std::logic_error,
        "Error - the BC Strategy called \"" << bc.strategy()
                                            << "\" is not supported by flujo::BCStrategyFactory.\n"
                                               "Relevant boundary condition:\n\n"
                                            << bc << std::endl);

    return bc_tm;
  }
};

}  // namespace flujo

#endif /** __Flujo_BCStrategy_Factory_HPP__ */
