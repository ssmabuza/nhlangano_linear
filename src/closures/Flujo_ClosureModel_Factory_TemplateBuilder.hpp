// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_ClosureModel_Factory_TemplateBuilder_HPP__
#define __Flujo_ClosureModel_Factory_TemplateBuilder_HPP__

#include "Flujo_ClosureModel_Factory.hpp"
#include "Teuchos_RCP.hpp"

namespace flujo {

class ClosureModelFactory_TemplateBuilder {
public:
  template <typename EvalT>
  Teuchos::RCP<panzer::ClosureModelFactoryBase> build() const {
    return Teuchos::rcp(static_cast<panzer::ClosureModelFactoryBase*>(new ClosureModelFactory<EvalT>));
  }
};

}  // namespace flujo

#endif /** __Flujo_ClosureModel_Factory_TemplateBuilder_HPP__ */
