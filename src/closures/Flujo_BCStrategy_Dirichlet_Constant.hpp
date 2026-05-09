// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_BCStrategy_Dirichlet_Constant_HPP__
#define __Flujo_BCStrategy_Dirichlet_Constant_HPP__

#include "Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Teuchos_RCP.hpp"

namespace flujo {

template <typename EvalT>
class BCStrategy_Dirichlet_Constant
    : public panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT> {
public:
  BCStrategy_Dirichlet_Constant(const panzer::BC& bc,
                                const Teuchos::RCP<panzer::GlobalData>& global_data);

  void setup(const panzer::PhysicsBlock& side_pb,
             const Teuchos::ParameterList& user_data) override;

  void buildAndRegisterEvaluators(
      PHX::FieldManager<panzer::Traits>& fm, const panzer::PhysicsBlock& pb,
      const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
      const Teuchos::ParameterList& models,
      const Teuchos::ParameterList& user_data) const override;

private:
  Teuchos::RCP<panzer::PureBasis> basis_;
};

}  // namespace flujo

#include "Flujo_BCStrategy_Dirichlet_Constant_impl.hpp"

#endif /** __Flujo_BCStrategy_Dirichlet_Constant_HPP__ */
