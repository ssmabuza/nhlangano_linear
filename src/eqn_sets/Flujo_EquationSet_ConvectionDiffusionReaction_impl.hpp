// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_EquationSet_ConvectionDiffusionReaction_impl_HPP__
#define __Flujo_EquationSet_ConvectionDiffusionReaction_impl_HPP__

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"

namespace flujo {

template <typename EvalT>
EquationSet_ConvectionDiffusionReaction<EvalT>::EquationSet_ConvectionDiffusionReaction(
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const int& default_integration_order,
    const panzer::CellData& cell_data,
    const Teuchos::RCP<panzer::GlobalData>& global_data,
    const bool build_transient_support)
    : panzer::EquationSet_DefaultImpl<EvalT>(
          params, default_integration_order, cell_data, global_data, build_transient_support),
      concentration_name_("conc") {
  Teuchos::ParameterList valid_parameters;
  this->setDefaultValidParameters(valid_parameters);
  valid_parameters.set("Model ID", "", "Closure model id associated with this equation set");
  valid_parameters.set("Basis Order", 1, "Order of the basis");
  valid_parameters.set("Integration Order", 2, "Order of the integration");
  valid_parameters.set("Basis Type", "HGrad", "Basis for concentration field");
  params->validateParametersAndSetDefaults(valid_parameters);

  const int basis_order = params->get<int>("Basis Order");
  const int integration_order = params->get<int>("Integration Order");
  const std::string basis_type = params->get<std::string>("Basis Type");
  const std::string model_id = params->get<std::string>("Model ID");

  this->addDOF(concentration_name_, basis_type, basis_order, integration_order);
  this->addDOFGrad(concentration_name_);
  this->addDOFTimeDerivative(concentration_name_);
  this->addClosureModel(model_id);
  this->setupDOFs();
}

template <typename EvalT>
void EquationSet_ConvectionDiffusionReaction<EvalT>::buildAndRegisterEquationSetEvaluators(
    PHX::FieldManager<panzer::Traits>&,
    const panzer::FieldLibrary&,
    const Teuchos::ParameterList&) const {
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Flujo EquationSet_ConvectionDiffusionReaction evaluators are not implemented yet.");
}

}  // namespace flujo

#endif /** __Flujo_EquationSet_ConvectionDiffusionReaction_impl_HPP__ */
