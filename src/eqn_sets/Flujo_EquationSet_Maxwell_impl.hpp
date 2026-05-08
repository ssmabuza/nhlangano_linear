// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================


#ifndef __Flujo_EquationSet_Maxwell_impl_HPP__
#define __Flujo_EquationSet_Maxwell_impl_HPP__

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"

namespace flujo {

template <typename EvalT>
EquationSet_Maxwell<EvalT>::EquationSet_Maxwell(
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const int& default_integration_order,
    const panzer::CellData& cell_data,
    const Teuchos::RCP<panzer::GlobalData>& global_data,
    const bool build_transient_support)
    : panzer::EquationSet_DefaultImpl<EvalT>(
          params, default_integration_order, cell_data, global_data, build_transient_support),
      e_field_name_("E_edge"),
      b_field_name_("B_face") {
  Teuchos::ParameterList valid_parameters;
  this->setDefaultValidParameters(valid_parameters);
  valid_parameters.set("Model ID", "", "Closure model id associated with this equation set");
  valid_parameters.set("Basis Order", 1, "Order of the basis");
  valid_parameters.set("Integration Order", 2, "Order of the integration");
  params->validateParametersAndSetDefaults(valid_parameters);

  const int basis_order = params->get<int>("Basis Order");
  const int integration_order = params->get<int>("Integration Order");
  const std::string model_id = params->get<std::string>("Model ID");

  this->addDOF(e_field_name_, "HCurl", basis_order, integration_order);
  this->addDOFCurl(e_field_name_);
  this->addDOFTimeDerivative(e_field_name_);
  this->addDOF(b_field_name_, "HDiv", basis_order, integration_order);
  this->addDOFTimeDerivative(b_field_name_);
  this->addClosureModel(model_id);
  this->setupDOFs();
}

template <typename EvalT>
void EquationSet_Maxwell<EvalT>::buildAndRegisterEquationSetEvaluators(
    PHX::FieldManager<panzer::Traits>&,
    const panzer::FieldLibrary&,
    const Teuchos::ParameterList&) const {
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Flujo EquationSet_Maxwell evaluators are not implemented yet.");
}

}  // namespace flujo

#endif /** __Flujo_EquationSet_Maxwell_impl_HPP__ */
