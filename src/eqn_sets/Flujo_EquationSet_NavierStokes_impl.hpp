// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_EquationSet_NavierStokes_impl_HPP__
#define __Flujo_EquationSet_NavierStokes_impl_HPP__

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"

namespace flujo {

template <typename EvalT>
EquationSet_NavierStokes<EvalT>::EquationSet_NavierStokes(
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const int& default_integration_order,
    const panzer::CellData& cell_data,
    const Teuchos::RCP<panzer::GlobalData>& global_data,
    const bool build_transient_support)
    : panzer::EquationSet_DefaultImpl<EvalT>(
          params, default_integration_order, cell_data, global_data, build_transient_support),
      dimension_(cell_data.baseCellDimension()),
      velocity_name_("velocity"),
      pressure_name_("pressure") {
  Teuchos::ParameterList valid_parameters;
  this->setDefaultValidParameters(valid_parameters);
  valid_parameters.set("Model ID", "", "Closure model id associated with this equation set");
  valid_parameters.set("Basis Order", 1, "Order of the basis");
  valid_parameters.set("Integration Order", 2, "Order of the integration");
  valid_parameters.set("Velocity Basis Type", "HGrad", "Basis for velocity field");
  valid_parameters.set("Pressure Basis Type", "HGrad", "Basis for pressure field");
  params->validateParametersAndSetDefaults(valid_parameters);

  const int basis_order = params->get<int>("Basis Order");
  const int integration_order = params->get<int>("Integration Order");
  const std::string model_id = params->get<std::string>("Model ID");
  const std::string velocity_basis = params->get<std::string>("Velocity Basis Type");
  const std::string pressure_basis = params->get<std::string>("Pressure Basis Type");

  TEUCHOS_TEST_FOR_EXCEPTION(
      dimension_ < 2 || dimension_ > 3, std::logic_error,
      "EquationSet_NavierStokes only supports 2D/3D cell topologies.");

  this->addDOF(velocity_name_, velocity_basis, basis_order, integration_order);
  this->addDOFGrad(velocity_name_);
  this->addDOFTimeDerivative(velocity_name_);

  this->addDOF(pressure_name_, pressure_basis, basis_order, integration_order);
  this->addDOFGrad(pressure_name_);

  this->addClosureModel(model_id);
  this->setupDOFs();
}

template <typename EvalT>
void EquationSet_NavierStokes<EvalT>::buildAndRegisterEquationSetEvaluators(
    PHX::FieldManager<panzer::Traits>&,
    const panzer::FieldLibrary&,
    const Teuchos::ParameterList&) const {
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Flujo EquationSet_NavierStokes evaluators are not implemented yet.");
}

}  // namespace flujo

#endif /** __Flujo_EquationSet_NavierStokes_impl_HPP__ */
