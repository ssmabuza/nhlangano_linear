// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_EquationSet_EnergyTransport_impl_HPP__
#define __Flujo_EquationSet_EnergyTransport_impl_HPP__

#include "Flujo_Convection.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_TestForException.hpp"

namespace flujo {

template <typename EvalT>
EquationSet_EnergyTransport<EvalT>::EquationSet_EnergyTransport(
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const int& default_integration_order,
    const panzer::CellData& cell_data,
    const Teuchos::RCP<panzer::GlobalData>& global_data,
    const bool build_transient_support)
    : panzer::EquationSet_DefaultImpl<EvalT>(
          params, default_integration_order, cell_data, global_data, build_transient_support),
      prefix_(""),
      temperature_name_("temperature"),
      density_name_("density"),
      heat_capacity_name_("heat_capacity"),
      thermal_conductivity_name_("thermal_conductivity"),
      velocity_name_("velocity"),
      heat_source_name_("heat_source"),
      supg_stabilization_name_("energy_supg_stabilization"),
      convection_mode_("ON"),
      convection_in_conservation_form_(true),
      supg_enabled_(false) {
  Teuchos::ParameterList valid_parameters;
  this->setDefaultValidParameters(valid_parameters);
  valid_parameters.set("Model ID", "", "Closure model id associated with this equation set");
  valid_parameters.set("Prefix", "", "Prefix for multiple energy transport fields");
  valid_parameters.set("Basis Order", 1, "Order of the basis");
  valid_parameters.set("Integration Order", 2, "Order of the integration");
  valid_parameters.set("Basis Type", "HGrad", "Basis for temperature field");
  valid_parameters.set("Temperature Field Name", temperature_name_, "Primary temperature DOF name");
  valid_parameters.set("Density", density_name_, "Closure field for density");
  valid_parameters.set("Heat Capacity", heat_capacity_name_, "Closure field for heat capacity");
  valid_parameters.set("Thermal Conductivity", thermal_conductivity_name_,
                       "Closure field for thermal conductivity");
  valid_parameters.set("Velocity", velocity_name_, "Closure field for velocity");
  valid_parameters.set("Heat Source", heat_source_name_, "Closure field for energy source term");
  valid_parameters.set("SUPG Stabilization", supg_stabilization_name_,
                       "Closure field for SUPG stabilization contribution");
  valid_parameters.set("Convection", "ON",
                       "Enables or disables the advection term",
                       Teuchos::rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("ON", "OFF"))));
  valid_parameters.set("Convection Term is in Conservation Form", true,
                       "Use conservation form for the advection term");
  valid_parameters.set("SUPG", "OFF",
                       "Enables or disables residual-based SUPG stabilization",
                       Teuchos::rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("ON", "OFF"))));
  params->validateParametersAndSetDefaults(valid_parameters);

  const int basis_order = params->get<int>("Basis Order");
  const int integration_order = params->get<int>("Integration Order");
  const std::string basis_type = params->get<std::string>("Basis Type");
  const std::string model_id = params->get<std::string>("Model ID");

  prefix_ = params->get<std::string>("Prefix");
  temperature_name_ = prefix_ + params->get<std::string>("Temperature Field Name");
  density_name_ = params->get<std::string>("Density");
  heat_capacity_name_ = params->get<std::string>("Heat Capacity");
  thermal_conductivity_name_ = params->get<std::string>("Thermal Conductivity");
  velocity_name_ = params->get<std::string>("Velocity");
  heat_source_name_ = params->get<std::string>("Heat Source");
  supg_stabilization_name_ = params->get<std::string>("SUPG Stabilization");
  convection_mode_ = params->get<std::string>("Convection");
  convection_in_conservation_form_ = params->get<bool>("Convection Term is in Conservation Form");
  supg_enabled_ = params->get<std::string>("SUPG") == "ON";

  this->addDOF(temperature_name_, basis_type, basis_order, integration_order);
  this->addDOFGrad(temperature_name_);
  this->addDOFTimeDerivative(temperature_name_);
  this->addClosureModel(model_id);
  this->setupDOFs();
}

template <typename EvalT>
void EquationSet_EnergyTransport<EvalT>::buildAndRegisterEquationSetEvaluators(
    PHX::FieldManager<panzer::Traits>& fm,
    const panzer::FieldLibrary&,
    const Teuchos::ParameterList&) const {
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;
  using std::vector;

  RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(temperature_name_);
  RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(temperature_name_);

  vector<string> residual_operator_names;
  if (this->buildTransientSupport()) {
    const string resid = "RESIDUAL_" + temperature_name_ + "_TIME_OP";
    ParameterList p("Time Derivative " + temperature_name_);
    p.set("Residual Name", resid);
    p.set("Value Name", "DXDT_" + temperature_name_);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", 1.0);
    const vector<string> field_multiplier_names{density_name_, heat_capacity_name_};
    p.set("Field Multipliers", Teuchos::rcpFromRef(field_multiplier_names));
    this->template registerEvaluator<EvalT>(
        fm, rcp(new panzer::Integrator_BasisTimesScalar<EvalT, panzer::Traits>(p)));
    residual_operator_names.push_back(resid);
  }

  if (convection_mode_ == "ON") {
    if (convection_in_conservation_form_) {
      const string resid = "RESIDUAL_" + temperature_name_ + "_CONVECTION";
      ParameterList p("Convection " + temperature_name_);
      p.set("Residual Name", resid);
      p.set("Flux Name", velocity_name_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", -1.0);
      const vector<string> field_multiplier_names{density_name_, heat_capacity_name_, temperature_name_};
      p.set("Field Multipliers", Teuchos::rcpFromRef(field_multiplier_names));
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_GradBasisDotVector<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    } else {
      const string convection_flux_name = prefix_ + temperature_name_ + "_CONVECTION_OP";
      {
        ParameterList p("Convection Operator " + temperature_name_);
        p.set("IR", ir);
        p.set("Operator Name", convection_flux_name);
        p.set("Velocity Name", velocity_name_);
        p.set("Gradient Name", "GRAD_" + temperature_name_);
        p.set("Multiplier", 1.0);
        this->template registerEvaluator<EvalT>(fm, rcp(new Convection<EvalT, panzer::Traits>(p)));
      }

      const string resid = "RESIDUAL_" + temperature_name_ + "_CONVECTION";
      ParameterList p("Convection Integrator " + temperature_name_);
      p.set("Residual Name", resid);
      p.set("Value Name", convection_flux_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      const vector<string> field_multiplier_names{density_name_, heat_capacity_name_};
      p.set("Field Multipliers", Teuchos::rcpFromRef(field_multiplier_names));
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_BasisTimesScalar<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    }
  }

  {
    const string resid = "RESIDUAL_" + temperature_name_ + "_DIFFUSION";
    ParameterList p("Diffusion " + temperature_name_);
    p.set("Residual Name", resid);
    p.set("Flux Name", "GRAD_" + temperature_name_);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", -1.0);
    const vector<string> field_multiplier_names{thermal_conductivity_name_};
    p.set("Field Multipliers", Teuchos::rcpFromRef(field_multiplier_names));
    this->template registerEvaluator<EvalT>(
        fm, rcp(new panzer::Integrator_GradBasisDotVector<EvalT, panzer::Traits>(p)));
    residual_operator_names.push_back(resid);
  }

  {
    const string resid = "RESIDUAL_" + temperature_name_ + "_SOURCE";
    ParameterList p("Heat Source " + temperature_name_);
    p.set("Residual Name", resid);
    p.set("Value Name", heat_source_name_);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", -1.0);
    this->template registerEvaluator<EvalT>(
        fm, rcp(new panzer::Integrator_BasisTimesScalar<EvalT, panzer::Traits>(p)));
    residual_operator_names.push_back(resid);
  }

  if (supg_enabled_) {
    const string resid = "RESIDUAL_" + temperature_name_ + "_SUPG";
    ParameterList p("SUPG " + temperature_name_);
    p.set("Residual Name", resid);
    p.set("Value Name", supg_stabilization_name_);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", -1.0);
    this->template registerEvaluator<EvalT>(
        fm, rcp(new panzer::Integrator_BasisTimesScalar<EvalT, panzer::Traits>(p)));
    residual_operator_names.push_back(resid);
  }

  this->buildAndRegisterResidualSummationEvaluator(fm, temperature_name_, residual_operator_names);
}

}  // namespace flujo

#endif /** __Flujo_EquationSet_EnergyTransport_impl_HPP__ */
