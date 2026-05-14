// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_EquationSet_Biot_impl_HPP__
#define __Flujo_EquationSet_Biot_impl_HPP__

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Integrator_BasisTimesVector.hpp"
#include "Panzer_Integrator_DivBasisTimesScalar.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"

namespace flujo {

template <typename EvalT>
EquationSet_Biot<EvalT>::EquationSet_Biot(
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const int& default_integration_order,
    const panzer::CellData& cell_data,
    const Teuchos::RCP<panzer::GlobalData>& global_data,
    const bool build_transient_support)
    : panzer::EquationSet_DefaultImpl<EvalT>(
          params, default_integration_order, cell_data, global_data, build_transient_support),
      dimension_(cell_data.baseCellDimension()),
      displacement_name_("displacement"),
      pore_pressure_name_("pore_pressure"),
      filtration_velocity_name_("filtration_velocity"),
      inverse_permeability_name_("inverse_permeability"),
      fluid_source_name_("fluid_source"),
      body_force_name_("body_force"),
      poroelastic_stress_flux_name_("poroelastic_stress_flux"),
      coupled_fluid_traction_name_("coupled_fluid_traction"),
      structure_density_name_("structure_density"),
      storage_coefficient_name_("storage_coefficient"),
      biot_willis_name_("biot_willis"),
      spring_coefficient_name_("spring_coefficient") {
  Teuchos::ParameterList valid_parameters;
  this->setDefaultValidParameters(valid_parameters);
  valid_parameters.set("Model ID", "", "Closure model id associated with this equation set");
  valid_parameters.set("Basis Order", 1, "Order of the basis");
  valid_parameters.set("Integration Order", 2, "Order of the integration");
  valid_parameters.set("Displacement Basis Type", "HGrad", "Basis for solid displacement");
  valid_parameters.set("Pressure Basis Type", "HGrad", "Basis for pore pressure");
  valid_parameters.set("Filtration Basis Type", "HDiv", "Basis for Darcy filtration velocity");
  valid_parameters.set("Inverse Permeability", inverse_permeability_name_,
                       "Closure field for kappa^{-1}");
  valid_parameters.set("Fluid Source", fluid_source_name_, "Closure field for G_b");
  valid_parameters.set("Body Force", body_force_name_, "Closure field for F_b");
  valid_parameters.set("Poroelastic Stress Flux", poroelastic_stress_flux_name_,
                       "Closure field for the poroelastic stress flux");
  valid_parameters.set("Coupled Fluid Traction", coupled_fluid_traction_name_,
                       "Closure field for Stokes-Biot interface traction");
  valid_parameters.set("Structure Density", structure_density_name_,
                       "Closure field for rho_b");
  valid_parameters.set("Storage Coefficient", storage_coefficient_name_,
                       "Closure field for c_0");
  valid_parameters.set("Biot-Willis Parameter", biot_willis_name_, "Closure field for alpha");
  valid_parameters.set("Spring Coefficient", spring_coefficient_name_, "Closure field for gamma");
  params->validateParametersAndSetDefaults(valid_parameters);

  const int basis_order = params->get<int>("Basis Order");
  const int integration_order = params->get<int>("Integration Order");
  const std::string model_id = params->get<std::string>("Model ID");
  const std::string displacement_basis = params->get<std::string>("Displacement Basis Type");
  const std::string pressure_basis = params->get<std::string>("Pressure Basis Type");
  const std::string filtration_basis = params->get<std::string>("Filtration Basis Type");

  inverse_permeability_name_ = params->get<std::string>("Inverse Permeability");
  fluid_source_name_ = params->get<std::string>("Fluid Source");
  body_force_name_ = params->get<std::string>("Body Force");
  poroelastic_stress_flux_name_ = params->get<std::string>("Poroelastic Stress Flux");
  coupled_fluid_traction_name_ = params->get<std::string>("Coupled Fluid Traction");
  structure_density_name_ = params->get<std::string>("Structure Density");
  storage_coefficient_name_ = params->get<std::string>("Storage Coefficient");
  biot_willis_name_ = params->get<std::string>("Biot-Willis Parameter");
  spring_coefficient_name_ = params->get<std::string>("Spring Coefficient");

  TEUCHOS_TEST_FOR_EXCEPTION(
      dimension_ < 2 || dimension_ > 3, std::logic_error,
      "EquationSet_Biot only supports 2D/3D cell topologies.");

  this->addDOF(displacement_name_, displacement_basis, basis_order, integration_order);
  this->addDOFGrad(displacement_name_);
  this->addDOFTimeDerivative(displacement_name_);

  this->addDOF(pore_pressure_name_, pressure_basis, basis_order, integration_order);
  this->addDOFGrad(pore_pressure_name_);
  this->addDOFTimeDerivative(pore_pressure_name_);

  this->addDOF(filtration_velocity_name_, filtration_basis, basis_order, integration_order);
  this->addDOFDiv(filtration_velocity_name_);

  this->addClosureModel(model_id);
  this->setupDOFs();
}

template <typename EvalT>
void EquationSet_Biot<EvalT>::buildAndRegisterEquationSetEvaluators(
    PHX::FieldManager<panzer::Traits>& fm,
    const panzer::FieldLibrary&,
    const Teuchos::ParameterList&) const {
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;
  using std::vector;

  {
    RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(displacement_name_);
    RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(displacement_name_);

    vector<string> residual_operator_names;
    if (this->buildTransientSupport()) {
      const string resid = "RESIDUAL_" + displacement_name_ + "_TIME_OP";
      ParameterList p("Time Derivative " + displacement_name_);
      p.set("Residual Name", resid);
      p.set("Value Name", "DXDT_" + displacement_name_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      const vector<string> field_multiplier_names{structure_density_name_};
      p.set("Field Multipliers", Teuchos::rcpFromRef(field_multiplier_names));
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_BasisTimesVector<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    }

    {
      const string resid = "RESIDUAL_" + displacement_name_ + "_PRESSURE_GRAD";
      ParameterList p("Pore Pressure Gradient " + displacement_name_);
      p.set("Residual Name", resid);
      p.set("Value Name", "GRAD_" + pore_pressure_name_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      const vector<string> field_multiplier_names{biot_willis_name_};
      p.set("Field Multipliers", Teuchos::rcpFromRef(field_multiplier_names));
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_BasisTimesVector<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    }

    {
      const string resid = "RESIDUAL_" + displacement_name_ + "_STRESS_OP";
      ParameterList p("Poroelastic Stress " + displacement_name_);
      p.set("Residual Name", resid);
      p.set("Value Name", poroelastic_stress_flux_name_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", -1.0);
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_BasisTimesVector<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    }

    {
      const string resid = "RESIDUAL_" + displacement_name_ + "_SPRING_OP";
      ParameterList p("Spring " + displacement_name_);
      p.set("Residual Name", resid);
      p.set("Value Name", displacement_name_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      const vector<string> field_multiplier_names{spring_coefficient_name_};
      p.set("Field Multipliers", Teuchos::rcpFromRef(field_multiplier_names));
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_BasisTimesVector<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    }

    {
      const string resid = "RESIDUAL_" + displacement_name_ + "_BODY_FORCE";
      ParameterList p("Body Force " + displacement_name_);
      p.set("Residual Name", resid);
      p.set("Value Name", body_force_name_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", -1.0);
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_BasisTimesVector<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    }

    {
      const string resid = "RESIDUAL_" + displacement_name_ + "_COUPLED_TRACTION";
      ParameterList p("Coupled Fluid Traction " + displacement_name_);
      p.set("Residual Name", resid);
      p.set("Value Name", coupled_fluid_traction_name_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", -1.0);
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_BasisTimesVector<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    }

    this->buildAndRegisterResidualSummationEvaluator(fm, displacement_name_, residual_operator_names);
  }

  {
    RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(pore_pressure_name_);
    RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(pore_pressure_name_);

    vector<string> residual_operator_names;
    if (this->buildTransientSupport()) {
      const string resid = "RESIDUAL_" + pore_pressure_name_ + "_TIME_OP";
      ParameterList p("Time Derivative " + pore_pressure_name_);
      p.set("Residual Name", resid);
      p.set("Value Name", "DXDT_" + pore_pressure_name_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      const vector<string> field_multiplier_names{storage_coefficient_name_};
      p.set("Field Multipliers", Teuchos::rcpFromRef(field_multiplier_names));
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_BasisTimesScalar<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    }

    {
      const string resid = "RESIDUAL_" + pore_pressure_name_ + "_SOLID_DILATATION";
      ParameterList p("Solid Dilatation " + pore_pressure_name_);
      p.set("Residual Name", resid);
      p.set("Flux Name", "DXDT_" + displacement_name_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", -1.0);
      const vector<string> field_multiplier_names{biot_willis_name_};
      p.set("Field Multipliers", Teuchos::rcpFromRef(field_multiplier_names));
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_GradBasisDotVector<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    }

    {
      const string resid = "RESIDUAL_" + pore_pressure_name_ + "_FILTRATION_DIV";
      ParameterList p("Filtration Divergence " + pore_pressure_name_);
      p.set("Residual Name", resid);
      p.set("Value Name", "DIV_" + filtration_velocity_name_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_BasisTimesScalar<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    }

    {
      const string resid = "RESIDUAL_" + pore_pressure_name_ + "_FLUID_SOURCE";
      ParameterList p("Fluid Source " + pore_pressure_name_);
      p.set("Residual Name", resid);
      p.set("Value Name", fluid_source_name_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", -1.0);
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_BasisTimesScalar<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    }

    this->buildAndRegisterResidualSummationEvaluator(fm, pore_pressure_name_, residual_operator_names);
  }

  {
    RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(filtration_velocity_name_);
    RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(filtration_velocity_name_);

    vector<string> residual_operator_names;
    {
      const string resid = "RESIDUAL_" + filtration_velocity_name_ + "_GRAD_PRESSURE";
      ParameterList p("Pressure Gradient " + filtration_velocity_name_);
      p.set("Residual Name", resid);
      p.set("Value Name", pore_pressure_name_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_DivBasisTimesScalar<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    }

    {
      const string resid = "RESIDUAL_" + filtration_velocity_name_ + "_DARCY_MASS";
      ParameterList p("Darcy Mass " + filtration_velocity_name_);
      p.set("Residual Name", resid);
      p.set("Value Name", filtration_velocity_name_);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      const vector<string> field_multiplier_names{inverse_permeability_name_};
      p.set("Field Multipliers", Teuchos::rcpFromRef(field_multiplier_names));
      this->template registerEvaluator<EvalT>(
          fm, rcp(new panzer::Integrator_BasisTimesVector<EvalT, panzer::Traits>(p)));
      residual_operator_names.push_back(resid);
    }

    this->buildAndRegisterResidualSummationEvaluator(fm, filtration_velocity_name_, residual_operator_names);
  }
}

}  // namespace flujo

#endif /** __Flujo_EquationSet_Biot_impl_HPP__ */
