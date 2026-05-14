// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_ClosureModel_Factory_impl_HPP__
#define __Flujo_ClosureModel_Factory_impl_HPP__

#include <cmath>
#include <string>
#include <vector>

#include "Flujo_ArrheniusReactionSource.hpp"
#include "Flujo_BiotPoroelasticStress.hpp"
#include "Flujo_StokesBiotBeaversJosephSlip.hpp"
#include "Flujo_StokesBiotFluidTraction.hpp"
#include "Flujo_StokesBiotNormalFiltration.hpp"
#include "Flujo_SupgScalarTransport.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Constant.hpp"
#include "Panzer_ConstantVector.hpp"
#include "Teuchos_TestForException.hpp"

namespace flujo {

namespace {

std::vector<double> readInterfaceNormal(const Teuchos::ParameterList& plist, const int spatial_dimension) {
  std::vector<double> interface_normal(spatial_dimension, 0.0);
  if (plist.isType<std::string>("Interface Normal")) {
    const std::string normal_string = plist.get<std::string>("Interface Normal");
    if (normal_string == "X") {
      interface_normal[0] = 1.0;
    } else if (normal_string == "Y") {
      interface_normal[1] = 1.0;
    } else if (normal_string == "Z") {
      interface_normal[2] = 1.0;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                 "Unsupported interface normal \"" << normal_string << "\".");
    }
  } else {
    for (int i = 0; i < spatial_dimension; ++i) {
      const std::string component_name = "Interface Normal " + std::to_string(i);
      if (plist.isType<double>(component_name)) {
        interface_normal[i] = plist.get<double>(component_name);
      }
    }
  }

  double normal_magnitude = 0.0;
  for (const double component : interface_normal) {
    normal_magnitude += component * component;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
      normal_magnitude <= 0.0, std::logic_error,
      "Stokes-Biot coupling closures require a nonzero interface normal.");
  normal_magnitude = std::sqrt(normal_magnitude);
  for (double& component : interface_normal) {
    component /= normal_magnitude;
  }
  return interface_normal;
}

}  // namespace

template <typename EvalT>
Teuchos::RCP<std::vector<Teuchos::RCP<PHX::Evaluator<panzer::Traits>>>>
ClosureModelFactory<EvalT>::buildClosureModels(
    const std::string& model_id,
    const Teuchos::ParameterList& models,
    const panzer::FieldLayoutLibrary& fl,
    const Teuchos::RCP<panzer::IntegrationRule>& ir,
    const Teuchos::ParameterList&,
    const Teuchos::ParameterList&,
    const Teuchos::RCP<panzer::GlobalData>&,
    PHX::FieldManager<panzer::Traits>&) const {
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<std::vector<RCP<PHX::Evaluator<panzer::Traits>>>> evaluators =
      rcp(new std::vector<RCP<PHX::Evaluator<panzer::Traits>>>);

  TEUCHOS_TEST_FOR_EXCEPTION(
      !models.isSublist(model_id), std::logic_error,
      "Failed to find requested closure model \"" << model_id << "\".");

  std::vector<Teuchos::RCP<const panzer::PureBasis>> bases;
  fl.uniqueBases(bases);

  const ParameterList& model_entries = models.sublist(model_id);
  for (ParameterList::ConstIterator model_it = model_entries.begin();
       model_it != model_entries.end(); ++model_it) {
    const std::string key = model_it->first;
    const ParameterList& plist = Teuchos::getValue<ParameterList>(model_it->second);
    bool found = false;

    if (plist.isType<double>("Value")) {
      {
        ParameterList input;
        input.set("Name", key);
        input.set("Value", plist.get<double>("Value"));
        input.set("Data Layout", ir->dl_scalar);
        evaluators->push_back(rcp(new panzer::Constant<EvalT, panzer::Traits>(input)));
      }

      for (const auto& basis : bases) {
        ParameterList input;
        input.set("Name", key);
        input.set("Value", plist.get<double>("Value"));
        const Teuchos::RCP<const panzer::BasisIRLayout> basis_layout = basisIRLayout(*basis, *ir);
        input.set("Data Layout", basis_layout->functional);
        evaluators->push_back(rcp(new panzer::Constant<EvalT, panzer::Traits>(input)));
      }
      found = true;
    }

    if (plist.isType<std::string>("Type")) {
      const std::string type = plist.get<std::string>("Type");

      if (type == "BIOT POROELASTIC STRESS") {
        const double lame_mu = plist.get<double>("Lame Mu");
        const double lame_lambda = plist.get<double>("Lame Lambda");
        const double biot_willis = plist.get<double>("Biot-Willis Parameter");
        const std::string displacement_grad_name = plist.get<std::string>("Displacement Gradient Name");
        const std::string pore_pressure_name = plist.get<std::string>("Pore Pressure Name");
        evaluators->push_back(rcp(new BiotPoroelasticStress<EvalT, panzer::Traits>(
            key, displacement_grad_name, pore_pressure_name, *ir, lame_mu, lame_lambda,
            biot_willis)));
        found = true;
      }

      if (type == "STOKES BIOT FLUID TRACTION") {
        const double fluid_viscosity = plist.get<double>("Fluid Viscosity");
        const std::string fluid_pressure_name = plist.get<std::string>("Fluid Pressure Name");
        const std::string fluid_velocity_grad_name = plist.get<std::string>("Fluid Velocity Gradient Name");
        const std::string normal_traction_name = plist.get<std::string>("Normal Traction Name");
        const std::string tangential_traction_name = plist.get<std::string>("Tangential Traction Name");
        const std::vector<double> interface_normal = readInterfaceNormal(plist, ir->spatial_dimension);
        evaluators->push_back(rcp(new StokesBiotFluidTraction<EvalT, panzer::Traits>(
            normal_traction_name, tangential_traction_name, fluid_pressure_name,
            fluid_velocity_grad_name, *ir, fluid_viscosity, interface_normal)));
        found = true;
      }

      if (type == "STOKES BIOT NORMAL FILTRATION") {
        const std::string fluid_velocity_name = plist.get<std::string>("Fluid Velocity Name");
        const std::string structure_velocity_name = plist.get<std::string>("Structure Velocity Name");
        const std::vector<double> interface_normal = readInterfaceNormal(plist, ir->spatial_dimension);
        evaluators->push_back(rcp(new StokesBiotNormalFiltration<EvalT, panzer::Traits>(
            key, fluid_velocity_name, structure_velocity_name, *ir, interface_normal)));
        found = true;
      }

      if (type == "STOKES BIOT BJS SLIP") {
        const double slip_coefficient = plist.get<double>("Slip Coefficient");
        const std::string fluid_velocity_name = plist.get<std::string>("Fluid Velocity Name");
        const std::string tangential_traction_name = plist.get<std::string>("Tangential Traction Name");
        evaluators->push_back(rcp(new StokesBiotBeaversJosephSlip<EvalT, panzer::Traits>(
            key, fluid_velocity_name, tangential_traction_name, *ir, slip_coefficient)));
        found = true;
      }

      if (type == "ARRHENIUS REACTION SOURCE") {
        const double pre_exponential_factor = plist.get<double>("Pre-Exponential Factor");
        const double activation_energy = plist.get<double>("Activation Energy");
        const double gas_constant = plist.get<double>("Gas Constant");
        const double reactant_order = plist.get<double>("Reactant Order");
        const std::string temperature_name = plist.get<std::string>("Temperature Name");
        const std::string reactant_concentration_name =
            plist.get<std::string>("Reactant Concentration Name");
        evaluators->push_back(rcp(new ArrheniusReactionSource<EvalT, panzer::Traits>(
            key, temperature_name, reactant_concentration_name, *ir, pre_exponential_factor,
            activation_energy, gas_constant, reactant_order)));
        found = true;
      }

      if (type == "SUPG SCALAR TRANSPORT") {
        const double tau_scale = plist.get<double>("Tau Scale");
        const double reference_length = plist.get<double>("Reference Length");
        const std::string velocity_name = plist.get<std::string>("Velocity Name");
        const std::string scalar_gradient_name = plist.get<std::string>("Scalar Gradient Name");
        const std::string source_name = plist.get<std::string>("Source Name");
        evaluators->push_back(rcp(new SupgScalarTransport<EvalT, panzer::Traits>(
            key, velocity_name, scalar_gradient_name, source_name, *ir, tau_scale,
            reference_length)));
        found = true;
      }
    }

    if (plist.isType<double>("Value X")) {
      ParameterList input;
      input.set("Name", key);
      input.set("Value X", plist.get<double>("Value X"));
      if (ir->spatial_dimension > 1) {
        input.set("Value Y", plist.get<double>("Value Y"));
      }
      if (ir->spatial_dimension > 2) {
        input.set("Value Z", plist.get<double>("Value Z"));
      }
      input.set("Data Layout", ir->dl_vector);
      evaluators->push_back(rcp(new panzer::ConstantVector<EvalT, panzer::Traits>(input)));
      found = true;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
        !found, std::logic_error,
        "Closure model \"" << key << "\" in model \"" << model_id << "\" is not supported.");
  }

  return evaluators;
}

}  // namespace flujo

#endif /** __Flujo_ClosureModel_Factory_impl_HPP__ */
