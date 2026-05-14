// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_EquationSet_EnergyTransport_HPP__
#define __Flujo_EquationSet_EnergyTransport_HPP__

#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Teuchos_RCP.hpp"

namespace flujo {

template <typename EvalT>
class EquationSet_EnergyTransport : public panzer::EquationSet_DefaultImpl<EvalT> {
public:
  EquationSet_EnergyTransport(const Teuchos::RCP<Teuchos::ParameterList>& params,
                              const int& default_integration_order,
                              const panzer::CellData& cell_data,
                              const Teuchos::RCP<panzer::GlobalData>& global_data,
                              const bool build_transient_support);

  void buildAndRegisterEquationSetEvaluators(
      PHX::FieldManager<panzer::Traits>& fm,
      const panzer::FieldLibrary& field_library,
      const Teuchos::ParameterList& user_data) const override;

private:
  std::string prefix_;
  std::string temperature_name_;
  std::string density_name_;
  std::string heat_capacity_name_;
  std::string thermal_conductivity_name_;
  std::string velocity_name_;
  std::string heat_source_name_;
  std::string supg_stabilization_name_;
  std::string convection_mode_;
  bool convection_in_conservation_form_;
  bool supg_enabled_;
};

}  // namespace flujo

#include "Flujo_EquationSet_EnergyTransport_impl.hpp"

#endif /** __Flujo_EquationSet_EnergyTransport_HPP__ */
