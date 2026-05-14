// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_EquationSet_Biot_HPP__
#define __Flujo_EquationSet_Biot_HPP__

#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Teuchos_RCP.hpp"

namespace flujo {

template <typename EvalT>
class EquationSet_Biot : public panzer::EquationSet_DefaultImpl<EvalT> {
public:
  EquationSet_Biot(const Teuchos::RCP<Teuchos::ParameterList>& params,
                   const int& default_integration_order,
                   const panzer::CellData& cell_data,
                   const Teuchos::RCP<panzer::GlobalData>& global_data,
                   const bool build_transient_support);

  void buildAndRegisterEquationSetEvaluators(
      PHX::FieldManager<panzer::Traits>& fm,
      const panzer::FieldLibrary& field_library,
      const Teuchos::ParameterList& user_data) const override;

private:
  int dimension_;
  std::string displacement_name_;
  std::string pore_pressure_name_;
  std::string filtration_velocity_name_;
  std::string inverse_permeability_name_;
  std::string fluid_source_name_;
  std::string body_force_name_;
  std::string poroelastic_stress_flux_name_;
  std::string coupled_fluid_traction_name_;
  std::string structure_density_name_;
  std::string storage_coefficient_name_;
  std::string biot_willis_name_;
  std::string spring_coefficient_name_;
};

}  // namespace flujo

#include "Flujo_EquationSet_Biot_impl.hpp"

#endif /** __Flujo_EquationSet_Biot_HPP__ */
