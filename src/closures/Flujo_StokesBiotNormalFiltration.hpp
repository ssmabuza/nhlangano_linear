// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_StokesBiotNormalFiltration_HPP__
#define __Flujo_StokesBiotNormalFiltration_HPP__

#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

namespace flujo {

template <typename EvalT, typename Traits>
class StokesBiotNormalFiltration : public panzer::EvaluatorWithBaseImpl<Traits>,
                                     public PHX::EvaluatorDerived<EvalT, Traits> {
public:
  StokesBiotNormalFiltration(const std::string& filtration_jump_name,
                             const std::string& fluid_velocity_name,
                             const std::string& structure_velocity_name,
                             const panzer::IntegrationRule& ir,
                             const std::vector<double>& interface_normal);

  void evaluateFields(typename Traits::EvalData workset);

private:
  using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<ScalarT, panzer::Cell, panzer::Point> filtration_jump_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim> fluid_velocity_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim> structure_velocity_;
  Kokkos::View<double*, Kokkos::HostSpace> interface_normal_;
  int spatial_dimension_;
};

}  // namespace flujo

#include "Flujo_StokesBiotNormalFiltration_impl.hpp"

#endif /** __Flujo_StokesBiotNormalFiltration_HPP__ */
