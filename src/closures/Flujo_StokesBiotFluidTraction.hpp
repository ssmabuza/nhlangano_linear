// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_StokesBiotFluidTraction_HPP__
#define __Flujo_StokesBiotFluidTraction_HPP__

#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

namespace flujo {

template <typename EvalT, typename Traits>
class StokesBiotFluidTraction : public panzer::EvaluatorWithBaseImpl<Traits>,
                                  public PHX::EvaluatorDerived<EvalT, Traits> {
public:
  StokesBiotFluidTraction(const std::string& normal_traction_name,
                          const std::string& tangential_traction_name,
                          const std::string& fluid_pressure_name,
                          const std::string& fluid_velocity_grad_name,
                          const panzer::IntegrationRule& ir,
                          const double fluid_viscosity,
                          const std::vector<double>& interface_normal);

  void evaluateFields(typename Traits::EvalData workset);

private:
  using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<ScalarT, panzer::Cell, panzer::Point> normal_traction_;
  PHX::MDField<ScalarT, panzer::Cell, panzer::Point, panzer::Dim> tangential_traction_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point> fluid_pressure_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim, panzer::Dim> fluid_velocity_grad_;
  double fluid_viscosity_;
  Kokkos::View<double*, Kokkos::HostSpace> interface_normal_;
  int spatial_dimension_;
};

}  // namespace flujo

#include "Flujo_StokesBiotFluidTraction_impl.hpp"

#endif /** __Flujo_StokesBiotFluidTraction_HPP__ */
