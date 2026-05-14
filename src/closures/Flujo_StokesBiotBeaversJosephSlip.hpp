// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_StokesBiotBeaversJosephSlip_HPP__
#define __Flujo_StokesBiotBeaversJosephSlip_HPP__

#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

namespace flujo {

template <typename EvalT, typename Traits>
class StokesBiotBeaversJosephSlip : public panzer::EvaluatorWithBaseImpl<Traits>,
                                       public PHX::EvaluatorDerived<EvalT, Traits> {
public:
  StokesBiotBeaversJosephSlip(const std::string& slip_residual_name,
                              const std::string& fluid_velocity_name,
                              const std::string& tangential_traction_name,
                              const panzer::IntegrationRule& ir,
                              const double slip_coefficient);

  void evaluateFields(typename Traits::EvalData workset);

private:
  using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<ScalarT, panzer::Cell, panzer::Point, panzer::Dim> slip_residual_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim> fluid_velocity_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim> tangential_traction_;
  double slip_coefficient_;
  int spatial_dimension_;
};

}  // namespace flujo

#include "Flujo_StokesBiotBeaversJosephSlip_impl.hpp"

#endif /** __Flujo_StokesBiotBeaversJosephSlip_HPP__ */
