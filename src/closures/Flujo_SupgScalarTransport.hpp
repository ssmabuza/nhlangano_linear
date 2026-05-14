// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_SupgScalarTransport_HPP__
#define __Flujo_SupgScalarTransport_HPP__

#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

namespace flujo {

template <typename EvalT, typename Traits>
class SupgScalarTransport : public panzer::EvaluatorWithBaseImpl<Traits>,
                            public PHX::EvaluatorDerived<EvalT, Traits> {
public:
  SupgScalarTransport(const std::string& stabilization_name,
                      const std::string& velocity_name,
                      const std::string& scalar_gradient_name,
                      const std::string& source_name,
                      const panzer::IntegrationRule& ir,
                      const double tau_scale,
                      const double reference_length);

  void evaluateFields(typename Traits::EvalData workset);

private:
  using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<ScalarT, panzer::Cell, panzer::Point> stabilization_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim> velocity_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim> scalar_gradient_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point> source_;
  double tau_scale_;
  double reference_length_;
};

}  // namespace flujo

#include "Flujo_SupgScalarTransport_impl.hpp"

#endif /** __Flujo_SupgScalarTransport_HPP__ */
