// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_ArrheniusReactionSource_HPP__
#define __Flujo_ArrheniusReactionSource_HPP__

#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

namespace flujo {

template <typename EvalT, typename Traits>
class ArrheniusReactionSource : public panzer::EvaluatorWithBaseImpl<Traits>,
                                public PHX::EvaluatorDerived<EvalT, Traits> {
public:
  ArrheniusReactionSource(const std::string& source_name,
                          const std::string& temperature_name,
                          const std::string& reactant_concentration_name,
                          const panzer::IntegrationRule& ir,
                          const double pre_exponential_factor,
                          const double activation_energy,
                          const double gas_constant,
                          const double reactant_order);

  void evaluateFields(typename Traits::EvalData workset);

private:
  using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<ScalarT, panzer::Cell, panzer::Point> source_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point> temperature_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point> reactant_concentration_;
  double pre_exponential_factor_;
  double activation_energy_;
  double gas_constant_;
  double reactant_order_;
};

}  // namespace flujo

#include "Flujo_ArrheniusReactionSource_impl.hpp"

#endif /** __Flujo_ArrheniusReactionSource_HPP__ */
