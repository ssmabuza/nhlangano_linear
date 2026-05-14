// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_Convection_HPP__
#define __Flujo_Convection_HPP__

#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"
#include "Teuchos_ParameterList.hpp"

namespace flujo {

template <typename EvalT, typename Traits>
class Convection : public panzer::EvaluatorWithBaseImpl<Traits>,
                   public PHX::EvaluatorDerived<EvalT, Traits> {
public:
  Convection(const Teuchos::ParameterList& params);

  void evaluateFields(typename Traits::EvalData workset);

private:
  using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<ScalarT, panzer::Cell, panzer::Point> convection_flux_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim> velocity_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim> scalar_gradient_;
  double multiplier_;
};

}  // namespace flujo

#include "Flujo_Convection_impl.hpp"

#endif /** __Flujo_Convection_HPP__ */
