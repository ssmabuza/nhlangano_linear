// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_BiotPoroelasticStress_HPP__
#define __Flujo_BiotPoroelasticStress_HPP__

#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

namespace flujo {

template <typename EvalT, typename Traits>
class BiotPoroelasticStress : public panzer::EvaluatorWithBaseImpl<Traits>,
                                public PHX::EvaluatorDerived<EvalT, Traits> {
public:
  BiotPoroelasticStress(const std::string& stress_flux_name,
                        const std::string& displacement_grad_name,
                        const std::string& pore_pressure_name,
                        const panzer::IntegrationRule& ir,
                        const double lame_mu,
                        const double lame_lambda,
                        const double biot_willis);

  void evaluateFields(typename Traits::EvalData workset);

private:
  using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<ScalarT, panzer::Cell, panzer::Point, panzer::Dim> stress_flux_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim, panzer::Dim> displacement_grad_;
  PHX::MDField<const ScalarT, panzer::Cell, panzer::Point> pore_pressure_;
  double lame_mu_;
  double lame_lambda_;
  double biot_willis_;
  int spatial_dimension_;
};

}  // namespace flujo

#include "Flujo_BiotPoroelasticStress_impl.hpp"

#endif /** __Flujo_BiotPoroelasticStress_HPP__ */
