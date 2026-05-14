// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_ArrheniusReactionSource_impl_HPP__
#define __Flujo_ArrheniusReactionSource_impl_HPP__

#include <cmath>

namespace flujo {

template <typename EvalT, typename Traits>
ArrheniusReactionSource<EvalT, Traits>::ArrheniusReactionSource(
    const std::string& source_name,
    const std::string& temperature_name,
    const std::string& reactant_concentration_name,
    const panzer::IntegrationRule& ir,
    const double pre_exponential_factor,
    const double activation_energy,
    const double gas_constant,
    const double reactant_order)
    : pre_exponential_factor_(pre_exponential_factor),
      activation_energy_(activation_energy),
      gas_constant_(gas_constant),
      reactant_order_(reactant_order) {
  source_ = PHX::MDField<ScalarT, panzer::Cell, panzer::Point>(source_name, ir.dl_scalar);
  temperature_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point>(temperature_name, ir.dl_scalar);
  reactant_concentration_ =
      PHX::MDField<const ScalarT, panzer::Cell, panzer::Point>(reactant_concentration_name, ir.dl_scalar);

  this->addEvaluatedField(source_);
  this->addDependentField(temperature_);
  this->addDependentField(reactant_concentration_);
  this->setName("ArrheniusReactionSource");
}

template <typename EvalT, typename Traits>
void ArrheniusReactionSource<EvalT, Traits>::evaluateFields(typename Traits::EvalData workset) {
  using panzer::index_t;

  const auto source = source_;
  const auto temperature = temperature_;
  const auto reactant_concentration = reactant_concentration_;
  const ScalarT pre_exponential_factor = pre_exponential_factor_;
  const ScalarT activation_energy = activation_energy_;
  const ScalarT gas_constant = gas_constant_;
  const ScalarT reactant_order = reactant_order_;

  Kokkos::MDRangePolicy<PHX::exec_space, Kokkos::Rank<2>> policy(
      {0, 0}, {workset.num_cells, source.extent_int(1)});

  Kokkos::parallel_for(
      "flujo:ArrheniusReactionSource", policy,
      KOKKOS_LAMBDA(const index_t cell, const index_t point) {
        const ScalarT temperature_value = temperature(cell, point);
        const ScalarT concentration_value = reactant_concentration(cell, point);
        const ScalarT exponent = -activation_energy / (gas_constant * temperature_value);
        const ScalarT rate = pre_exponential_factor * Kokkos::exp(exponent);
        source(cell, point) = rate * Kokkos::pow(concentration_value, reactant_order);
      });
}

}  // namespace flujo

#endif /** __Flujo_ArrheniusReactionSource_impl_HPP__ */
