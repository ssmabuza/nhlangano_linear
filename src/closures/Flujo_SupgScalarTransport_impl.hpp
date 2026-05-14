// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_SupgScalarTransport_impl_HPP__
#define __Flujo_SupgScalarTransport_impl_HPP__

#include <cmath>

namespace flujo {

template <typename EvalT, typename Traits>
SupgScalarTransport<EvalT, Traits>::SupgScalarTransport(
    const std::string& stabilization_name,
    const std::string& velocity_name,
    const std::string& scalar_gradient_name,
    const std::string& source_name,
    const panzer::IntegrationRule& ir,
    const double tau_scale,
    const double reference_length)
    : tau_scale_(tau_scale), reference_length_(reference_length) {
  stabilization_ = PHX::MDField<ScalarT, panzer::Cell, panzer::Point>(stabilization_name, ir.dl_scalar);
  velocity_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim>(
      velocity_name, ir.dl_vector);
  scalar_gradient_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim>(
      scalar_gradient_name, ir.dl_vector);
  source_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point>(source_name, ir.dl_scalar);

  this->addEvaluatedField(stabilization_);
  this->addDependentField(velocity_);
  this->addDependentField(scalar_gradient_);
  this->addDependentField(source_);
  this->setName("SupgScalarTransport");
}

template <typename EvalT, typename Traits>
void SupgScalarTransport<EvalT, Traits>::evaluateFields(typename Traits::EvalData workset) {
  using panzer::index_t;

  const auto stabilization = stabilization_;
  const auto velocity = velocity_;
  const auto scalar_gradient = scalar_gradient_;
  const auto source = source_;
  const ScalarT tau_scale = tau_scale_;
  const ScalarT reference_length = reference_length_;

  Kokkos::MDRangePolicy<PHX::exec_space, Kokkos::Rank<2>> policy(
      {0, 0}, {workset.num_cells, stabilization.extent_int(1)});

  Kokkos::parallel_for(
      "flujo:SupgScalarTransport", policy,
      KOKKOS_LAMBDA(const index_t cell, const index_t point) {
        ScalarT velocity_magnitude = ScalarT(0.0);
        ScalarT advection_residual = ScalarT(0.0);
        for (int dim = 0; dim < static_cast<int>(velocity.extent_int(2)); ++dim) {
          const ScalarT velocity_component = velocity(cell, point, dim);
          velocity_magnitude += velocity_component * velocity_component;
          advection_residual += velocity_component * scalar_gradient(cell, point, dim);
        }
        velocity_magnitude = Kokkos::sqrt(velocity_magnitude);
        const ScalarT tau =
            tau_scale * reference_length / (ScalarT(2.0) * velocity_magnitude + ScalarT(1.0e-12));
        stabilization(cell, point) = tau * (advection_residual - source(cell, point));
      });
}

}  // namespace flujo

#endif /** __Flujo_SupgScalarTransport_impl_HPP__ */
