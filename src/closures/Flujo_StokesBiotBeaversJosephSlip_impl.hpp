// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_StokesBiotBeaversJosephSlip_impl_HPP__
#define __Flujo_StokesBiotBeaversJosephSlip_impl_HPP__

namespace flujo {

template <typename EvalT, typename Traits>
StokesBiotBeaversJosephSlip<EvalT, Traits>::StokesBiotBeaversJosephSlip(
    const std::string& slip_residual_name,
    const std::string& fluid_velocity_name,
    const std::string& tangential_traction_name,
    const panzer::IntegrationRule& ir,
    const double slip_coefficient)
    : slip_coefficient_(slip_coefficient), spatial_dimension_(ir.spatial_dimension) {
  slip_residual_ = PHX::MDField<ScalarT, panzer::Cell, panzer::Point, panzer::Dim>(
      slip_residual_name, ir.dl_vector);
  fluid_velocity_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim>(
      fluid_velocity_name, ir.dl_vector);
  tangential_traction_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim>(
      tangential_traction_name, ir.dl_vector);

  this->addEvaluatedField(slip_residual_);
  this->addDependentField(fluid_velocity_);
  this->addDependentField(tangential_traction_);
  this->setName("StokesBiotBeaversJosephSlip");
}

template <typename EvalT, typename Traits>
void StokesBiotBeaversJosephSlip<EvalT, Traits>::evaluateFields(typename Traits::EvalData workset) {
  using panzer::index_t;

  const auto slip_residual = slip_residual_;
  const auto fluid_velocity = fluid_velocity_;
  const auto tangential_traction = tangential_traction_;
  const ScalarT slip_coefficient = slip_coefficient_;
  const int spatial_dimension = spatial_dimension_;

  Kokkos::MDRangePolicy<PHX::exec_space, Kokkos::Rank<2>> policy(
      {0, 0}, {workset.num_cells, slip_residual.extent_int(1)});

  Kokkos::parallel_for(
      "flujo:StokesBiotBeaversJosephSlip", policy,
      KOKKOS_LAMBDA(const index_t cell, const index_t point) {
        for (int i = 0; i < spatial_dimension; ++i) {
          slip_residual(cell, point, i) =
              tangential_traction(cell, point, i) + slip_coefficient * fluid_velocity(cell, point, i);
        }
      });
}

}  // namespace flujo

#endif /** __Flujo_StokesBiotBeaversJosephSlip_impl_HPP__ */
