// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_StokesBiotNormalFiltration_impl_HPP__
#define __Flujo_StokesBiotNormalFiltration_impl_HPP__

namespace flujo {

template <typename EvalT, typename Traits>
StokesBiotNormalFiltration<EvalT, Traits>::StokesBiotNormalFiltration(
    const std::string& filtration_jump_name,
    const std::string& fluid_velocity_name,
    const std::string& structure_velocity_name,
    const panzer::IntegrationRule& ir,
    const std::vector<double>& interface_normal)
    : spatial_dimension_(ir.spatial_dimension) {
  filtration_jump_ = PHX::MDField<ScalarT, panzer::Cell, panzer::Point>(
      filtration_jump_name, ir.dl_scalar);
  fluid_velocity_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim>(
      fluid_velocity_name, ir.dl_vector);
  structure_velocity_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim>(
      structure_velocity_name, ir.dl_vector);

  interface_normal_ = Kokkos::View<double*, Kokkos::HostSpace>("interface_normal", spatial_dimension_);
  for (int i = 0; i < spatial_dimension_; ++i) {
    interface_normal_(i) = interface_normal.at(i);
  }

  this->addEvaluatedField(filtration_jump_);
  this->addDependentField(fluid_velocity_);
  this->addDependentField(structure_velocity_);
  this->setName("StokesBiotNormalFiltration");
}

template <typename EvalT, typename Traits>
void StokesBiotNormalFiltration<EvalT, Traits>::evaluateFields(typename Traits::EvalData workset) {
  using panzer::index_t;

  const auto filtration_jump = filtration_jump_;
  const auto fluid_velocity = fluid_velocity_;
  const auto structure_velocity = structure_velocity_;
  const int spatial_dimension = spatial_dimension_;
  const double nx = interface_normal_(0);
  const double ny = spatial_dimension > 1 ? interface_normal_(1) : 0.0;
  const double nz = spatial_dimension > 2 ? interface_normal_(2) : 0.0;

  Kokkos::MDRangePolicy<PHX::exec_space, Kokkos::Rank<2>> policy(
      {0, 0}, {workset.num_cells, filtration_jump.extent_int(1)});

  Kokkos::parallel_for(
      "flujo:StokesBiotNormalFiltration", policy,
      KOKKOS_LAMBDA(const index_t cell, const index_t point) {
        ScalarT fluid_normal = ScalarT(0.0);
        ScalarT structure_normal = ScalarT(0.0);
        for (int i = 0; i < spatial_dimension; ++i) {
          const ScalarT normal_component = (i == 0 ? nx : (i == 1 ? ny : nz));
          fluid_normal += fluid_velocity(cell, point, i) * normal_component;
          structure_normal += structure_velocity(cell, point, i) * normal_component;
        }
        filtration_jump(cell, point) = fluid_normal - structure_normal;
      });
}

}  // namespace flujo

#endif /** __Flujo_StokesBiotNormalFiltration_impl_HPP__ */
