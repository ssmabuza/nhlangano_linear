// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_StokesBiotFluidTraction_impl_HPP__
#define __Flujo_StokesBiotFluidTraction_impl_HPP__

#include <cmath>

namespace flujo {

template <typename EvalT, typename Traits>
StokesBiotFluidTraction<EvalT, Traits>::StokesBiotFluidTraction(
    const std::string& normal_traction_name,
    const std::string& tangential_traction_name,
    const std::string& fluid_pressure_name,
    const std::string& fluid_velocity_grad_name,
    const panzer::IntegrationRule& ir,
    const double fluid_viscosity,
    const std::vector<double>& interface_normal)
    : fluid_viscosity_(fluid_viscosity), spatial_dimension_(ir.spatial_dimension) {
  normal_traction_ = PHX::MDField<ScalarT, panzer::Cell, panzer::Point>(
      normal_traction_name, ir.dl_scalar);
  tangential_traction_ = PHX::MDField<ScalarT, panzer::Cell, panzer::Point, panzer::Dim>(
      tangential_traction_name, ir.dl_vector);
  fluid_pressure_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point>(
      fluid_pressure_name, ir.dl_scalar);
  fluid_velocity_grad_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim, panzer::Dim>(
      fluid_velocity_grad_name, ir.dl_tensor);

  interface_normal_ = Kokkos::View<double*, Kokkos::HostSpace>("interface_normal", spatial_dimension_);
  for (int i = 0; i < spatial_dimension_; ++i) {
    interface_normal_(i) = interface_normal.at(i);
  }

  this->addEvaluatedField(normal_traction_);
  this->addEvaluatedField(tangential_traction_);
  this->addDependentField(fluid_pressure_);
  this->addDependentField(fluid_velocity_grad_);
  this->setName("StokesBiotFluidTraction");
}

template <typename EvalT, typename Traits>
void StokesBiotFluidTraction<EvalT, Traits>::evaluateFields(typename Traits::EvalData workset) {
  using panzer::index_t;

  const auto normal_traction = normal_traction_;
  const auto tangential_traction = tangential_traction_;
  const auto fluid_pressure = fluid_pressure_;
  const auto fluid_velocity_grad = fluid_velocity_grad_;
  const ScalarT fluid_viscosity = fluid_viscosity_;
  const int spatial_dimension = spatial_dimension_;
  const double nx = interface_normal_(0);
  const double ny = spatial_dimension > 1 ? interface_normal_(1) : 0.0;
  const double nz = spatial_dimension > 2 ? interface_normal_(2) : 0.0;

  Kokkos::MDRangePolicy<PHX::exec_space, Kokkos::Rank<2>> policy(
      {0, 0}, {workset.num_cells, normal_traction.extent_int(1)});

  if (spatial_dimension == 3) {
    Kokkos::parallel_for(
        "flujo:StokesBiotFluidTraction3D", policy,
        KOKKOS_LAMBDA(const index_t cell, const index_t point) {

          ScalarT rate_normal = ScalarT(0.0);
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
              const ScalarT rate_ij = ScalarT(0.5) *
                                      (fluid_velocity_grad(cell, point, i, j) +
                                       fluid_velocity_grad(cell, point, j, i));
              rate_normal += nx * rate_ij * (j == 0 ? nx : (j == 1 ? ny : nz));
            }
          }

          const ScalarT traction_normal =
              -fluid_pressure(cell, point) + ScalarT(2.0) * fluid_viscosity * rate_normal;
          normal_traction(cell, point) = traction_normal;

          ScalarT traction_tangent[3] = {ScalarT(0.0), ScalarT(0.0), ScalarT(0.0)};
          for (int i = 0; i < 3; ++i) {
            ScalarT traction_i = ScalarT(0.0);
            for (int j = 0; j < 3; ++j) {
              const ScalarT rate_ij = ScalarT(0.5) *
                                      (fluid_velocity_grad(cell, point, i, j) +
                                       fluid_velocity_grad(cell, point, j, i));
              traction_i += (-fluid_pressure(cell, point) * (i == j ? ScalarT(1.0) : ScalarT(0.0)) +
                             ScalarT(2.0) * fluid_viscosity * rate_ij) *
                            (j == 0 ? nx : (j == 1 ? ny : nz));
            }
            traction_tangent[i] = traction_i - traction_normal * (i == 0 ? nx : (i == 1 ? ny : nz));
          }

          tangential_traction(cell, point, 0) = traction_tangent[0];
          tangential_traction(cell, point, 1) = traction_tangent[1];
          tangential_traction(cell, point, 2) = traction_tangent[2];
        });
  } else {
    Kokkos::parallel_for(
        "flujo:StokesBiotFluidTraction2D", policy,
        KOKKOS_LAMBDA(const index_t cell, const index_t point) {
          ScalarT rate_normal = ScalarT(0.0);
          for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
              const ScalarT rate_ij = ScalarT(0.5) *
                                      (fluid_velocity_grad(cell, point, i, j) +
                                       fluid_velocity_grad(cell, point, j, i));
              rate_normal += (i == 0 ? nx : ny) * rate_ij * (j == 0 ? nx : ny);
            }
          }

          const ScalarT traction_normal =
              -fluid_pressure(cell, point) + ScalarT(2.0) * fluid_viscosity * rate_normal;
          normal_traction(cell, point) = traction_normal;

          for (int i = 0; i < 2; ++i) {
            ScalarT traction_i = ScalarT(0.0);
            for (int j = 0; j < 2; ++j) {
              const ScalarT rate_ij = ScalarT(0.5) *
                                      (fluid_velocity_grad(cell, point, i, j) +
                                       fluid_velocity_grad(cell, point, j, i));
              traction_i += (-fluid_pressure(cell, point) * (i == j ? ScalarT(1.0) : ScalarT(0.0)) +
                             ScalarT(2.0) * fluid_viscosity * rate_ij) *
                            (j == 0 ? nx : ny);
            }
            tangential_traction(cell, point, i) = traction_i - traction_normal * (i == 0 ? nx : ny);
          }
        });
  }
}

}  // namespace flujo

#endif /** __Flujo_StokesBiotFluidTraction_impl_HPP__ */
