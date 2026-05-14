// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_BiotPoroelasticStress_impl_HPP__
#define __Flujo_BiotPoroelasticStress_impl_HPP__

namespace flujo {

template <typename EvalT, typename Traits>
BiotPoroelasticStress<EvalT, Traits>::BiotPoroelasticStress(
    const std::string& stress_flux_name,
    const std::string& displacement_grad_name,
    const std::string& pore_pressure_name,
    const panzer::IntegrationRule& ir,
    const double lame_mu,
    const double lame_lambda,
    const double biot_willis)
    : lame_mu_(lame_mu),
      lame_lambda_(lame_lambda),
      biot_willis_(biot_willis),
      spatial_dimension_(ir.spatial_dimension) {
  stress_flux_ = PHX::MDField<ScalarT, panzer::Cell, panzer::Point, panzer::Dim>(
      stress_flux_name, ir.dl_vector);
  displacement_grad_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim, panzer::Dim>(
      displacement_grad_name, ir.dl_tensor);
  pore_pressure_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point>(
      pore_pressure_name, ir.dl_scalar);

  this->addEvaluatedField(stress_flux_);
  this->addDependentField(displacement_grad_);
  this->addDependentField(pore_pressure_);
  this->setName("BiotPoroelasticStress");
}

template <typename EvalT, typename Traits>
void BiotPoroelasticStress<EvalT, Traits>::evaluateFields(typename Traits::EvalData workset) {
  using panzer::index_t;

  const auto stress_flux = stress_flux_;
  const auto displacement_grad = displacement_grad_;
  const auto pore_pressure = pore_pressure_;
  const ScalarT lame_mu = lame_mu_;
  const ScalarT lame_lambda = lame_lambda_;
  const ScalarT biot_willis = biot_willis_;
  const int spatial_dimension = spatial_dimension_;

  Kokkos::MDRangePolicy<PHX::exec_space, Kokkos::Rank<2>> policy(
      {0, 0}, {workset.num_cells, stress_flux.extent_int(1)});

  if (spatial_dimension == 3) {
    Kokkos::parallel_for(
        "flujo:BiotPoroelasticStress3D", policy,
        KOKKOS_LAMBDA(const index_t cell, const index_t point) {
          ScalarT trace = displacement_grad(cell, point, 0, 0) +
                          displacement_grad(cell, point, 1, 1) +
                          displacement_grad(cell, point, 2, 2);
          const ScalarT pressure = pore_pressure(cell, point);

          for (int i = 0; i < 3; ++i) {
            ScalarT row_sum = ScalarT(0.0);
            for (int j = 0; j < 3; ++j) {
              const ScalarT strain_ij =
                  ScalarT(0.5) * (displacement_grad(cell, point, i, j) +
                                  displacement_grad(cell, point, j, i));
              const ScalarT stress_ij =
                  ScalarT(2.0) * lame_mu * strain_ij +
                  (j == i ? lame_lambda * trace - biot_willis * pressure : ScalarT(0.0));
              row_sum += stress_ij;
            }
            stress_flux(cell, point, i) = row_sum;
          }
        });
  } else {
    Kokkos::parallel_for(
        "flujo:BiotPoroelasticStress2D", policy,
        KOKKOS_LAMBDA(const index_t cell, const index_t point) {
          const ScalarT trace =
              displacement_grad(cell, point, 0, 0) + displacement_grad(cell, point, 1, 1);
          const ScalarT pressure = pore_pressure(cell, point);

          for (int i = 0; i < 2; ++i) {
            ScalarT row_sum = ScalarT(0.0);
            for (int j = 0; j < 2; ++j) {
              const ScalarT strain_ij =
                  ScalarT(0.5) * (displacement_grad(cell, point, i, j) +
                                  displacement_grad(cell, point, j, i));
              const ScalarT stress_ij =
                  ScalarT(2.0) * lame_mu * strain_ij +
                  (j == i ? lame_lambda * trace - biot_willis * pressure : ScalarT(0.0));
              row_sum += stress_ij;
            }
            stress_flux(cell, point, i) = row_sum;
          }
        });
  }
}

}  // namespace flujo

#endif /** __Flujo_BiotPoroelasticStress_impl_HPP__ */
