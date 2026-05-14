// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_Convection_impl_HPP__
#define __Flujo_Convection_impl_HPP__

namespace flujo {

template <typename EvalT, typename Traits>
Convection<EvalT, Traits>::Convection(const Teuchos::ParameterList& params)
    : multiplier_(params.get<double>("Multiplier")) {
  const Teuchos::RCP<panzer::IntegrationRule> ir = params.get<Teuchos::RCP<panzer::IntegrationRule>>("IR");

  convection_flux_ = PHX::MDField<ScalarT, panzer::Cell, panzer::Point>(
      params.get<std::string>("Operator Name"), ir->dl_scalar);
  velocity_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim>(
      params.get<std::string>("Velocity Name"), ir->dl_vector);
  scalar_gradient_ = PHX::MDField<const ScalarT, panzer::Cell, panzer::Point, panzer::Dim>(
      params.get<std::string>("Gradient Name"), ir->dl_vector);

  this->addEvaluatedField(convection_flux_);
  this->addDependentField(velocity_);
  this->addDependentField(scalar_gradient_);
  this->setName("Convection: " + convection_flux_.fieldTag().name());
}

template <typename EvalT, typename Traits>
void Convection<EvalT, Traits>::evaluateFields(typename Traits::EvalData workset) {
  using panzer::index_t;

  const auto convection_flux = convection_flux_;
  const auto velocity = velocity_;
  const auto scalar_gradient = scalar_gradient_;
  const ScalarT multiplier = multiplier_;

  Kokkos::MDRangePolicy<PHX::exec_space, Kokkos::Rank<2>> policy(
      {0, 0}, {workset.num_cells, convection_flux.extent_int(1)});

  Kokkos::parallel_for(
      "flujo:Convection", policy, KOKKOS_LAMBDA(const index_t cell, const index_t point) {
        ScalarT value = ScalarT(0.0);
        for (int dim = 0; dim < static_cast<int>(velocity.extent_int(2)); ++dim) {
          value += velocity(cell, point, dim) * scalar_gradient(cell, point, dim);
        }
        convection_flux(cell, point) = multiplier * value;
      });
}

}  // namespace flujo

#endif /** __Flujo_Convection_impl_HPP__ */
