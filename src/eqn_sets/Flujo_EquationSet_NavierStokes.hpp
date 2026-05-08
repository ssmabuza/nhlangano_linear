// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_EquationSet_NavierStokes_HPP__
#define __Flujo_EquationSet_NavierStokes_HPP__

#include "Panzer_EquationSet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Teuchos_RCP.hpp"

namespace flujo {

template <typename EvalT>
class EquationSet_NavierStokes : public panzer::EquationSet_DefaultImpl<EvalT> {
public:
  EquationSet_NavierStokes(const Teuchos::RCP<Teuchos::ParameterList>& params,
                           const int& default_integration_order,
                           const panzer::CellData& cell_data,
                           const Teuchos::RCP<panzer::GlobalData>& global_data,
                           const bool build_transient_support);

  void buildAndRegisterEquationSetEvaluators(
      PHX::FieldManager<panzer::Traits>& fm,
      const panzer::FieldLibrary& field_library,
      const Teuchos::ParameterList& user_data) const override;

private:
  int dimension_;
  std::string velocity_name_;
  std::string pressure_name_;
};

}  // namespace flujo

#include "Flujo_EquationSet_NavierStokes_impl.hpp"

#endif /** __Flujo_EquationSet_NavierStokes_HPP__ */
