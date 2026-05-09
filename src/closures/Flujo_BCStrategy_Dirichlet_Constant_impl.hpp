// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_BCStrategy_Dirichlet_Constant_impl_HPP__
#define __Flujo_BCStrategy_Dirichlet_Constant_impl_HPP__

#include "Panzer_Constant.hpp"
#include "Panzer_ConstantVector.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"

namespace flujo {

template <typename EvalT>
BCStrategy_Dirichlet_Constant<EvalT>::BCStrategy_Dirichlet_Constant(
    const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data)
    : panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc, global_data) {
  TEUCHOS_ASSERT(this->m_bc.strategy() == "Constant");
}

template <typename EvalT>
void BCStrategy_Dirichlet_Constant<EvalT>::setup(
    const panzer::PhysicsBlock& side_pb, const Teuchos::ParameterList&) {
  this->addDOF(this->m_bc.equationSetName());
  this->addTarget("Constant_" + this->m_bc.equationSetName(), this->m_bc.equationSetName(),
                  "Residual_" + this->m_bc.identifier());

  const auto& dofs = side_pb.getProvidedDOFs();
  for (const auto& dof : dofs) {
    if (dof.first == this->m_bc.equationSetName()) {
      basis_ = dof.second;
      break;
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
      Teuchos::is_null(basis_), std::runtime_error,
      "Error the name \"" << this->m_bc.equationSetName()
                          << "\" is not a valid DOF for the boundary condition:\n"
                          << this->m_bc << "\n");
}

template <typename EvalT>
void BCStrategy_Dirichlet_Constant<EvalT>::buildAndRegisterEvaluators(
    PHX::FieldManager<panzer::Traits>& fm, const panzer::PhysicsBlock&,
    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>&,
    const Teuchos::ParameterList&, const Teuchos::ParameterList&) const {
  using Teuchos::ParameterList;

  if (basis_->isScalarBasis()) {
    ParameterList p("BC Constant Dirichlet");
    p.set("Name", "Constant_" + this->m_bc.equationSetName());
    p.set("Data Layout", basis_->functional);
    p.set("Value", this->m_bc.params()->template get<double>("Value"));
    this->template registerEvaluator(fm, Teuchos::rcp(new panzer::Constant<EvalT>(p)));
  } else if (basis_->isVectorBasis()) {
    ParameterList p("BC Constant Vector Dirichlet");
    p.set("Name", "Constant_" + this->m_bc.equationSetName());
    p.set("Data Layout", basis_->functional_grad);
    p.set("Value X", this->m_bc.params()->template get<double>("Value X"));
    if (basis_->dimension() > 1) p.set("Value Y", this->m_bc.params()->template get<double>("Value Y"));
    if (basis_->dimension() > 2) p.set("Value Z", this->m_bc.params()->template get<double>("Value Z"));
    this->template registerEvaluator(fm, Teuchos::rcp(new panzer::ConstantVector<EvalT>(p)));
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "Unsupported basis type for flujo::BCStrategy_Dirichlet_Constant");
  }
}

}  // namespace flujo

#endif /** __Flujo_BCStrategy_Dirichlet_Constant_impl_HPP__ */
