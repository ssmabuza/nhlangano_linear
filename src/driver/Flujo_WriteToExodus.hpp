// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_WriteToExodus_HPP__
#define __Flujo_WriteToExodus_HPP__

#include <Teuchos_RCP.hpp>
#include <Thyra_VectorBase.hpp>
#include <Panzer_ModelEvaluator.hpp>
#include <Panzer_ResponseLibrary.hpp>
#include <Panzer_STK_Interface.hpp>

namespace flujo
{

  template <class Scalar>
  inline
  void WriteToExodus(double time_stamp,
                     const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &x,
                     const panzer::ModelEvaluator<Scalar> &model,
                     panzer::ResponseLibrary<panzer::Traits> &stkIOResponseLibrary,
                     panzer_stk::STK_Interface &mesh)
  {
    // fill STK mesh objects
    Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = model.createInArgs();
    inArgs.set_x(x);
    inArgs.set_t(time_stamp);

    panzer::AssemblyEngineInArgs respInput;
    model.setupAssemblyInArgs(inArgs, respInput);

    stkIOResponseLibrary.addResponsesToInArgs<panzer::Traits::Residual>(respInput);
    stkIOResponseLibrary.evaluate<panzer::Traits::Residual>(respInput);

    mesh.writeToExodus(time_stamp);
  }

}
// end namespace flujo

#endif /* __Flujo_WriteToExodus_HPP__ */