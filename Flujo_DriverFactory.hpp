// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_DriverFactory_HPP__
#define __Flujo_DriverFactory_HPP__

#include <Teuchos_RCP.hpp>
#include "Flujo_Driver.hpp"


namespace flujo {

class DriverFactory {

public:
  /**
   * @brief Build a driver based on the input parameters and the communicator.
   * 
   * @param input_params The input parameters.
   * @param comm The communicator.
   * @return Teuchos::RCP<Driver> The driver.
   */
  static Teuchos::RCP<Driver> build(Teuchos::RCP<Teuchos::ParameterList> input_params, 
                                    Teuchos::RCP<const Teuchos::Comm<int>> comm) 
  { return Teuchos::null; }

};

}
// end namespace flujo

#endif /** __Flujo_DriverFactory_HPP__ */