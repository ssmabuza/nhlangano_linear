// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_DriverFractory_HPP__
#define __Flujo_DriverFractory_HPP__

#include <Teuchos_RCP.hpp>
#include "Flujo_Driver.hpp"


namespace flujo {

class DriverFactory {

  static Teuchos::RCP<Driver> build() {}



};

}
// end namespace flujo