// ==============================================================================
//                Flujo: Copyright Valiant Scientific
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ==============================================================================

#ifndef __Flujo_Driver_HPP__
#define __Flujo_Driver_HPP__

#include <stdexcept>


namespace flujo {

/**
 * @brief A base class for all driver classes.
 * 
 */
class Driver {
public:

    /**
     * @brief Set up the driver.
     */
    virtual void setup() {
        throw std::runtime_error("Setup method not implemented in Driver");
    }

    /**
     * @brief Solve the driver.
     */
    virtual void solve() const {
        throw std::runtime_error("Solve method not implemented in Driver");
    }

    /**
     * @brief Destructor for the driver.
     */
    virtual ~Driver() {}
};

} // end namespace flujo

#endif /** __Flujo_Driver_HPP__ */