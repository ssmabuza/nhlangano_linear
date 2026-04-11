// driver.hpp

#ifndef DRIVER_HPP
#define DRIVER_HPP

#include <stdexcept>

class Driver {
public:
    virtual void start() {
        throw std::runtime_error("Start method not implemented");
    }

    virtual void stop() {
        throw std::runtime_error("Stop method not implemented");
    }

    virtual ~Driver() {}
};

#endif // DRIVER_HPP