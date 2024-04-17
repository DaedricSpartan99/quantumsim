#pragma once

namespace qsim {
    
    // use a specific coordinate system 
    template<typename ... Coords>
    class potential {
    public:
        virtual double operator()(const Coords&... access) const = 0;
    };
}
