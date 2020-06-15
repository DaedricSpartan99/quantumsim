#pragma once

#include "potential.hpp"
#include <list>

namespace qsim::pot {
       
    template <typename Coords> 
    class composed : public std::list<std::shared_ptr<potential<Coords>>>, virtual public potential<Coords> {
    public:

        composed(std::initializer_list<potential<Coords>*> list = {}) {
            reserve(list.size());
            for (auto * pot : list)
                push_back(std::shared_ptr<potential<Coords>>(pot));;
        }
        
        virtual double operator()(const Coords& access) const override {
            double out(0);

            for (const auto& pot : *this)
                out += pot(access);

            return out;
        }
    }; 
}
