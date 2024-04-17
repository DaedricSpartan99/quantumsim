#pragma once 

#include "wave_packet.hpp"

namespace qsim::emu {

    class wave : public std::vector<std::pair<double, wave_packet>> {

        std::shared_ptr<region> map; 
        std::shared_ptr<potential<size_t>> V; 
        
        using std::vector<std::pair<double, wave_packet>>::push_back;

    public:

        wave(std::shared_ptr<region> map, std::shared_ptr<potential<size_t>> V);
        
        void push_back(double energy);
        
        /*
         * normalize with respect to the 
         * region length
         */
        void normalize();

        /*
         * determine A and B in function of the potential
         */
        void eval_coefficients();
        
        /*
         * integrals
         */
        double position() const;
        double momentum() const;
        double energy() const;
    };
}
