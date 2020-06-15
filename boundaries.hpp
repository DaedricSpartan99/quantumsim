#pragma once

namespace qsim {

    /*
     * Interval class, useful for the 1D case
     */
    class interval {
    protected:
        double lower, upper; 
    public:

        interval(double low, double up)
            : lower(low), upper(up) {}

        // continuity condition
        virtual std::pair<wave_t, wave_t> continuity() const = 0;

        // wave function value
        virtual wave_t before()(double) const = 0;
        virtual wave_t after()(double) const = 0;

        /*
         * Integrals
         */
        
        // return the normalization factor before + after the boundary
        virtual double normalization() const = 0;
        virtual double position() const = 0;
        virtual double momentum() const = 0;
        virtual double energy() const = 0;
    };
}
