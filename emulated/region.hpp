#pragma once


namespace qsim::emu {
    
    /*
     * mapping domain (x1, x2, ...., xM)
     */
    class region : private std::vector<double> {
        friend class wave_packet;
    public:
        region(std::initializer_list<double> init = {}) 
            : region(init) {
            std::sort(begin(), end());
        }

        inline size_t size() const {
            return std::vector<double>::size() + 1;
        }

        void emplace(double x) {
            auto place = std::lower_bound(begin(), end(), x);
            std::vector<double>::emplace(place, x);
        }

        size_t find(double x) {
            auto place = std::lower_bound(begin(), end(), x);
            return place - begin();
        }

        using std::vector<double>::erase;
    };
}
