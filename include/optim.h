#ifndef OPTIM_H
#define OPTIM_H

#include <vector>

#include "config.h"

class Optim
{
    public:
        Optim(Config& input);
        virtual ~Optim();
        int run();
        double en();
        void min(Config& output);
    protected:
    private:
        std::vector<Atom> coords;
        double m_en;
};

#endif // OPTIM_H
