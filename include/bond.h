/****************************************************
**  dislocMMC
**
**  A Metropolis Monte Carlo algorithm to model the
**  evolution of structure in defective sp2 bonded
**  carbon systems
**
**  The code calls the LAMMPS atomistic simulation
**  program for structural optimisation with a
**  reactive force-field
**
**  T.Trevethan 2016
*****************************************************/

#ifndef BOND_H
#define BOND_H

#include "atom.h"

//class to store bond objects (formed from two atoms)

class Bond
{
    public:
        Bond(Atom at1, Atom at2, long i1, long i2);
        virtual ~Bond();
        void coords(Atom& at1, Atom& at2);
        void index(long& i1, long& i2);
        void setindex(long i1, long i2);
        Atom midpoint(double bounds[5]);
    private:
        Atom m_atom1; //coordinates of atom 1
        Atom m_atom2; //coordinates of atom 2
        long m_i1; //index of atom 1
        long m_i2; //index of atom 2
};

#endif // BOND_H
