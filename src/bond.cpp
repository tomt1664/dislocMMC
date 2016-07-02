#include "atom.h"
#include "bond.h"

// create bond
Bond::Bond(Atom at1, Atom at2, long i1, long i2)
{
    m_atom1 = at1;
    m_atom2 = at2;
    m_i1 = i1;
    m_i2 = i2;
}

//return coordinates
void Bond::coords(Atom& at1, Atom& at2)
{
    at1 = m_atom1;
    at2 = m_atom2;
}

//return indicies
void Bond::index(long& i1, long& i2)
{
    i1 = m_i1;
    i2 = m_i2;
}

Bond::~Bond()
{
    //dtor
}
