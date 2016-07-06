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

//return an atom object containing the midpoint of the bond
Atom Bond::midpoint()
{
    double xt,yt,zt;
    xt = m_atom2.getx() + (m_atom1.getx() - m_atom2.getx())/2.0;
    yt = m_atom2.gety() + (m_atom1.gety() - m_atom2.gety())/2.0;
    zt = m_atom2.getz() + (m_atom1.getz() - m_atom2.getz())/2.0;
    Atom atom(xt,yt,zt);
    return atom;
}

Bond::~Bond()
{
    //dtor
}
