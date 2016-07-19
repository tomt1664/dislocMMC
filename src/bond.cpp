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

//set indicies
void Bond::setindex(long i1, long i2)
{
    m_i1 = i1;
    m_i2 = i2;
}

//return an atom object containing the midpoint of the bond
Atom Bond::midpoint(double bounds[5])
{
    //implement periodic boundary conditions
    double xadd = 0.0;
    double yadd = 0.0;
    double zadd = 0.0;

    //get cell dimensions
    double perx = bounds[1] - bounds[0];
    double pery = bounds[3] - bounds[2];
    double perz = bounds[5] - bounds[4];

    //check for boundary crossing
    if(m_atom1.getx() < (bounds[0]+0.25*perx) && m_atom2.getx() > (bounds[0]+0.75*perx)) xadd = -perx;
    if(m_atom1.getx() > (bounds[0]+0.75*perx) && m_atom2.getx() < (bounds[0]+0.25*perx)) xadd = perx;

    if(m_atom1.gety() < (bounds[2]+0.25*pery) && m_atom2.gety() > (bounds[2]+0.75*pery)) yadd = -pery;
    if(m_atom1.gety() > (bounds[2]+0.75*pery) && m_atom2.gety() < (bounds[2]+0.25*pery)) yadd = pery;

    if(m_atom1.getz() < (bounds[4]+0.25*perz) && m_atom2.getz() > (bounds[4]+0.75*perz)) zadd = -perz;
    if(m_atom1.getz() > (bounds[4]+0.75*perz) && m_atom2.getz() < (bounds[4]+0.25*perz)) zadd = perz;

    double xt,yt,zt;
    xt = m_atom2.getx() + (m_atom1.getx() - xadd - m_atom2.getx())/2.0;
    yt = m_atom2.gety() + (m_atom1.gety() - yadd - m_atom2.gety())/2.0;
    zt = m_atom2.getz() + (m_atom1.getz() - zadd - m_atom2.getz())/2.0;
    Atom atom(xt,yt,zt);
    return atom;
}

Bond::~Bond()
{
    //dtor
}
