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

#include <cmath>
#include <iostream>
#include "atom.h"

Atom::Atom() //default constructor to create atom centered at the origin
{
    m_x = 0.0;
    m_y = 0.0;
    m_z = 0.0;
}

Atom::Atom(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

Atom::~Atom()
{
    //dtor
}

void Atom::setPos(double x, double y, double z) //change the position of an existing atom
{
    m_x = x;
    m_y = y;
    m_z = z;
}

double Atom::dist(Atom at, double bounds[5]) //calculate distance to atom at - applying the periodic boundary conditions
{
    double bdist;
    double xadd = 0.0;
    double yadd = 0.0;
    double zadd = 0.0;

    //get cell dimensions
    double perx = bounds[1] - bounds[0];
    double pery = bounds[3] - bounds[2];
    double perz = bounds[5] - bounds[4];



    //check for boundary crossing
    if(m_x < (bounds[0]+0.25*perx) && at.getx() > (bounds[0]+0.75*perx)) xadd = -perx;
    if(m_x > (bounds[0]+0.75*perx) && at.getx() < (bounds[0]+0.25*perx)) xadd = perx;

    if(m_y < (bounds[2]+0.25*pery) && at.gety() > (bounds[2]+0.75*pery)) yadd = -pery;
    if(m_y > (bounds[2]+0.75*pery) && at.gety() < (bounds[2]+0.25*pery)) yadd = pery;

    if(m_z < (bounds[4]+0.25*perz) && at.getz() > (bounds[4]+0.75*perz)) zadd = -perz;
    if(m_z > (bounds[4]+0.75*perz) && at.getz() < (bounds[4]+0.25*perz)) zadd = perz;

    double xd = m_x - (at.getx()+xadd);
    double yd = m_y - (at.gety()+yadd);
    double zd = m_z - (at.getz()+zadd);
    bdist = sqrt(xd*xd+yd*yd+zd*zd);
    return bdist;
}
