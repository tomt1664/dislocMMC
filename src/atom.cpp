#include <cmath>
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

double Atom::dist(Atom at) //calculate distance to atom at
{
    double bdist;
    double xd = m_x - at.getx();
    double yd = m_y - at.gety();
    double zd = m_z - at.getz();
    bdist = sqrt(xd*xd+yd*yd+zd*zd);
    return bdist;
}
