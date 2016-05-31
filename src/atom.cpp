#include "atom.h"

Atom::Atom()
{

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

void Atom::setPos(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}
