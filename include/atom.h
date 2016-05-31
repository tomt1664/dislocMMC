#ifndef ATOM_H
#define ATOM_H


class Atom
{
    public:
        Atom();
        Atom(double x, double y, double z);
        virtual ~Atom();
        void setPos(double x, double y, double z);
        double getx() { return m_x; }
        double gety() { return m_y; }
        double getz() { return m_z; }
    private:
        double m_x;
        double m_y;
        double m_z;
};

#endif // ATOM_H
