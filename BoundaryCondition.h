#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include "Mesh.h"

class BCondition
{
public:

    BCondition( const Mesh &in_mesh ) : mesh( in_mesh ) {}
    virtual ~BCondition(void) {};

    virtual real operator()
    (
      const real * const v,
      const int id,
      const int bi, const int bj, const int bk
    ) const = 0;

protected:
    const Mesh &mesh;
};

class Dirichlet : public BCondition
{
public:

    //------------------------------------------------------------------
    // constructor & destructor
    //------------------------------------------------------------------
    Dirichlet ( const Mesh & in_mesh, const real in_value )
        : BCondition( in_mesh ), value( in_value ) {}

    ~Dirichlet( void ) {}

    //------------------------------------------------------------------
    // operator
    //------------------------------------------------------------------
    real operator()
    (
      const real * const v,
      const int id,
      const int bi, const int bj, const int bk
    ) const
    {
        return value;
    }

private:
    const real value;
};


class Neumann_plus : public BCondition
{
public:

    //------------------------------------------------------------------
    // constructor & destructor
    //------------------------------------------------------------------
    Neumann_plus( const Mesh & in_mesh, const int in_d )
        : BCondition( in_mesh ), d( in_d ) {}

    ~Neumann_plus( void ) {}

    //------------------------------------------------------------------
    // operator
    //------------------------------------------------------------------
    real operator()
    (
      const real * const v,
      const int id,
      const int bi, const int bj, const int bk
    ) const
    {
        return v[id + mesh.ST(d)];
    }

private:
    const  int    d;
};

class Neumann_minus : public BCondition
{
public:

    //------------------------------------------------------------------
    // constructor & destructor
    //------------------------------------------------------------------
    Neumann_minus( const Mesh & in_mesh, const int in_d )
        : BCondition( in_mesh ), d( in_d ) {}

    ~Neumann_minus( void ) {}

    //------------------------------------------------------------------
    // operator
    //------------------------------------------------------------------
    real operator()
    (
      const real * const v,
      const int id,
      const int bi, const int bj, const int bk
    ) const
    {
        return v[id - mesh.ST(d)];
    }

private:
    const  int    d;
};



#endif //BOUNDARYCONDITION_H
