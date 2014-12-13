#ifndef SOR_H
#define SOR_H

#include <cmath>
#include "Variable.h"
#include "Stencils.h"

typedef double real;

class LHSofLinearSimultaneousEQ
{
public:

    LHSofLinearSimultaneousEQ( void ) {}
    virtual ~LHSofLinearSimultaneousEQ( void ) {}
    virtual real operator () ( const real * const phi, const int i, const int j, const int k ) const = 0;
};


class LHSofPoissonEQ : public LHSofLinearSimultaneousEQ
{
public:
     LHSofPoissonEQ( const Mesh & in_mesh )
         : LHSofLinearSimultaneousEQ(), mesh( in_mesh ),
           inv_dxdx( 1.0 / ( mesh.delta(0) * mesh.delta(0) ) ),
           inv_dydy( 1.0 / ( mesh.delta(1) * mesh.delta(1) ) ) {}

    ~LHSofPoissonEQ( void ) {}

    real
    operator () ( const real * const phi, const int i, const int j, const int k ) const
    {
        const int id = mesh.id(i,j,k);

        const int st0 = mesh.ST(0);
        const int st1 = mesh.ST(1);

        const real lhs =  ( phi[id-st0] - 2.0 * phi[id] + phi[id+st0] ) * inv_dxdx
                        + ( phi[id-st1] - 2.0 * phi[id] + phi[id+st1] ) * inv_dydy;

        return lhs;
    }

private:
    const Mesh & mesh;

    const real inv_dxdx;
    const real inv_dydy;
};


class RHSofLinearSimultaneousEQ
{
public:

    RHSofLinearSimultaneousEQ( void ) {}
    virtual ~RHSofLinearSimultaneousEQ( void ) {}

    virtual real operator () ( const int i, const int j, const int k ) const = 0;
};


class RHSofPoissonEQ : public RHSofLinearSimultaneousEQ
{
public:
     RHSofPoissonEQ( const Mesh & in_mesh, const real * const * const in_velocity, const real in_dt )
         : RHSofLinearSimultaneousEQ(), mesh( in_mesh ), velocity( in_velocity ), dt( in_dt ),
           inv_dx( 1.0 / mesh.delta(0) ), inv_dy( 1.0 / mesh.delta(1) ) {}
    ~RHSofPoissonEQ( void ) {}

    real
    operator () ( const int i, const int j, const int k ) const
    {
        const int id = mesh.id(i,j,k);

        const int st0 = mesh.ST(0);
        const int st1 = mesh.ST(1);

        const real * u = velocity[0];
        const real * v = velocity[1];

        const real rhs = ( u[id] - u[id-st0] ) * inv_dx + ( v[id] - v[id-st1] ) * inv_dy;

        //return rhs * dt;
        return rhs / dt;
    }

private:
    const Mesh & mesh;
    const real * const * const velocity;
    const real dt;

    const real inv_dx;
    const real inv_dy;
};

template <int D>
class SOR {

public:

    SOR
    (
        real * const                        in_phi,
        const BCondition * const * const    in_bc_set,
        const LHSofLinearSimultaneousEQ    &in_lhs,
        const RHSofLinearSimultaneousEQ    &in_rhs,
        const Mesh &                        in_mesh,
        const unsigned                      in_max_iteration = 10000,
        const real                          in_alpha         = 1.7
    ) : phi( in_phi ), bc_set( in_bc_set ), lhs( in_lhs ), rhs( in_rhs ), mesh( in_mesh ),
        max_iteration( in_max_iteration ), alpha( in_alpha ){}

    void
    solve( const real error_limit = 1e-6 )
    {
        numberOfIter = 0;

        const real B_ij = get_B_ij();

        do {
                error = 0.0;

                const int k = mesh.HN() ;
                //for ( int k = mesh.HN() ; k < mesh.HNF(2) ; k++ ) {
                for ( int j = mesh.HN() ; j < mesh.HNF(1) ; j++ ) {
                for ( int i = mesh.HN() ; i < mesh.HNF(0) ; i++ ) {

                    const int   id = mesh.id(i,j,k);

                    const real adjust =  1.0 / B_ij * ( lhs( phi, i,j,k ) - rhs(i,j,k) );

                    phi[id] += alpha * adjust;

                    if ( error < fabs(adjust) ) {
                        error = fabs(adjust);
                    }

                }}

                CellCenter::update_boundary<D>( phi, bc_set, mesh );
                numberOfIter++;

                if ( numberOfIter % 10000 == 0 ) {
                    std::cout << "iter  = " << numberOfIter;
                    std::cout << "\terror = " << error        << std::endl;
                }

        } while ( error > error_limit && numberOfIter <= max_iteration );

        //std::cout << "iter  = " << numberOfIter;
        //std::cout << "\terror = " << error        << std::endl;
    }

    real
    get_NumberOfIter() { return numberOfIter; }


private:
    
    real * const  phi;

    const BCondition * const * const bc_set;

    const LHSofLinearSimultaneousEQ &lhs;
    const RHSofLinearSimultaneousEQ &rhs;

    const Mesh &        mesh;

    const unsigned      max_iteration;
    const real          alpha;

    unsigned numberOfIter;
    real     error;

    real
    get_B_ij() const
    {
        real B_ij = 0.0;

        for ( int d = 0 ; d < D ; d++ ) {
            B_ij += 2.0 / ( mesh.delta(d) * mesh.delta(d) );
        }

        return B_ij;
    }

};


#endif //SOR_H
