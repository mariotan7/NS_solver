#ifndef SUMUPRHS_H
#define SUMUPRHS_H

#include "Mesh.h"

template < int N >
class SumUpRHS
{
public:

    SumUpRHS
        (
            const real in_dt,
            const real * const * const in_rhs,
            const real * const in_coef

        ) : dt( in_dt ), rhs( in_rhs ), coef( in_coef ){};


    real operator() ( const int id ) const
    {
        real sum_rhs = 0.0;

        for ( int n = 0 ; n < N ; n++ ) {
            sum_rhs += dt * ( coef[n] ) * ( rhs[n][id] );
        }

        return sum_rhs;
    }

private:

    const real dt;
    const real * const * const rhs;
    const real * const coef;
};

#endif //SUMUPRHS_H
