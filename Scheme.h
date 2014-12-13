#ifndef SCHEME_H
#define SCHEME_H

#include "Mesh.h"

class Scheme
{
public:

    Scheme( const Mesh &in_mesh ) : mesh( in_mesh ) {}
    virtual ~Scheme(void) {}

    virtual real operator()
        ( const real * const * const f, const int id ) const = 0;

protected:
    const Mesh &mesh;
};


template < int D >
class SecondCenter : public Scheme
{
public:
    SecondCenter ( const Mesh &in_mesh ) : Scheme( in_mesh )
    {
        for ( int d = 0 ; d < D ; d++ ) {
            inv_dxdx[d] = 1.0 / ( mesh.delta(d) * mesh.delta(d) );
        }
    }
    ~SecondCenter( void ) {}

    real operator() ( const real * const * const a_f, const int id ) const
    {
        real dfdt = 0.0;
        const real * const f = a_f[0];

        for ( int d = 0 ; d < D ; d++ ) {

            const int   c = id;

            const int   p1 = id + 1 * mesh.ST(d);
            const int   m1 = id - 1 * mesh.ST(d);

            //const real dxdx = mesh.delta(d) * mesh.delta(d);
            //const real inv_dxdx = 1.0 / mesh.delta(d) * mesh.delta(d);

            //dfdt += ( f[p1] - 2.0 * f[c] + f[m1] ) / dxdx;
            dfdt += ( f[p1] - 2.0 * f[c] + f[m1] ) * inv_dxdx[d];
        }

        return  dfdt;
    }

private:
    real inv_dxdx[D];
};

template < int D >
class ForthCenter : public Scheme
{
public:
    ForthCenter ( const Mesh &in_mesh ) : Scheme( in_mesh )
    {
        for ( int d = 0 ; d < D ; d++ ) {
            inv_dxdx[d] = 1.0 / ( 24.0 * 24.0 * mesh.delta(d) * mesh.delta(d) );
        }
    }

    ~ ForthCenter( void ) {}

    real operator() ( const real * const * const a_f, const int id ) const
    {
        real dfdt = 0.0;
        const real * const f = a_f[0];

        for ( int d = 0 ; d < D ; d++ ) {

            const int   c = id;

            const int   p1 = id + 1 * mesh.ST(d);
            const int   m1 = id - 1 * mesh.ST(d);

            const int   p2 = id + 2 * mesh.ST(d);
            const int   m2 = id - 2 * mesh.ST(d);

            const int   p3 = id + 3 * mesh.ST(d);
            const int   m3 = id - 3 * mesh.ST(d);


            //const real inv_dxdx = 1.0 / (24.0 * 24.0 * mesh.delta(d) * mesh.delta(d));

            dfdt += (
                      + 1.0    * f[m3] - 54.0 * f[m2] + 783.0 * f[m1]
                      - 1460.0 * f[c]
                      + 1.0    * f[p3] - 54.0 * f[p2] + 783.0 * f[p1]

                    ) * inv_dxdx[d];

        }

        return  dfdt;
    }

private:
    real inv_dxdx[D];
};

#endif
