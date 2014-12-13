#ifndef NEWSCHEME_H
#define NEWSCHEME_H
template < int D >
class Diffusion
{
public:
    Diffusion ( const real in_mu, const Mesh &in_mesh )
        : mu( in_mu ), mesh( in_mesh )
    {
        for ( int d = 0 ; d < D ; d++ ) {
            inv_dxdx[d] = 1.0 / ( mesh.delta(d) * mesh.delta(d) );
        }
    }
    ~Diffusion ( void ) {}

    real operator() ( const real * const f ) const
    {
        real dfdt = 0.0;

        for ( int d = 0 ; d < D ; d++ ) {

            const int   p1 = + mesh.ST(d);
            const int   m1 = - mesh.ST(d);

            dfdt += mu * ( f[p1] - 2.0 * f[0] + f[m1] ) * inv_dxdx[d];
        }
        return  dfdt;
    }

private:
    const real mu;
    const Mesh & mesh;
    real inv_dxdx[D];
};

class Pressure
{
public:
    Pressure( const int in_d, const Mesh & in_mesh )
        : d( in_d ), mesh( in_mesh ), inv_dx( 1.0 / mesh.delta(d) ) {}

    real operator() ( const real * const p ) const
    {
        return - ( p[mesh.ST(d)] - p[0] ) * inv_dx;
    }

private:
    const int d;
    const Mesh mesh;
    const real inv_dx;
};

class CorrectVelocity
{
public:
    CorrectVelocity
    (
      const real * const in_phi,
      const int    in_d,
      const real   in_dt,
      const Mesh & in_mesh
    ) : phi ( in_phi  ), d( in_d ), dt( in_dt ),
        mesh( in_mesh ), inv_dx( 1.0 / mesh.delta(d) ) {}

    real operator () ( const real * const vel, const int id ) const
    {
        const real dvel = - dt * ( phi[id+mesh.ST(d)] - phi[id] ) * inv_dx;

        return vel[id] + dvel;
    }

private:

    const real * const phi;
    const int  d;
    const real dt;
    const Mesh & mesh;
    const real inv_dx;
};

#endif //NEWSCHEME_H
