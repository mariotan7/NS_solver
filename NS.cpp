#include "NS.h"

#include "initial_value.h"
#include "SOR.h"
#include "NewScheme.h"

#include <cmath>

Iter
time_marching_RK1_cavity
(
    real * const nu,
    real * const nv,
    real * const p,
    const real * const u,
    const real * const v,
    const real dt,
    const Mesh & mesh
){
    const real Re = 1000.0;
    //const real Re = 1.0;

    calc_nu( nu, u, v, p, dt, Re, u, mesh );
    calc_nv( nv, u, v, p, dt, Re, v, mesh );

    boundary_update_u_v( nu, nv, mesh );

    static int i = 0;
    real * phi;
    if ( i == 0 ) {
        phi = new real[mesh.size()];
        make_initial_value_0( phi, mesh );
        i++;
    }

    calc_phi( phi, nu, nv, dt, mesh );

    const Iter iter = calc_phi( phi, nu, nv, dt, 1.5, 10000, 1.0e-4, mesh );
    correct_u_v_p( nu, nv, p, phi, dt, mesh );

    //delete[] phi; phi = NULL;

    boundary_update_u_v( nu, nv, mesh );

    return iter;
}

void
time_marching_RK2_cavity
(
    real * * const pu,
    real * * const pv,
    real * * const pp,
    const real dt,
    const Mesh & mesh
){
    const real Re = 3200.0;

    real * & u_0 = *pu;
    real * & v_0 = *pv;
    real * & p   = *pp;

    //--------------------------
    // RK 1
    //--------------------------
    const real dt_1 = dt / 2.0;
    real * u_1 = new real[mesh.size()];
    real * v_1 = new real[mesh.size()];

    calc_nu( u_1, u_0, v_0, p, dt_1, Re, u_0, mesh );
    calc_nv( v_1, u_0, v_0, p, dt_1, Re, v_0, mesh );
    boundary_update_u_v( u_1, v_1, mesh );

    real * phi = new real[mesh.size()];

    calc_phi( phi, u_1, v_1, dt_1, mesh );
    correct_u_v_p( u_1, v_1, p, phi, dt_1, mesh );

    boundary_update_u_v( u_1, v_1, mesh );

    //--------------------------
    // RK 2
    //--------------------------
    const real dt_2 = dt;
    real * u_2 = new real[mesh.size()];
    real * v_2 = new real[mesh.size()];

    calc_nu( u_2, u_1, v_1, p, dt_2, Re, u_0, mesh );
    calc_nv( v_2, u_1, v_1, p, dt_2, Re, v_0, mesh );

    boundary_update_u_v( u_2, v_2, mesh );

    calc_phi( phi, u_2, v_2, dt_2, mesh );
    correct_u_v_p( u_2, v_2, p, phi, dt_2, mesh );

    boundary_update_u_v( u_2, v_2, mesh );

    delete[] u_0;
    delete[] v_0;
    delete[] u_1; u_1 = NULL;
    delete[] v_1; v_1 = NULL;
    delete[] phi; phi = NULL;

    u_0 = u_2; u_2 = NULL;
    v_0 = v_2; v_2 = NULL;
}

//------------------------
// calc nu
//------------------------
void
calc_nu
(
          real * const nu,
    const real * const  u,
    const real * const  v,
    const real * const  p,
    const real dt,
    const real Re,
    const real * const  u0,
    const Mesh & mesh
){
    //Diffusion<2> sc( 1.0 / Re, mesh );
    //Pressure pr( 0, mesh );

    const int k = mesh.hn() ;
    for ( int j = mesh.hn() ; j < mesh.hnf(1)     ; j++ ) {
    for ( int i = mesh.hn() ; i < mesh.hnf(0) - 1 ; i++ ) {

        const int id = mesh.id(i,j,k);

        //const real adv = calc_adv( &u[id], &v[id], mesh.ST(0), mesh.ST(1), mesh.delta(0), mesh.delta(1) );
        //const real prs = pr( &p[id] );
        //const real dif = sc( &u[id] );

        const real adv = calc_adv_u( u, v, id, mesh );
        const real dif = calc_dif_u( u, Re, id, mesh);
        const real prs = calc_prs_u( p, id, mesh );

        const real rhs = adv + prs + dif;
        //const real rhs = prs + dif;
        nu[id] = u0[id] + dt * rhs;

    }}
}

//------------------------
// calc nv
//------------------------
void
calc_nv
(
          real * const nv,
    const real * const  u,
    const real * const  v,
    const real * const  p,
    const real dt,
    const real Re,
    const real * const  v0,
    const Mesh & mesh
){
    //Diffusion<2> sc( 1.0 / Re, mesh );
    //Pressure     pr( 1, mesh );

    const int k = mesh.hn() ;
    for ( int j = mesh.hn() ; j < mesh.hnf(1) - 1 ; j++ ) {
    for ( int i = mesh.hn() ; i < mesh.hnf(0)     ; i++ ) {

        const int id = mesh.id(i,j,k);

        //const real adv = calc_adv( &v[id], &u[id], mesh.ST(1), mesh.ST(0), mesh.delta(1), mesh.delta(0) );
        //const real prs = pr( &p[id] );
        //const real dif = sc( &v[id] );

        const real adv = calc_adv_v( u, v, id, mesh );
        const real prs = calc_prs_v( p, id, mesh );
        const real dif = calc_dif_v( v, Re, id, mesh );

        const real rhs = adv + prs + dif;
        //const real rhs = prs + dif;
        nv[id] = v0[id] + dt * rhs;
    }}
}

//------------------------
// calc advection term
//------------------------
real
calc_adv
(
    const real * const u,
    const real * const v,
    const int st0,
    const int st1,
    const real dx,
    const real dy
){
    const real term_0 = intrp( &u[0]   , +st0 );
    const real term_1 = intrp( &u[0]   , -st0 );

    const real term_2 = intrp( &v[0]   , +st0 );
    const real term_3 = intrp( &v[-st1], +st0 );

    const real term_4 = intrp( &u[0]   , +st1 );
    const real term_5 = intrp( &u[0]   , -st1 );

    const real adv_test =  1.0 / dx * ( ( term_0 * term_0 ) - ( term_1 * term_1 ) )
                         + 1.0 / dy * ( ( term_2 * term_4 ) - ( term_3 * term_5 ) );
    return -adv_test;
}


//------------------------
// calc phi
//------------------------
Iter
calc_phi
(
          real * const phi,
    const real * const u,
    const real * const v,
    const real dt,
    const real beta,
    const int  iter_max,
    const real error_condition,
    const Mesh & mesh
){
    int  iter;
    real  err = 1e+10;

    for ( iter = 1; iter < iter_max && err > error_condition ; ++iter ) {
        err = calc_phi_SOR( phi, u, v, dt, beta, mesh );
    }

    return Iter(iter, err);
}

real
calc_phi_SOR
(
          real * const phi,
    const real * const u,
    const real * const v,
    const real dt,
    const real beta,
    const Mesh & mesh
){
    real sum_for_numerator   = 0.0;
    real sum_for_denominator = 0.0;

    boundary_update_phi( phi, mesh );

    const int k = mesh.hn() ;
    for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
    for ( int i = mesh.hn() ; i < mesh.hnf(0) ; i++ ) {

        const int id = mesh.id(i,j,k);

        //if ( i != mesh.hn() || j != mesh.hn() ) {


            const int id_e = mesh.id(i+1,j,k);
            const int id_w = mesh.id(i-1,j,k);

            const int id_n = mesh.id(i,j+1,k);
            const int id_s = mesh.id(i,j-1,k);

            const real B_x  = 1.0 / std::pow(mesh.delta(0), 2);
            const real B_y  = 1.0 / std::pow(mesh.delta(1), 2);

            const real B_ij = 2.0 * ( B_x + B_y );

            const real psi = psi_ij( u, v, id, dt, mesh );

            const real nabla_phi =   B_y  * ( phi[id_s] + phi[id_n] )
                                   + B_x  * ( phi[id_w] + phi[id_e] )
                                   - B_ij * phi[id];

            const real adjust = nabla_phi - psi;
            
            phi[id] += beta / B_ij * adjust;

            sum_for_numerator   += std::pow(nabla_phi - psi, 2);
            sum_for_denominator += std::pow(psi, 2);

        //} else {
            //phi[id] = 0.0;
        //}
    }}

    boundary_update_phi( phi, mesh );

    const int nn = mesh.n(0) * mesh.n(1);
    const real numerator   = calc_norm( sum_for_numerator, nn );
    const real denominator = calc_norm( sum_for_denominator, nn );

    return numerator / denominator;
}

real calc_norm
(
    const real value,
    const int  n
){
    return std::sqrt(value / static_cast<real>(n));
}


real calc_B
(
    const real dx,
    const int  i_e,
    const int pi_e,
    const int  i_w,
    const int pi_w
){
    return ( i_e != pi_e && i_w != pi_w ) ? (1.0 / dx / dx) : 0.0 ;
}


void
boundary_update_phi
(
    real * const phi,
    const Mesh & mesh
){
    const int k = mesh.hn();

    // west wall
    const int i_w = mesh.hn() - 1;
    for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
        const int id_w = mesh.id(i_w,j,k);
        const int id_i = mesh.id(i_w+1,j,k);
        phi[id_w] = phi[id_i];
    }

    // east wall
    const int i_e = mesh.hnf(0);
    for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
        const int id_e = mesh.id(i_e,j,k);
        const int id_i = mesh.id(i_e-1,j,k);
        phi[id_e] = phi[id_i];
    }

    // south wall
    const int j_s = mesh.hn() - 1;
    for ( int i = mesh.hn() ; i < mesh.hnf(0) ; i++ ) {
        const int id_s = mesh.id(i,j_s,k);
        const int id_i = mesh.id(i,j_s+1,k);
        phi[id_s] = phi[id_i];
    }

    // north wall
    const int j_n = mesh.hnf(1);
    for ( int i = mesh.hn() ; i < mesh.hnf(0) ; i++ ) {
        const int id_n = mesh.id(i,j_n,k);
        const int id_i = mesh.id(i,j_n-1,k);
        phi[id_n] = phi[id_i];
    }

}

real psi_ij
(
    const real * const u,
    const real * const v,
    const int id,
    const real dt,
    const Mesh & mesh
){
    return 1.0 / dt * calc_div_locol( &u[id], &v[id], mesh.st(0), mesh.st(1), mesh.delta(0), mesh.delta(1));
}

real calc_div_locol
(
    const real * const u,
    const real * const v,
    const int st0,
    const int st1,
    const real d0,
    const real d1
){
    return ( (u[0] - u[-st0]) / d0 + (v[0] - v[-st1]) / d1 );
}


void
calc_phi
(
          real * const phi,
    const real * const u,
    const real * const v,
    const real dt,
    const Mesh & mesh
){
    make_initial_value_0( phi, mesh );

    const real * const vel[2] = { u, v };

    LHSofPoissonEQ lhs( mesh );
    RHSofPoissonEQ rhs( mesh, vel, dt );

    Neumann_plus  neu_p_x( mesh, 0 );
    Neumann_minus neu_m_x( mesh, 0 );

    Neumann_plus  neu_p_y( mesh, 1 );
    Neumann_minus neu_m_y( mesh, 1 );

    BCondition * bc_set[6];

    bc_set[ 0] = &neu_p_x; bc_set[ 1] = &neu_m_x;
    bc_set[ 2] = &neu_p_y; bc_set[ 3] = &neu_m_y;
    bc_set[ 4] = NULL;     bc_set[ 5] = NULL;

    SOR<2> sor( phi, bc_set, lhs, rhs, mesh, 10000, 1.2 );
    sor.solve(1e-9);
}

//------------------------
// correct u, v, p
//------------------------
void
correct_u_v_p
(
    real * const u,
    real * const v,
    real * const p,
    const real * const phi,
    const real dt,
    const Mesh & mesh
){
    //CorrectVelocity correct_u_func( phi, 0, dt, mesh );
    //CellFaceX::for_each<CorrectVelocity>( u, correct_u_func, mesh );

    //CorrectVelocity correct_v_func( phi, 1, dt, mesh );
    //CellFaceY::for_each<CorrectVelocity>( v, correct_v_func, mesh );

    correct_u( u, phi, dt, mesh );
    correct_v( v, phi, dt, mesh );
    correct_p( p, phi    , mesh );
}

void copy_array
(
    real * const       cpy_u,
    const real * const u,
    const Mesh & mesh
){
    const int k = mesh.hn() ;
    for ( int j = 0 ; j < mesh.NN(1) ; j++ ) {
    for ( int i = 0 ; i < mesh.NN(0) ; i++ ) {
        const int id = mesh.id(i,j,k);
        cpy_u[id] = u[id];
    }}
}

void check_array
(
    const real * const a,
    const real * const b,
    const Mesh & mesh
){
    real max_error = 0.0;

    bool flag = true;

    const int k = mesh.hn() ;
    for ( int j = 0 ; j < mesh.NN(1) ; j++ ) {
    for ( int i = 0 ; i < mesh.NN(0) ; i++ ) {
        const int id = mesh.id(i,j,k);

        const real error = fabs( a[id] - b[id] );
        if ( max_error < error ) {
            max_error = error;
        }

        if ( a[id] != b[id] ) {
            flag = false;
        }
    }}

    if ( flag == false ) {
        std::cout << "incorrect!\t";
        std::cout << "max error = " << max_error << std::endl;
    } else {
        std::cout << "OK!" << std::endl;
    }

}


void
correct_u
(
    real * const u,
    const real * const phi,
    const real dt,
    const Mesh & mesh
){
    const int st0 = mesh.ST(0);

    const int k = mesh.hn() ;
    for ( int j = mesh.hn() ; j < mesh.hnf(1)     ; j++ ) {
    for ( int i = mesh.hn() ; i < mesh.hnf(0) - 1 ; i++ ) {
        const int id = mesh.id(i,j,k);
        u[id] = u[id] - dt * ( phi[id+st0] - phi[id] ) / mesh.delta(0);
    }}
}

void
correct_v
(
    real * const v,
    const real * const phi,
    const real dt,
    const Mesh & mesh
){
    const int st1 = mesh.ST(1);

    const int k = mesh.hn();
    for ( int j = mesh.hn() ; j < mesh.hnf(1) - 1 ; j++ ) {
    for ( int i = mesh.hn() ; i < mesh.hnf(0)     ; i++ ) {
        const int id = mesh.id(i,j,k);
        v[id] = v[id] - dt * ( phi[id+st1] - phi[id] ) / mesh.delta(1);
    }}
}

void
correct_p
(
    real * const p,
    const real * const phi,
    const Mesh & mesh
){
    const int k = mesh.hn();
    for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
    for ( int i = mesh.hn() ; i < mesh.hnf(0) ; i++ ) {
        const int id = mesh.id(i,j,k);

        p[id] = p[id] + phi[id];
        //p[id] -= phi[id];

    }}
}

void
boundary_update_only_v
(
    real * const v,
    const Mesh & mesh
){
    const int k = mesh.hn();

    // north and south wall
    for ( int i = mesh.hn(); i < mesh.hnf(0) ; i++ ) {
        const int j_b = mesh.hn()   - 1;
        const int j_t = mesh.hnf(1) - 1;

        const int id_b = mesh.id(i,j_b,k);
        const int id_t = mesh.id(i,j_t,k);

        v[id_b] = 0.0;
        v[id_t] = 0.0;
    }

    // west and east wall
    for ( int j = mesh.hn() ; j < mesh.hnf(1) - 1 ; j++ ) {
    //for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
        const int i_l = mesh.hn() - 1;
        const int i_r = mesh.hnf(0);

        const int id_l = mesh.id(i_l,j,k);
        const int id_r = mesh.id(i_r,j,k);

        v[id_l] = - v[id_l+mesh.ST(0)];
        v[id_r] = - v[id_r-mesh.ST(0)];
    }

}

void
boundary_update_only_u
(
    real * const u,
    const Mesh & mesh
){
    const int k = mesh.hn();

    // east wall
    const int i_e = mesh.hn() - 1;
    for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
        const int id = mesh.id(i_e,j,k);
        u[id] = 0.0;
    }

    // west wall
    const int i_w = mesh.hnf(0) - 1;
    for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
        const int id = mesh.id(i_w,j,k);
        u[id] = 0.0;
    }
 
    // south wall
    const int j_s  = mesh.hn() - 1;
    const int j_is = mesh.hn();
    for ( int i = mesh.hn() ; i < mesh.hnf(0) - 1; i++ ) {
        const int id_s  = mesh.id(i,j_s ,k);
        const int id_is = mesh.id(i,j_is,k);

        u[id_s] = - u[id_is];
    }

    // north wall
    const int j_n  = mesh.hnf(1);
    const int j_in = mesh.hnf(1) - 1;
    for ( int i = mesh.hn() ; i < mesh.hnf(0) - 1; i++ ) {
        const int id_n  = mesh.id(i,j_n ,k);
        const int id_in = mesh.id(i,j_in,k);

        //u[id_n] = - u[id_in];
        u[id_n] = 2.0 - u[id_in];
        //const real x = mesh.x(i,j_n,k);
        //u[id_n] = - 4.0 * x * (x - 1.0) - u[id_in];
    }

}

//------------------------
// boundary update u, v, p
//------------------------
void
boundary_update_u_v
(
    real * const u,
    real * const v,
    const Mesh &mesh
){
    boundary_update_only_u( u, mesh );
    boundary_update_only_v( v, mesh );
}

/*
void
boundary_update_u_v
(
    real * const u,
    real * const v,
    const Mesh &mesh
){
    const int k = mesh.hn();

    //--------------------
    // top and bottom of the wall
    //--------------------
    for ( int i = mesh.hn(); i < mesh.hnf(0) - 1 ; i++ ) {
    //for ( int i = mesh.hn(); i < mesh.hnf(0) ; i++ ) {
        const int j_b = mesh.hn() - 1;
        const int j_t = mesh.hnf(1);

        const int id_b = mesh.id(i,j_b,k);
        const int id_t = mesh.id(i,j_t,k);

        //u[id_b] = 2.0 - u[id_b+mesh.ST(1)];
        //u[id_t] = 0.0 - u[id_t-mesh.ST(1)];
        u[id_b] = 0.0 - u[id_b+mesh.ST(1)];
        u[id_t] = 2.0 - u[id_t-mesh.ST(1)];
    }

    for ( int i = mesh.hn(); i < mesh.hnf(0) ; i++ ) {
        const int j_b = mesh.hn()   - 1;
        const int j_t = mesh.hnf(1) - 1;

        const int id_b = mesh.id(i,j_b,k);
        const int id_t = mesh.id(i,j_t,k);

        v[id_b] = 0.0;
        v[id_t] = 0.0;
    }

    //--------------------
    // left and right of the wall
    //--------------------
    for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
        const int i_l = mesh.hn()   - 1;
        const int i_r = mesh.hnf(0) - 1;

        const int id_l = mesh.id(i_l,j,k);
        const int id_r = mesh.id(i_r,j,k);

        u[id_l] = 0.0;
        u[id_r] = 0.0;
    }

    for ( int j = mesh.hn() ; j < mesh.hnf(1) - 1 ; j++ ) {
    //for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
        const int i_l = mesh.hn() - 1;
        const int i_r = mesh.hnf(0);

        const int id_l = mesh.id(i_l,j,k);
        const int id_r = mesh.id(i_r,j,k);

        v[id_l] = - v[id_l+mesh.ST(0)];
        v[id_r] = - v[id_r-mesh.ST(0)];
    }
}
*/

real
calc_adv_u
(
    const real * const u,
    const real * const v,
    const int id,
    const Mesh & mesh
){
    const int st0 = mesh.ST(0);
    const int st1 = mesh.ST(1);

    const real adv   =  1.0 / mesh.delta(0) *
                              (   ( ( u[id+st0] + u[id] ) * 0.5 ) * ( ( u[id+st0] + u[id] ) * 0.5 )
                                - ( ( u[id-st0] + u[id] ) * 0.5 ) * ( ( u[id-st0] + u[id] ) * 0.5 )
                              )
                      + 1.0 / mesh.delta(1) *
                              (   ( ( v[id+st0    ] + v[id    ] ) * 0.5 ) * ( ( u[id+st1] + u[id] ) * 0.5 )
                                - ( ( v[id+st0-st1] + v[id-st1] ) * 0.5 ) * ( ( u[id-st1] + u[id] ) * 0.5 )  );

    return -adv;

}

real
calc_adv_v
(
    const real * const u,
    const real * const v,
    const int id,
    const Mesh & mesh
){
    const int st0 = mesh.ST(0);
    const int st1 = mesh.ST(1);

    const real adv   =  1.0 / mesh.delta(1) *
                            (   ( ( v[id    ] + v[id+st1] ) * 0.5 ) * ( ( v[id    ] + v[id+st1] ) * 0.5 )
                              - ( ( v[id-st1] + v[id    ] ) * 0.5 ) * ( ( v[id-st1] + v[id    ] ) * 0.5 )
                            )
                      + 1.0 / mesh.delta(0) *
                            (   ( ( u[id    ] + u[id+st1    ] ) * 0.5 ) * ( ( v[id    ] + v[id+st0] ) * 0.5 )
                              - ( ( u[id-st0] + u[id-st0+st1] ) * 0.5 ) * ( ( v[id-st0] + v[id    ] ) * 0.5 ) );



    return -adv;
}


real
calc_prs_u
(
    const real * const p,
    const int id,
    const Mesh & mesh
){
    const int st0 = mesh.ST(0);

    const real prs = ( p[id+st0] - p[id] ) / mesh.delta(0);

    return - prs;
}

real
calc_prs_v
(
    const real * const p,
    const int id,
    const Mesh & mesh
){
    const int st1 = mesh.ST(1);

    const real prs = ( p[id+st1] - p[id] ) / mesh.delta(1);

    return - prs;
}

real
calc_dif_u
(
    const real * const u,
    const real Re,
    const int id,
    const Mesh & mesh
){
    const int st0 = mesh.ST(0);
    const int st1 = mesh.ST(1);

    const real dif =   ( u[id-st0] - 2.0 * u[id] + u[id+st0] ) / ( mesh.delta(0) * mesh.delta(0) )
                     + ( u[id-st1] - 2.0 * u[id] + u[id+st1] ) / ( mesh.delta(1) * mesh.delta(1) );

    return 1.0 / Re * dif;
}

real
calc_dif_v
(
    const real * const v,
    const real Re,
    const int id,
    const Mesh & mesh
){
    const int st0 = mesh.ST(0);
    const int st1 = mesh.ST(1);

    const real dif =   ( v[id-st0] - 2.0 * v[id] + v[id+st0] ) / ( mesh.delta(0) * mesh.delta(0) )
                     + ( v[id-st1] - 2.0 * v[id] + v[id+st1] ) / ( mesh.delta(1) * mesh.delta(1) );

    return 1.0 / Re * dif;
}

/*
real
calc_dif
(
    const real * const f,
    const real mu,
    const int st0,
    const int st1,
    const real dx0,
    const real dx1
){
    const real dif =   ( difference( &f[0], st0, dx0 ) - difference( &f[-st0], st0, dx0 ) ) / dx0
                     + ( difference( &f[0], st1, dx1 ) - difference( &f[-st1], st1, dx1 ) ) / dx1;
    return mu * dif;
}

real
calc_prs
(
    const real * const p,
    const int stride,
    const real dx
){
    const real prs = ( p[stride] - p[0] ) / dx;
    return - prs;
}
*/
