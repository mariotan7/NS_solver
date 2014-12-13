#ifndef NS_H
#define NS_H

#include "Mesh.h"
#include "Iter.h"

inline
real intrp
(
    const real * const f,
    const int stride
){
    return ( f[stride] + f[0] ) * 0.5;
}

inline
real difference
(
    const real * const f,
    const int stride,
    const real dx
){
    return ( f[stride] - f[0] ) / dx;
}

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
);

void
time_marching_RK2_cavity
(
    real * * const pu,
    real * * const pv,
    real * * const pp,
    const real dt,
    const Mesh & mesh
);

void copy_array
(
    real * const       cpy_u,
    const real * const u,
    const Mesh & mesh
);

void check_array
(
    const real * const a,
    const real * const b,
    const Mesh & mesh
);

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
);

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
);

//------------------------
// calc adv
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
);

real
calc_adv_u
(
    const real * const u,
    const real * const v,
    const int id,
    const Mesh & mesh
);

real
calc_adv_v
(
    const real * const u,
    const real * const v,
    const int id,
    const Mesh & mesh
);

//------------------------
// calc phi
//------------------------
void
calc_phi
(
          real * const phi,
    const real * const u,
    const real * const v,
    const real dt,
    const Mesh & mesh
);

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
);

real
calc_phi_SOR
(
          real * const phi,
    const real * const u,
    const real * const v,
    const real dt,
    const real beta,
    const Mesh & mesh
);

real calc_norm
(
    const real value,
    const int  n
);

real calc_B
(
    const real dx,
    const int  i_e,
    const int pi_e,
    const int  i_w,
    const int pi_w
);

real psi_ij
(
    const real * const u,
    const real * const v,
    const int id,
    const real dt,
    const Mesh & mesh
);

void
boundary_update_phi
(
    real * const phi,
    const Mesh & mesh
);

real calc_div_locol
(
    const real * const u,
    const real * const v,
    const int st0,
    const int st1,
    const real d0,
    const real d1
);

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
);

void
correct_u
(
    real * const u,
    const real * const phi,
    const real dt,
    const Mesh & mesh
);

void
correct_v
(
    real * const v,
    const real * const phi,
    const real dt,
    const Mesh & mesh
);

void
correct_p
(
    real * const p,
    const real * const phi,
    const Mesh & mesh
);


//------------------------
// boundary update u, v, p
//------------------------
void
boundary_update_u_v
(
    real * const u,
    real * const v,
    const Mesh &mesh
);

void
boundary_update_only_u
(
    real * const u,
    const Mesh & mesh
);

void
boundary_update_only_v
(
    real * const v,
    const Mesh & mesh
);

/*
real
calc_prs
(
    const real * const p,
    const int stride,
    const real dx
);

real
calc_dif
(
    const real * const f,
    const real mu,
    const int st0,
    const int st1,
    const real dx0,
    const real dx1
);

*/

real
calc_dif_u
(
    const real * const u,
    const real Re,
    const int id,
    const Mesh & mesh
);

real
calc_dif_v
(
    const real * const v,
    const real Re,
    const int id,
    const Mesh & mesh
);

real
calc_prs_u
(
    const real * const p,
    const int id,
    const Mesh & mesh
);

real
calc_prs_v
(
    const real * const p,
    const int id,
    const Mesh & mesh
);


#endif // NS_H
