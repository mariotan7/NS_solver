#pragma once
#include "real.h"
#include "Mesh.h"

void output_result
(
    const char * const u_filename,
    const char * const v_filename,
    const real * const u,
    const real * const v, 
    const Mesh & mesh
);

void calc_u_result
(
    real * const u_result,
    real * const y,
    const real * const u,
    const Mesh & mesh
);

void calc_v_result
(
    real * const v_result,
    real * const x,
    const real * const v,
    const Mesh & mesh
);

void save_result_file
(
    const char * const u_filename,
    const char * const v_filename,
    const real * const u,
    const real * const y,
    const real * const v,
    const real * const x,
    const Mesh & mesh
);


real average( const real a, const real b );
