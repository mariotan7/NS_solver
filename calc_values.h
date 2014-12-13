#pragma once
#include "Mesh.h"
#include "real.h"

void calc_div
(
          real * const div,
    const real * const u,
    const real * const v,
    const Mesh & mesh
);

void abs_array
(
    real * const array,
    const Mesh & mesh
);

real integration
(
    const real * const array,
    const Mesh & mesh
);

real max_from_array
(   
    const real * const array,
    const Mesh & mesh
);
