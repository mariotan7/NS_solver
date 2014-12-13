#pragma once

#include "Mesh.h"

void make_initial_value_circle
(
    real * const f,
    const Mesh &mesh
);

void make_initial_value_1
(
    real * const u,
    const Mesh &mesh
);

void make_initial_value_0
(
    real * const u,
    const Mesh &mesh
);

void uniform_value
(
    real * const u,
    const real value,
    const Mesh &mesh
);
