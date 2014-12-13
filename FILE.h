#ifndef FILE_H
#define FILE_H

#include <string>

#include "Mesh.h"

void save_array
(
    const int d,
    const real * const vof,
    const std::string filename,
    const std::string position,
    const Mesh &mesh
);

#endif //FILE_H
