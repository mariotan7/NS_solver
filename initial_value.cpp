#include "initial_value.h"

void make_initial_value_circle
(
    real * const f,
    const Mesh &mesh
){
    for ( int k = mesh.HN() ; k < mesh.HNF(2) ; k++ ) {
    for ( int j = mesh.HN() ; j < mesh.HNF(1) ; j++ ) {
    for ( int i = mesh.HN() ; i < mesh.HNF(0) ; i++ ) {

        const int id = mesh.id(i,j,k);

        const real dist_sq = mesh.dist_sq( i, j, k, 0.5, 0.5, 0.5 );

        if ( dist_sq < 0.2*0.2 ) {
            f[id] = 1.0;
        } else {
            f[id] = 0.0;
        }

    }}}
}

void uniform_value
(
    real * const u,
    const real value,
    const Mesh &mesh
){
    for ( int id = 0 ; id < mesh.size() ; id++ ) {
        u[id] = value;
    }
}

void make_initial_value_1
(
    real * const u,
    const Mesh &mesh
){
    uniform_value( u, 1.0, mesh );
}

void make_initial_value_0
(
    real * const u,
    const Mesh &mesh
){
    uniform_value( u, 0.0, mesh );
}
