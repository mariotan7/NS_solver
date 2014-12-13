#include "calc_values.h"
#include "NS.h"
#include <cmath>
#include <iostream>

void abs_array
(
    real * const array,
    const Mesh & mesh
){
    const int k = mesh.hn() ;
    for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
    for ( int i = mesh.hn() ; i < mesh.hnf(0) ; i++ ) {

        const int id = mesh.id(i,j,k);

        array[id] = std::abs( array[id] );

    }}
}

void calc_div
(
          real * const div,
    const real * const u,
    const real * const v,
    const Mesh & mesh
){
    const int k = mesh.hn() ;
    for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
    for ( int i = mesh.hn() ; i < mesh.hnf(0) ; i++ ) {

        const int id = mesh.id(i,j,k);

        div[id] = calc_div_locol( &u[id], &v[id], mesh.st(0), mesh.st(1), mesh.delta(0), mesh.delta(1));
    }}
}

real integration
(
    const real * const array,
    const Mesh & mesh
){
    real value = 0.0;

    const int k = mesh.hn() ;
    for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
    for ( int i = mesh.hn() ; i < mesh.hnf(0) ; i++ ) {

        const int id = mesh.id(i,j,k);
        value += array[id];
    }}

    const real vol = mesh.delta(0) * mesh.delta(1);

    return value * vol;
}

real max_from_array
(   
    const real * const array,
    const Mesh & mesh
){
    real max_value = 0.0;

    const int k = mesh.hn() ;
    for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
    for ( int i = mesh.hn() ; i < mesh.hnf(0) ; i++ ) {
        const int id = mesh.id(i,j,k);

        max_value = std::max( max_value, array[id] );

    }}

    return max_value;
}
