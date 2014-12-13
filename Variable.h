#ifndef VARIABLE_H
#define VARIABLE_H

#include "BoundaryCondition.h"
#include "Scheme.h"
#include "SumUpRHS.h"

//-----------------------------------------------
// Loop for cell center variables
//-----------------------------------------------
namespace CellCenter{

real *
calc_rhs
(
    const real * const * const a_f,
    const Scheme &scheme,
    const Mesh & mesh
){
    real * rhs = new real[mesh.size()];

    for ( int k = mesh.HN() ; k < mesh.HNF(2) ; k++ ) {
    for ( int j = mesh.HN() ; j < mesh.HNF(1) ; j++ ) {
    for ( int i = mesh.HN() ; i < mesh.HNF(0) ; i++ ) {

        const int id = mesh.id(i,j,k);
        rhs[id] = scheme( a_f, id );

    }}}

    return rhs;
}

template <int N>
real *
time_marching
(
    const real * const f,
    const SumUpRHS<N> sum_Up_RHS,
    const Mesh & mesh
){
    real * nf = new real[mesh.size()];

    for ( int k = mesh.HN() ; k < mesh.HNF(2) ; k++ ) {
    for ( int j = mesh.HN() ; j < mesh.HNF(1) ; j++ ) {
    for ( int i = mesh.HN() ; i < mesh.HNF(0) ; i++ ) {

        const int id = mesh.id(i,j,k);
        nf[id] = f[id] + sum_Up_RHS( id );

    }}}

    return nf;
}

} // namespace CellCenter

//-----------------------------------------------
// Loop for x cell face variables
//-----------------------------------------------
namespace CellFaceX{

real *
calc_rhs
(
    const real * const * const a_f,
    const Scheme &scheme,
    const Mesh & mesh
){
    real * rhs = new real[mesh.size()];

    for ( int k = mesh.HN() ; k < mesh.HNF(2)     ; k++ ) {
    for ( int j = mesh.HN() ; j < mesh.HNF(1)     ; j++ ) {
    for ( int i = mesh.HN() ; i < mesh.HNF(0) - 1 ; i++ ) {

        const int id = mesh.id(i,j,k);
        rhs[id] = scheme( a_f, id );

    }}}

    return rhs;
}

template <int N>
real *
time_marching
(
    const real * const f,
    const SumUpRHS<N> sum_Up_RHS,
    const Mesh & mesh
) {
    real * nf = new real[mesh.size()];

    for ( int k = mesh.HN() ; k < mesh.HNF(2)     ; k++ ) {
    for ( int j = mesh.HN() ; j < mesh.HNF(1)     ; j++ ) {
    for ( int i = mesh.HN() ; i < mesh.HNF(0) - 1 ; i++ ) {

        const int id = mesh.id(i,j,k);
        nf[id] = f[id] + sum_Up_RHS( id );

    }}}

    return nf;
}

template <class FUNCTOR>
void
for_each
(
    real * const f,
    const FUNCTOR func,
    const Mesh & mesh
){

    for ( int k = mesh.HN() ; k < mesh.HNF(2)     ; k++ ) {
    for ( int j = mesh.HN() ; j < mesh.HNF(1)     ; j++ ) {
    for ( int i = mesh.HN() ; i < mesh.HNF(0) - 1 ; i++ ) {

        const int id = mesh.id(i,j,k);
        f[id] = func( f, id );

    }}}
}

} // namespace CellFaceX

//-----------------------------------------------
// Loop for y cell face variables
//-----------------------------------------------
namespace CellFaceY{

real *
calc_rhs
(
    const real * const * const a_f,
    const Scheme &scheme,
    const Mesh & mesh
){
    real * rhs = new real[mesh.size()];

    for ( int k = mesh.HN() ; k < mesh.HNF(2)     ; k++ ) {
    for ( int j = mesh.HN() ; j < mesh.HNF(1) - 1 ; j++ ) {
    for ( int i = mesh.HN() ; i < mesh.HNF(0)     ; i++ ) {

        const int id = mesh.id(i,j,k);
        rhs[id] = scheme( a_f, id );

    }}}

    return rhs;
}

template <int N>
real *
time_marching
(
    const real * const f,
    const SumUpRHS<N> sum_Up_RHS,
    const Mesh & mesh
) {
    real * nf = new real[mesh.size()];

    for ( int k = mesh.HN() ; k < mesh.HNF(2)     ; k++ ) {
    for ( int j = mesh.HN() ; j < mesh.HNF(1) - 1 ; j++ ) {
    for ( int i = mesh.HN() ; i < mesh.HNF(0)     ; i++ ) {

        const int id = mesh.id(i,j,k);
        nf[id] = f[id] + sum_Up_RHS( id );

    }}}

    return nf;
}

template <class FUNCTOR>
void
for_each
(
    real * const f,
    const FUNCTOR func,
    const Mesh & mesh
){
    for ( int k = mesh.HN() ; k < mesh.HNF(2)     ; k++ ) {
    for ( int j = mesh.HN() ; j < mesh.HNF(1) - 1 ; j++ ) {
    for ( int i = mesh.HN() ; i < mesh.HNF(0)     ; i++ ) {

        const int id = mesh.id(i,j,k);
        f[id] = func( f, id );

    }}}
}

} // namespace CellFaceY

//-----------------------------------------------
// Loop for z cell face variables
//-----------------------------------------------
namespace CellFaceZ{

real *
calc_rhs
(
    const real * const * const a_f,
    const Scheme &scheme,
    const Mesh & mesh
){
    real * rhs = new real[mesh.size()];

    for ( int k = mesh.HN() ; k < mesh.HNF(2) - 1 ; k++ ) {
    for ( int j = mesh.HN() ; j < mesh.HNF(1)     ; j++ ) {
    for ( int i = mesh.HN() ; i < mesh.HNF(0)     ; i++ ) {

        const int id = mesh.id(i,j,k);
        rhs[id] = scheme( a_f, id );

    }}}

    return rhs;
}

template <int N>
real *
time_marching
(
    const real * const f,
    const SumUpRHS<N> sum_Up_RHS,
    const Mesh & mesh
) {
    real * nf = new real[mesh.size()];

    for ( int k = mesh.HN() ; k < mesh.HNF(2) - 1 ; k++ ) {
    for ( int j = mesh.HN() ; j < mesh.HNF(1)     ; j++ ) {
    for ( int i = mesh.HN() ; i < mesh.HNF(0)     ; i++ ) {

        const int id = mesh.id(i,j,k);
        nf[id] = f[id] + sum_Up_RHS( id );

    }}}

    return nf;
}

} // namespace CellFaceZ

//-----------------------------------------------
// boundary update
//-----------------------------------------------
namespace CellCenter{

void
part_boundary_update
(
    real * const f,
    const BCondition &bc,
    const int cf,
    const int s,
    const Mesh & mesh
){

    for ( int bk = 0 ; bk < mesh.b_N( cf, 2 ) ; bk++ ) {
    for ( int bj = 0 ; bj < mesh.b_N( cf, 1 ) ; bj++ ) {
    for ( int bi = 0 ; bi < mesh.b_N( cf, 0 ) ; bi++ ) {
        const int id = mesh.b_id( cf, s, bi, bj, bk );
        f[id] =  bc( f, id, bi, bj, bk );
    }}}
}

template<int D>
void
update_boundary
(
    real * const f,
    const BCondition * const * const bc,
    const Mesh & mesh
){

    for ( int cf = 0 ; cf < D ; cf++ ) {
    for ( int  s = 0 ;  s < 2 ;  s++ ) {
        const int cfs = cf * 2 + s;
        CellCenter::part_boundary_update( f, *(bc[cfs]), cf, s, mesh );
    }}
}

} // namespace CellCenter

#endif //VARIABLE_H
