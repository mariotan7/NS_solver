#include "output.h"
#include <cstdio>
#include <iostream>

void output_result
(
    const char * const u_filename,
    const char * const v_filename,
    const real * const u,
    const real * const v, 
    const Mesh & mesh
){
    real * u_result = new real[mesh.n(0)];
    real * v_result = new real[mesh.n(1)];

    real * y = new real[mesh.n(0)];
    real * x = new real[mesh.n(1)];

    calc_u_result( u_result, y, u, mesh );
    calc_v_result( v_result, x, v, mesh );

    save_result_file( u_filename, v_filename, u_result, y,
                                              v_result, x, mesh );

    delete[] u_result; delete[] y;
    delete[] v_result; delete[] x;
}


void calc_u_result
(
    real * const u_result,
    real * const y,
    const real * const u,
    const Mesh & mesh
){
    for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) {
        const int k = mesh.hn();

        const int i_e = mesh.nn(0) / 2 - 1;
        const int i_w = mesh.nn(0) / 2;

        const int j_shift = j - mesh.hn();

        const int id_e = mesh.id(i_e,j,k);
        const int id_w = mesh.id(i_w,j,k);

        u_result[j_shift] = average(u[id_e], u[id_w]);
               y[j_shift] = mesh.y(i_w,j,k);
    }
}

void calc_v_result
(
    real * const v_result,
    real * const x,
    const real * const v,
    const Mesh & mesh
){
    for ( int i = mesh.hn() ; i < mesh.hnf(0) ; i++ ) {
        const int k = mesh.hn();

        const int j_s = mesh.nn(1) / 2 - 1;
        const int j_n = mesh.nn(1) / 2;

        const int i_shift = i - mesh.hn();

        const int id_s = mesh.id(i,j_s,k);
        const int id_n = mesh.id(i,j_n,k);

        v_result[i_shift] = average(v[id_s], v[id_n]);
               x[i_shift] = mesh.x(i,j_n,k);
    }
}


void save_result_file
(
    const char * const u_filename,
    const char * const v_filename,
    const real * const u,
    const real * const y,
    const real * const v,
    const real * const x,
    const Mesh & mesh
){
    const real U = 1.0;
    const real H = 1.0;

    FILE * u_fp = fopen( u_filename, "w" );
    for ( int i = 0 ; i < mesh.n(1) ; i++ ) {

        const real uu = u[i] / U;
        const real yy = (y[i] - 0.5) / H;

        fprintf( u_fp, "%e\t%e\n", uu, yy );
    }
    fclose(u_fp);
    std::cout << u_filename << " is saved." << std::endl;

    FILE * v_fp = fopen( v_filename, "w" );
    for ( int j = 0 ; j < mesh.n(0) ; j++ ) {

        const real xx = (x[j] - 0.5) / H;
        const real vv = v[j] / U;

        fprintf( v_fp, "%e\t%e\n", xx, vv );
    }
    fclose(v_fp);
    std::cout << v_filename << " is saved." << std::endl;
}

real average( const real a, const real b )
{
    return (a+b) / 2.0;
}

/*
void calc_result
(
    real * const u,
    real * const y,
    const real * const u_reslt,
    const Mesh & mesh,
    const int d0, // u -> 1, v -> 0
    const int d1  // u -> 0, v -> 1
){
    for ( int j = mesh.hn() ; j < mesh.n(d0) ; j++ ) {
        const int k = mesh.nn(2) / 2;

        const int i_e = mesh.nn(d1) / 2 - 1;
        const int i_w = mesh.nn(d1) / 2;

        const int j_shift = j - mesh.hn();

        u_result[j_shift] = average(u[i_e], u[i_w]);
    }
}
*/
