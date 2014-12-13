#include "FILE.h"

#include <iostream>
#include <fstream>
#include <sstream>

void save_array
(
    const int d,
    const real * const vof,
    const std::string filename,
    const std::string position,
    const Mesh &mesh
){
    std::ofstream file( filename );

    int k = mesh.HN();

    if ( d == 3 ) {
        k = mesh.NN(2)/2;
    }

    //for ( int i = mesh.HN() ; i < mesh.HNF(0) ; i++ ) {
    //for ( int j = mesh.HN() ; j < mesh.HNF(1) ; j++ ) {
    for ( int i = 0 ; i < mesh.NN(0) ; i++ ) {
        for ( int j = 0 ; j < mesh.NN(1) ; j++ ) {

            const int id = mesh.id(i,j,k);

            real x, y;
            if ( position == "cell_center" ) {
                x = mesh.x(i,j,k);
                y = mesh.y(i,j,k);
            } else if ( position == "staggered_x" ) {
                x = mesh.x(i,j,k) + 0.5 * mesh.delta(0);
                y = mesh.y(i,j,k);
            } else if ( position == "staggered_y" ) {
                x = mesh.x(i,j,k);
                y = mesh.y(i,j,k) + 0.5 * mesh.delta(1);
            } else {
                std::cout << "error at " << filename << std::endl;
                x = mesh.x(i,j,k);
                y = mesh.y(i,j,k);
            }

            file << x << "\t" << y << "\t" << vof[id] << std::endl;

        }

        file << std::endl;
    }

    std::cout << filename << " is saved." << std::endl;
}
