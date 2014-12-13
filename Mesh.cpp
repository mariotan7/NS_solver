#include "Mesh.h"
//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
Mesh::Mesh
(
    const int  in_HN,
    const real in_L0, const int  in_N0,
    const real in_L1, const int  in_N1,
    const real in_L2, const int  in_N2
){
    m_HN = in_HN;

    m_L[0] = in_L0; m_L[1] = in_L1; m_L[2] = in_L2;
    m_N[0] = in_N0; m_N[1] = in_N1; m_N[2] = in_N2;

    /*
    if ( m_N[1] == 1 && m_N[2] == 1 ) {
        m_ST[0] = 1; m_ST[1] = 0     ; m_ST[2] = 0;
    } else if ( m_N[2] == 1 ) {
        m_ST[0] = 1; m_ST[1] = NN(0) ; m_ST[2] = 0;
    } else {
        m_ST[0] = 1; m_ST[1] = NN(0) ; m_ST[2] = NN(0)*NN(1);
    }
    */

    m_ST[0] = 1; m_ST[1] = NN(0) ; m_ST[2] = NN(0)*NN(1);

    f_array_b_id[0][0] = &Mesh::b_id_0_0; f_array_b_id[0][1] = &Mesh::b_id_0_1;
    f_array_b_id[1][0] = &Mesh::b_id_1_0; f_array_b_id[1][1] = &Mesh::b_id_1_1;
    f_array_b_id[2][0] = &Mesh::b_id_2_0; f_array_b_id[2][1] = &Mesh::b_id_2_1;

    f_array_b_N[0][0] = &Mesh::b_N_0_0; f_array_b_N[0][1] = &Mesh::b_N_0_1; f_array_b_N[0][2] = &Mesh::b_N_0_2;
    f_array_b_N[1][0] = &Mesh::b_N_1_0; f_array_b_N[1][1] = &Mesh::b_N_1_1; f_array_b_N[1][2] = &Mesh::b_N_1_2;
    f_array_b_N[2][0] = &Mesh::b_N_2_0; f_array_b_N[2][1] = &Mesh::b_N_2_1; f_array_b_N[2][2] = &Mesh::b_N_2_2;
}

//------------------------------------------------------------------
// for printing the information of this mesh
//------------------------------------------------------------------
void
Mesh::print_mesh_data( void ) const
{
    std::cout << "*--------------- Mesh information -------------***" << std::endl;
    std::cout << "*--\t  L  :\t" <<   L(0)  << "\t" <<   L(1) << "\t" <<   L(2) << std::endl;
    std::cout << "*--\t  N  :\t" <<   N(0)  << "\t" <<   N(1) << "\t" <<   N(2) << std::endl;
    std::cout << "*--\t HN  :\t" <<  HN( )  << "\t" <<  HN( ) << "\t" <<  HN( ) << std::endl;
    std::cout << "*--\t HNF :\t" << HNF(0)  << "\t" << HNF(1) << "\t" << HNF(2) << std::endl;
    std::cout << "*--\t ST  :\t" <<  ST(0)  << "\t" <<  ST(1) << "\t" <<  ST(2) << std::endl;
    std::cout << "*--\t NN  :\t" <<  NN(0)  << "\t" <<  NN(1) << "\t" <<  NN(2) << std::endl;
    std::cout << "*--\t size:\t" <<  size() << "\t" << std::endl;

    std::cout << std::endl;

    for ( int cf = 0 ; cf < 3 ; cf++ ) {
    for ( int ct = 0 ; ct < 3 ; ct++ ) {
        std::cout << "*--\t b_N(" << cf <<", "<< ct
            << ") = " <<  b_N(cf,ct) << "\t" << std::endl;
    } std::cout << std::endl;
    }

    std::cout << "*----------------------------------------------***" << std::endl;
}
