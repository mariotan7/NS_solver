#include "NS.h"
#include "initial_value.h"
#include "FILE.h"
#include "output.h"
#include "calc_values.h"

void swap( real * * a, real * * b )
{
    real * temp = *a;
    *a = *b;
    *b = temp;
}


int test(int argc, char* argv[])
{
    if ( argc != 3 ) {
        std::cout << "error."    << std::endl;
        std::cout << "./run <n> <endtimestep>" << std::endl;
        return 1;
    }

    const int n = std::atoi(argv[1]); 
    
    //----------------------------------------------
    // define space domain
    //----------------------------------------------
    const int  D = 2;
    const real l = 1.0;
    const int hn = 3;
    //const int hn = 2;

    //----------------------------------------------
    // define mesh
    //----------------------------------------------
    const Mesh mesh( hn, l, n, l, n );
    //mesh.print_mesh_data();

    //----------------------------------------------
    // define time domain
    //----------------------------------------------
    const real      dt = 0.1 * mesh.delta(0) / 1.0; 
    //const real dt = 1.0e-5;
    const real endtime = 200 * 1.0 / 1.0;
    const int  endstep = static_cast<int>(endtime / dt);
    //const int  endstep = std::atoi(argv[2]);

    real * u = new real[mesh.size()];
    real * v = new real[mesh.size()];
    real * p = new real[mesh.size()];

    real * nu = new real[mesh.size()];
    real * nv = new real[mesh.size()];

    uniform_value(  u, 0.0, mesh );
    uniform_value(  v, 0.0, mesh );
    uniform_value( nu, 0.0, mesh );
    uniform_value( nv, 0.0, mesh );
    uniform_value(  p, 0.0, mesh );

    boundary_update_u_v( u, v, mesh );

    for ( int timestep = 0 ; timestep < endstep ; timestep++ ) {

        const Iter iter = time_marching_RK1_cavity( nu, nv, p, u, v, dt, mesh );

        swap( &nu, &u );
        swap( &nv, &v );
        //const Iter iter = time_marching_RK2_cavity( &u, &v, &p, dt, mesh );

        if ( timestep % 100 == 0 ) {

            real * const div = new real[mesh.size()];
            uniform_value( div, 0.0, mesh );

            calc_div( div, u, v, mesh );

            const std::string div_filename( "div.dat" );
            const std::string div_position( "cell_center" );
            //save_array( D, div, div_filename, div_position, mesh );

            abs_array( div, mesh );
            const real max_abs_div_value = max_from_array( div, mesh );
            delete[] div;
            
            std::cout << timestep / static_cast<real>(endstep) * 100 << " %";
            std::cout << "\t timestep = " << timestep;
            std::cout << "\titer = " << iter.get_iter() << "\terr = " << iter.get_error();
            std::cout << "\tdiv = " << std::scientific << std::showpos << max_abs_div_value << std::endl;
        }
    }

    const std::string u_filename( "u.dat" );
    const std::string u_position( "staggered_x" );
    save_array( D, u, u_filename, u_position, mesh );

    const std::string v_filename( "v.dat" );
    const std::string v_position( "staggered_y" );
    save_array( D, v, v_filename, v_position, mesh );

    const std::string p_filename( "p.dat" );
    const std::string p_position( "cell_center" );
    save_array( D, p, p_filename, p_position, mesh );

    output_result( "result_u.dat", "result_v.dat", u, v, mesh );

    delete[] u;
    delete[] v;
    delete[] p;
    delete[] nu;
    delete[] nv;

    return 0;
}

int test1(int argc, char* argv[])
{
    if ( argc != 2 ) {
        std::cout << "error."    << std::endl;
        std::cout << "./run <n>" << std::endl;
        return 1;
    }

    const int n = std::atoi(argv[1]); 
    
    //----------------------------------------------
    // define space domain
    //----------------------------------------------
    const int  D = 2;
    const real l = 1.0;
    const int hn = 3;

    //----------------------------------------------
    // define mesh
    //----------------------------------------------
    const Mesh mesh( hn, l, n, l, n );
    mesh.print_mesh_data();

    real * u = new real[mesh.size()];

    uniform_value( u, -100.0, mesh );
    const std::string u_filename( "u.dat" );
    const std::string u_position( "cell_center" );
    save_array( D, u, u_filename, u_position, mesh );

    delete[] u;

    return 0;
}

int test2(int argc, char * argv[])
{
    if ( argc != 2 ) {
        std::cout << "error."    << std::endl;
        std::cout << "./run <n>" << std::endl;
        return 1;
    }

    const int n = std::atoi(argv[1]); 
    
    //----------------------------------------------
    // define space domain
    //----------------------------------------------
    const int  D = 2;
    const real l = 1.0;
    const int hn = 3;

    //----------------------------------------------
    // define mesh
    //----------------------------------------------
    const Mesh mesh( hn, l, n, l, n );
    //mesh.print_mesh_data();

    real * phi = new real[mesh.size()];
    uniform_value( phi, 200.0, mesh );

    const int k = mesh.hn() ;
    for ( int j = mesh.hn() ; j < mesh.hnf(1) ; j++ ) { 
    for ( int i = mesh.hn() ; i < mesh.hnf(0) ; i++ ) {
        const int id = mesh.id(i,j,k);

        phi[id] = 100;
    }}

    const std::string u_init_filename( "phi_init.dat" );
    const std::string u_position( "cell_center" );
    save_array( D, phi, u_init_filename, u_position, mesh );

    boundary_update_phi( phi, mesh );

    const std::string u_filename( "phi.dat" );
    save_array( D, phi, u_filename, u_position, mesh );

    delete[] phi;
    return 0;
}

int test3(int argc, char * argv[])
{
    if ( argc != 2 ) {
        std::cout << "error."    << std::endl;
        std::cout << "./run <n>" << std::endl;
        return 1;
    }

    const int n = std::atoi(argv[1]); 
    
    //----------------------------------------------
    // define space domain
    //----------------------------------------------
    const int  D = 2;
    const real l = 1.0;
    const int hn = 3;

    //----------------------------------------------
    // define mesh
    //----------------------------------------------
    const Mesh mesh( hn, l, n, l, n );

    real * u = new real[mesh.size()];
    real * v = new real[mesh.size()];
    real * p = new real[mesh.size()];

    uniform_value( u, 2.0, mesh );
    uniform_value( v, 2.0, mesh );
    uniform_value( p, 0.0, mesh );

    const std::string u_init_filename( "u_init.dat" );
    const std::string u_position( "staggered_x" );
    save_array( D, u, u_init_filename, u_position, mesh );

    const std::string v_filename( "v.dat" );
    const std::string v_position( "staggered_y" );
    save_array( D, v, v_filename, v_position, mesh );

    const std::string p_filename( "p.dat" );
    const std::string p_position( "cell_center" );
    save_array( D, p, p_filename, p_position, mesh );

    //boundary_update_u_v( u, v, mesh );
    boundary_update_only_u(u, mesh);
    boundary_update_only_v(v, mesh);

    const std::string u_filename( "u.dat" );
    save_array( D, u, u_filename, u_position, mesh );

    const std::string v_init_filename( "v_init.dat" );
    save_array( D, v, v_init_filename, v_position, mesh );



    delete[] u;
    delete[] v;
    delete[] p;

    return 0;
}

int test5()
{
    real * a = new real[3];
    real * b = new real[3];

    a[0] = 1;
    a[1] = 1;
    a[2] = 1;

    b[0] = 2;
    b[1] = 2;
    b[2] = 2;

    for ( int i = 0 ; i < 3 ; i++ ) {
        std::cout << a[i] << "\t" << b[i] << std::endl;
    }

    swap(&a,&b);

    for ( int i = 0 ; i < 3 ; i++ ) {
        std::cout << a[i] << "\t" << b[i] << std::endl;
    }

    delete[] a, delete[] b;

    return 0;
}

int main(int argc, char* argv[])
{
    return test(argc, argv);
    //return test1(argc, argv);
    //return test2(argc, argv);
    //return test3(argc, argv);
    //return test5();
}
