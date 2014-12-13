#ifndef MESH_H
#define MESH_H

#include <iostream>

typedef double real;

//-----------------------------------------------------
// Mesh mesh class
//-----------------------------------------------------
class Mesh {

public:

    //------------------------------------------------------------------
    // constructor & destructor
    //------------------------------------------------------------------
    Mesh
    (
        const int  in_HN,
        const real in_L0,       const int  in_N0,
        const real in_L1 = 0.0, const int  in_N1 = 1,
        const real in_L2 = 0.0, const int  in_N2 = 1
    );

    ~Mesh(void) {}

    //------------------------------------------------------------------
    // getter
    //------------------------------------------------------------------
    real L   ( const int d ) const { return l(d);    }
    real l   ( const int d ) const { return m_L[d];  }

    int  N   ( const int d ) const { return n(d);    }
    int  n   ( const int d ) const { return m_N[d];  }

    int  ST  ( const int d ) const { return st(d);   }
    int  st  ( const int d ) const { return m_ST[d]; }

    int b_id ( const int cf, const int s, const int bi, const int bj, const int bk ) const
                { return (this->*(f_array_b_id[cf][s]))(bi, bj, bk);  };

    int b_N  ( const int cf, const int ct ) const
                { return (this->*(f_array_b_N[cf][ct]))();  };

    //------------------------------------------------------------------
    // Calculater
    //------------------------------------------------------------------
    int  HNF    ( const int d ) const { return hnf(d);}
    int  NN     ( const int d ) const { return nn(d) ;}
    int  HN     ( void )        const { return hn()  ;}

    int  hnf    ( const int d ) const { return m_N[d] + 1 * m_HN;                    }
    int  nn     ( const int d ) const { return m_N[d] + 2 * m_HN;                    }
    int  hn     ( void )        const { return m_HN;                                 }

    int  size   ( void )        const { return NN(0) * NN(1) * NN(2);                }
    real delta  ( const int d ) const { return m_L[d] / static_cast<real>( m_N[d] ); }
    real dist_sq
    (
      const int i,         const int j,         const int k,
      const real px = 0.0, const real py = 0.0, const real pz = 0.0  ) const;

    //------------------------------------------------------------------
    // get caluculation region id
    //------------------------------------------------------------------
    int id ( const int i, const int j, const int k = 3 ) const { return m_ST[2] * k + m_ST[1] * j + m_ST[0] * i; }


    //------------------------------------------------------------------
    // return the cell center position and
    //------------------------------------------------------------------
    real x ( const int i, const int j, const int k ) const { return delta(0) * ( static_cast<real>(i - m_HN) + 0.5 ); }
    real y ( const int i, const int j, const int k ) const { return delta(1) * ( static_cast<real>(j - m_HN) + 0.5 ); }
    real z ( const int i, const int j, const int k ) const { return delta(2) * ( static_cast<real>(k - m_HN) + 0.5 ); }

    //------------------------------------------------------------------
    // for printing the information of this mesh
    //------------------------------------------------------------------
    void print_mesh_data ( void ) const;


private:

    //------------------------------------------------------------------
    // for printing the information of this Mesh
    //------------------------------------------------------------------
    real  m_L[3];
    int   m_N[3];
    int   m_ST[3];

    int   m_HN;

    //------------------------------------------------------------------
    // get boudary region id and width ( the array of funcitons )
    //------------------------------------------------------------------
    int (Mesh::* f_array_b_id[3][2]) ( const int bi, const int bj, const int bk ) const;
    int (Mesh::* f_array_b_N [3][3]) () const;

    //------------------------------------------------------------------
    // get boudary id functions
    //------------------------------------------------------------------
    int b_id_0_0 ( const int bi, const int bj, const int bk ) const;
    int b_id_0_1 ( const int bi, const int bj, const int bk ) const;

    int b_id_1_0 ( const int bi, const int bj, const int bk ) const;
    int b_id_1_1 ( const int bi, const int bj, const int bk ) const;

    int b_id_2_0 ( const int bi, const int bj, const int bk ) const;
    int b_id_2_1 ( const int bi, const int bj, const int bk ) const;

    //------------------------------------------------------------------
    // get boudary N
    //------------------------------------------------------------------
    int b_N_0_0 () const { return m_HN  ; }
    int b_N_0_1 () const { return m_N[1]; }
    int b_N_0_2 () const { return m_N[2]; }

    int b_N_1_0 () const { return m_N[0]; }
    int b_N_1_1 () const { return m_HN  ; }
    int b_N_1_2 () const { return m_N[2]; }

    int b_N_2_0 () const { return m_N[0]; }
    int b_N_2_1 () const { return m_N[1]; }
    int b_N_2_2 () const { return m_HN  ; }
};



//------------------------------------------------------------------
// Calculater
//------------------------------------------------------------------
inline real
Mesh::dist_sq
(
  const int i,   const int j,   const int k,
  const real px, const real py, const real pz ) const
{

    const real dx = px-x(i,j,k);
    const real dy = py-y(i,j,k);
    const real dz = pz-z(i,j,k);

    return dx*dx + dy*dy + dz*dz;
}


//------------------------------------------------------------------
// get boudary id functions
//------------------------------------------------------------------
inline
int
Mesh::b_id_0_0 ( const int bi, const int bj, const int bk ) const
{
    const int from = m_ST[2] * m_HN + m_ST[1] * m_HN;
    return    from + id(bi,bj,bk);
}

inline
int
Mesh::b_id_0_1 ( const int bi, const int bj, const int bk ) const
{
    const int from = m_ST[2] * m_HN + m_ST[1] * m_HN + HNF(0);
    return    from + id(bi,bj,bk);
}

inline
int
Mesh::b_id_1_0 ( const int bi, const int bj, const int bk ) const
{
    const int from = m_ST[2] * m_HN + m_HN;
    return    from + id(bi,bj,bk);
}

inline
int
Mesh::b_id_1_1 ( const int bi, const int bj, const int bk ) const
{
    const int from = m_ST[2] * m_HN + m_ST[1] * HNF(1) + m_HN;
    return    from + id(bi,bj,bk);
}

inline
int
Mesh::b_id_2_0 ( const int bi, const int bj, const int bk ) const
{
    const int from = m_ST[1] * m_HN + m_HN;
    return    from + id(bi,bj,bk);
}

inline
int
Mesh::b_id_2_1 ( const int bi, const int bj, const int bk ) const
{
    const int from = m_ST[2] * HNF(2) + m_ST[1] * m_HN + m_HN;
    return    from + id(bi,bj,bk);
}

#endif //MESH_H
