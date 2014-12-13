#ifndef STENCILS_H
#define STENCILS_H

template<int N, typename REAL>
class Stencils
{
public:

    //---------------------------
    // constructor & destructor
    //---------------------------
    Stencils( void ) : m_val( new REAL[N] ) {}
    virtual ~Stencils(void) { delete[] m_val; }

    //---------------------------
    // copy constructors
    //---------------------------
    Stencils( const Stencils<N, REAL> & sf )
    {
        m_val = new REAL[N];

        for ( int n = 0 ; n < N ; n++ ) {
            m_val[n] = sf(n);
        }
    }

    //---------------------------
    // operators
    //---------------------------
    const REAL operator () ( const int n ) const { return m_val[n]; }

    Stencils<N, REAL>   operator +  ( const Stencils<N, REAL> & rhs );
    Stencils<N, REAL>   operator -  ( const Stencils<N, REAL> & rhs );
    Stencils<N, REAL>   operator *  ( const Stencils<N, REAL> & rhs );
    Stencils<N, REAL>   operator /  ( const Stencils<N, REAL> & rhs );

    Stencils<N, REAL> & operator  = ( const Stencils<N, REAL> & rhs );
    Stencils<N, REAL> & operator += ( const Stencils<N, REAL> & rhs );
    Stencils<N, REAL> & operator -= ( const Stencils<N, REAL> & rhs );
    Stencils<N, REAL> & operator *= ( const Stencils<N, REAL> & rhs );
    Stencils<N, REAL> & operator /= ( const Stencils<N, REAL> & rhs );

protected:
    REAL & operator [] ( const int n ) { return m_val[n]; }

    //const int get_f_id( const int n, const int st ) const { return st * ( -(N-1)/2 + n ); }
    const int get_f_id( const int n, const int st ) const { return st * ( -(N)/2 + n ); }
    //const int get_f_id( const int n, const int st ) const { return st * ( n - N/2 ); }


private:
    REAL *m_val;
};

template<int N, typename REAL>
class Interpolations_1st : public Stencils<N, REAL>
{
public:

    typedef Stencils<N, REAL> BASE;

    //-----------------------
    // constructors
    //-----------------------
    Interpolations_1st( const REAL * const f, const int stride ) : Stencils<N, REAL>()
    {
        for ( int n = 0 ; n < N ; n++ ) {
            const int f_id = BASE::get_f_id( n, stride );
            (*this)[n] = ( f[ f_id + stride ] + f[f_id] ) * 0.5;
            //std::cout << f_id + stride << "\t" << f_id << std::endl;
            //abort();
        }

    }

    Interpolations_1st( const Stencils<N+1, REAL> & sf ) : Stencils<N, REAL>()
    {
        for ( int n = 0 ; n < N ; n++ ) {
            (*this)[n] = ( sf(n+1) + sf(n) ) * 0.5;
        }
    }

};

template<int N, typename REAL>
class Differences_1st : public Stencils<N, REAL>
{
public:

    typedef Stencils<N, REAL> BASE;

    Differences_1st( const REAL * const f, const int stride, const REAL dx ) : Stencils<N, REAL>()
    {
        //const REAL inv_dx = 1.0 / dx;

        for ( int n = 0 ; n < N ; n++ ) {
            const int f_id = BASE::get_f_id( n, stride );
            //(*this)[n] = ( f[ f_id + stride ] - f[f_id] ) * inv_dx;
            (*this)[n] = ( f[ f_id + stride ] - f[f_id] ) / dx;
        }
    }

    Differences_1st( const Stencils<N+1, REAL> & sf, const REAL dx ) : Stencils<N, REAL>()
    {
        //const REAL inv_dx = 1.0 / dx;

        for ( int n = 0 ; n < N ; n++ ) {
            //(*this)[n] = ( sf(n+1) - sf(n) ) * inv_dx;
            (*this)[n] = ( sf(n+1) - sf(n) ) / dx;
        }
    }
};


//-----------------------------------------------------------
// Stencils<N, REAL> class operators
//-----------------------------------------------------------
template<int N, typename REAL>
inline
Stencils<N, REAL>
Stencils<N, REAL>::operator + (const Stencils<N, REAL> & rhs )
{
    Stencils<N, REAL> sum;
    for ( int n = 0 ; n < N ; n++ ) {
        sum.m_val[n] = m_val[n] + rhs[n];
    }
    return sum;
}

template<int N, typename REAL>
inline
Stencils<N, REAL>
Stencils<N, REAL>::operator - (const Stencils<N, REAL> & rhs )
{
    Stencils<N, REAL> product;
    for ( int n = 0 ; n < N ; n++ ) {
        product.m_val[n] = m_val[n] - rhs[n];
    }
    return product;
}

template<int N, typename REAL>
inline
Stencils<N, REAL>
Stencils<N, REAL>::operator * (const Stencils<N, REAL> & rhs )
{
    Stencils<N, REAL> product;
    for ( int n = 0 ; n < N ; n++ ) {
        product.m_val[n] = m_val[n] * rhs(n);
    }
    return product;
}

template<int N, typename REAL>
inline
Stencils<N, REAL>
Stencils<N, REAL>::operator / (const Stencils<N, REAL> & rhs )
{
    Stencils<N, REAL> product;
    for ( int n = 0 ; n < N ; n++ ) {
        product.m_val[n] = m_val[n] / rhs[n];
    }
    return product;
}

template<int N, typename REAL>
inline
Stencils<N, REAL> &
Stencils<N, REAL>::operator = (const Stencils<N, REAL> & rhs )
{
    for ( int n = 0 ; n < N ; n++ ) {
        m_val[n] = rhs[n];
    }
    return *this;
}

template<int N, typename REAL>
inline
Stencils<N, REAL> &
Stencils<N, REAL>::operator += (const Stencils<N, REAL> & rhs )
{
    for ( int n = 0 ; n < N ; n++ ) {
        m_val[n] += rhs[n];
    }
    return *this;
}

template<int N, typename REAL>
inline
Stencils<N, REAL> &
Stencils<N, REAL>::operator -= (const Stencils<N, REAL> & rhs )
{
    for ( int n = 0 ; n < N ; n++ ) {
        m_val[n] -= rhs[n];
    }
    return *this;
}

template<int N, typename REAL>
inline
Stencils<N, REAL> &
Stencils<N, REAL>::operator *= (const Stencils<N, REAL> & rhs )
{
    for ( int n = 0 ; n < N ; n++ ) {
        m_val[n] *= rhs[n];
    }
    return *this;
}

template<int N, typename REAL>
inline
Stencils<N, REAL> &
Stencils<N, REAL>::operator /= (const Stencils<N, REAL> & rhs )
{
    for ( int n = 0 ; n < N ; n++ ) {
        m_val[n] /= rhs[n];
    }
    return *this;
}

#endif //STENCILS_H
