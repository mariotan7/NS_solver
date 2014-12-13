#pragma once
#include "real.h"

class Iter
{
public:
    Iter( const int iter, const real error )
        :iter(iter), error(error) {}

    int  get_iter()  const { return iter;  }
    real get_error() const { return error; }

private:
    const  int iter;
    const real error;
};
