//
//  initialization.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/11.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
#include <cmath>

void solver::initialization()
{
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        for (int itert=0; itert<Ntheta; ++itert)
        {
            for (int iterr=0; iterr<Nrp; ++iterr)
            {
                //gsl_matrix_set(Fields, iterf*Nrp+iterr, itert, exp(-10*r->data[iterr]*r->data[iterr]));
                gsl_matrix_set(Fields, iterf*Nrp+iterr, itert, sin(theta->data[itert]));
            }
        }
    }
}
