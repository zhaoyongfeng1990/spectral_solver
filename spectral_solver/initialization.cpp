//
//  initialization.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/11.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"

void solver::initialization()
{
    
//#ifdef MULTIPROCESS
//    omp_set_num_threads(8);
//#endif
    
//#ifdef MULTIPROCESS
//#pragma omp parallel for
//#endif
    for (int itert=0; itert<Ntheta; ++itert)
    {
            for (int iterr=0; iterr<Nrp; ++iterr)
            {
//                gsl_matrix_set(Fields, iterr, itert, 2*exp(-r->data[iterr]*r->data[iterr]/(0.04/radius/radius)));
//                gsl_matrix_set(Fields, 2*Nrp+iterr, itert, 15);
                gsl_matrix_set(Fields, iterr, itert, exp(-10*r->data[iterr]*r->data[iterr]));
                gsl_matrix_set(Fields, Nrp+iterr, itert, exp(-10*r->data[iterr]*r->data[iterr]));
                gsl_matrix_set(Fields, 2*Nrp+iterr, itert, exp(-10*r->data[iterr]*r->data[iterr]));
            }
    }
}
