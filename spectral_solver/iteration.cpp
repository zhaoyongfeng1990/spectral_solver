//
//  solve.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/10.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"

void solver::solve(int totaliter)
{
    setBoundary();
    //printstatus();
    //RK4
    
    for (int iter=0; iter<3; ++iter)
    {
        gsl_matrix_memcpy(HistoryFields[iter],Fields);
        RK4Step();
        setBoundary();
        ++timeIdx;
        //printstatus();
    }
    
    //BDF4
    for (int iter=3; iter<totaliter; ++iter)
    {
        BDF4Step();
        setBoundary();
        ++timeIdx;
        //printstatus();
    }
    printstatus();
    
}

void solver::setBoundary()
{
    dr(1);
#pragma omp parallel for
    for (int iter=0; iter<Ntheta; ++iter)
    {
        for (int iterField=0; iterField<NumField; ++iterField)
        {
            double d0=gsl_matrix_get(dFields, iterField*Nrp, iter);
            d0/=-Nr*(Nr-2);
            d0*=3;
            d0+=gsl_matrix_get(Fields, iterField*Nrp, iter);
            gsl_matrix_set(Fields, iterField*Nrp, iter, d0);
        }
    }
}