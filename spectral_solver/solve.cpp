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
    printstatus();
    //RK4
    k1=gsl_matrix_alloc(matrixH, Ntheta);
    k2=gsl_matrix_alloc(matrixH, Ntheta);
    k3=gsl_matrix_alloc(matrixH, Ntheta);
    k4=gsl_matrix_alloc(matrixH, Ntheta);
    
    odetempField=gsl_matrix_alloc(matrixH, Ntheta);
    
    for (int iter=0; iter<3; ++iter)
    {
        gsl_matrix_memcpy(HistoryFields[iter],Fields);
        RK4Step();
        setBoundary();
        ++timeIdx;
        printstatus();
    }
    
    gsl_matrix_free(k1);
    gsl_matrix_free(k2);
    gsl_matrix_free(k3);
    gsl_matrix_free(k4);
    
    //BDF4
    for (int iter=3; iter<totaliter; ++iter)
    {
        BDF4Step();
        setBoundary();
        ++timeIdx;
        printstatus();
    }
    
    gsl_matrix_free(odetempField);
}

void solver::setBoundary()
{
    dr(1);
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
    
    //debug
    dr(1);
    printdebugM(dFields, "dr.txt");
}
