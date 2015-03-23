//
//  solve.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/10.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
#include <iostream>

void solver::solve(double totaltime)
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
        printstatus();
    }
    
    //Increase time step
    for (int iter=0; iter<IncreaseTimes; ++iter)
    {
        DoubleTimeStep();
    }
    
    //BDF4
    while (time<totaltime)
    {
        BDF4Step();
        setBoundary();
        ++timeIdx;
        printstatus();
        //cout << timeIdx << endl;
        //        if (timeIdx%262144==0)
        //        {
        //            timeIdx=timeIdx/262144;
        //            printstatus();
        //            timeIdx=timeIdx*262144;
        //        }
    }
}

void solver::setBoundary()
{
    dr(1);
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
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

void solver::DoubleTimeStep()
{
    //stabilize
    for (int iter=0; iter<100; ++iter)
    {
        BDF4Step();
        setBoundary();
        ++timeIdx;
        printstatus();
    }
    for (int iterh=0; iterh<3; ++iterh)
    {
        gsl_matrix_memcpy(DoubledHistoryFields[iterh], Fields);
        for (int iter=0; iter<2; ++iter)
        {
            BDF4Step();
            setBoundary();
            ++timeIdx;
            printstatus();
        }
    }
    StepT*=2;
    gsl_matrix *swapTemp;
    for (int iter=0; iter<3; ++iter)
    {
        swapTemp=DoubledHistoryFields[iter];
        DoubledHistoryFields[iter]=HistoryFields[iter];
        HistoryFields[iter]=swapTemp;
    }
}