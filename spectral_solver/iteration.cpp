//
//  solve.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/10.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
#include <iostream>

void solver::solve(int totaliter)
{
    setBoundary();
    printstatus();
    //RK4
    
    for (int iter=0; iter<3; ++iter)
    {
        *HistoryFields[iter]=Fields;
        RK4Step();
        //setBoundary();
        ++timeIdx;
//        printstatus();
    }
    
    //Increase time step
    //for (int iter=0; iter<IncreaseTimes; ++iter)
    //{
    //    DoubleTimeStepBDF6();
    //}
    
    //BDF
    for (int iter=3; iter<totaliter; ++iter)
    {
        BDF4Step();
        //RK6Step();
        //setBoundary();
        ++timeIdx;
//        printstatus();
        //cout << timeIdx << endl;
        if (timeIdx%4194304==0)
        {
            timeIdx=timeIdx/4194304;
            printstatus();
            timeIdx=timeIdx*4194304;
        }
    }
}

void solver::setBoundary()
{
    //drWOA();
    dr(1);
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iterField=0; iterField<NumField; ++iterField)
        {
            for (int iter=0; iter<Ntheta; ++iter)
            {
            long double d0=dFields.ele(iterField*Nrp, iter);
            d0/=-Nr*(Nr-2);
            d0*=3;
            Fields.ele(iterField*Nrp, iter)+=d0;
        }
    }
}

void solver::DoubleTimeStepBDF4()
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
        *DoubledHistoryFields[iterh]=Fields;
        for (int iter=0; iter<2; ++iter)
        {
            BDF4Step();
            setBoundary();
            ++timeIdx;
            printstatus();
        }
    }
    StepT*=2;
    long double* swapTemp;
    for (int iter=0; iter<3; ++iter)
    {
        swapTemp=DoubledHistoryFields[iter]->data;
        DoubledHistoryFields[iter]->data=HistoryFields[iter]->data;
        HistoryFields[iter]->data=swapTemp;
    }
}

void solver::DoubleTimeStepBDF6()
{
    //stabilize
    for (int iter=0; iter<100; ++iter)
    {
        BDF6Step();
        setBoundary();
        ++timeIdx;
        printstatus();
    }
    for (int iterh=0; iterh<5; ++iterh)
    {
        *DoubledHistoryFields[iterh]=Fields;
        for (int iter=0; iter<2; ++iter)
        {
            BDF6Step();
            setBoundary();
            ++timeIdx;
            printstatus();
        }
    }
    StepT*=2;
    long double* swapTemp;
    for (int iter=0; iter<5; ++iter)
    {
        swapTemp=DoubledHistoryFields[iter]->data;
        DoubledHistoryFields[iter]->data=HistoryFields[iter]->data;
        HistoryFields[iter]->data=swapTemp;
    }
}
