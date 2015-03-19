//
//  solve.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/10.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
#include <iostream>
using namespace std;

void solver::solve(int totaliter)
{
    setBoundary();
    
    //RK4
    for (int iter=0; iter<3; ++iter)
    {
        if (0==cRank)
        {
            gsl_matrix_memcpy(HistoryFields[iter],Fields);
        }
        
        RK4Step();
        setBoundary();
        ++timeIdx;
//        if(0==cRank)
//            printstatus();
    }
    
    //BDF4
    for (int iter=3; iter<totaliter; ++iter)
    {
        BDF4Step();
        setBoundary();
        ++timeIdx;
        //printstatus();
//        if (cRank==0)
//        {
//            cout << timeIdx << endl;
//        }
        
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
    if(0==cRank)
    {
        for (int iterCPU=0; iterCPU<numOfProcessT; ++iterCPU)
        {
            MPI_Send(Fields->data, 1, TblockType[iterCPU], iterCPU*2+1, iterCPU*2+1, MPI_COMM_WORLD);
        }
    }
    if(cRank%2!=0)
    {
        MPI_Recv(dctr->data, jobPointsT, MPI_DOUBLE, 0, cRank, MPI_COMM_WORLD, &status);
        
        for (int iter=0; iter<jobT; ++iter)
        {
            tempFieldsLocal->data[iter]=dctr->data[iter];
        }
        dr(1);
        for (int iter=0; iter<jobT; ++iter)
        {
            double d0=gsl_matrix_get(dctr, 0, iter);
            d0/=-Nr*(Nr-2);
            d0*=3;
            tempFieldsLocal->data[iter]+=d0;
        }
        MPI_Send(tempFieldsLocal->data, jobT, MPI_DOUBLE, 0, 100+cRank, MPI_COMM_WORLD);
    }
    if (0==cRank)
    {
        for (int iterCPU=0; iterCPU<numOfProcessT; ++iterCPU)
        {
            MPI_Recv(Fields->data, 1, BoundaryType[iterCPU], iterCPU*2+1, 100+iterCPU*2+1, MPI_COMM_WORLD, &status);
        }
    }
}
