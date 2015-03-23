//
//  solve.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/10.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
#include <iostream>
#include <iomanip>
using namespace std;

void solver::solve(int totaliter)
{
    setBoundary();
    if (cRank==0)
    {
        //cout << timeIdx << endl;
        printstatus();
    }
    
    //RK4
    for (int iter=0; iter<3; ++iter)
    {
        if (cRank==0)
        {
            for (int iterCPU=1; iterCPU<numOfProcess; ++iterCPU)
            {
                MPI_Send(Fields->data, 1, BD4Type[iterCPU], iterCPU, iterCPU, MPI_COMM_WORLD);
            }
            MPI_Isend(Fields->data, 1, BD4Type[0], 0, 0, MPI_COMM_WORLD, &request);
            MPI_Irecv(HistoryFields[iter]->data, iterPoints, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
            
        }
        else
        {
            MPI_Recv(HistoryFields[iter]->data, iterPoints, MPI_DOUBLE, 0, cRank, MPI_COMM_WORLD, &status);
        }
        
        RK4Step();
        setBoundary();
        ++timeIdx;
//        if (cRank==0)
//        {
//            cout << timeIdx << endl;
//            printstatus();
//        }
    }
    
    //BDF4
    for (int iter=3; iter<totaliter; ++iter)
    {
        BDF4Step();
        setBoundary();
        ++timeIdx;
//        if (cRank==0)
//        {
//            cout << timeIdx << endl;
//            printstatus();
//        }
//        
        if (cRank==0 && timeIdx%262144==0)
        {
            timeIdx=timeIdx/262144;
            printstatus();
            timeIdx=timeIdx*262144;
        }
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
