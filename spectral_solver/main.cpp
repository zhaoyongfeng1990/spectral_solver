//
//  main.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include <iostream>
#include "solver.h"
#include <ctime>
using namespace std;



int main(int argc, const char  *argv[])
{
    MPI_Init(NULL, NULL);
    time_t ctime1, ctime2;
    time(&ctime1);
    solver* test1=new solver();
    test1->initialization();
    //if (test1->cRank==0)
    //{
    //    test1->readFile("0.txt");
        //test1->printstatus();
    //}
    //test1->setBoundary();
    
    //test1->timeIdx+=100;
    if (test1->cRank==0)
    {
        //test1->printstatus();
    }
    //test1->readFile("48.txt");
    //test1->timeIdx=16384*128*48;
    //test1->solve(16384*16);
    test1->solve(4000);
    //test1->RK4Step();
    
    //test1->setBoundary();
    //++test1->timeIdx;
    
    //test1->Fun(test1->Fields);
    
//    if(0==test1->cRank)
//    {
//        for (int iterCPU=0; iterCPU<test1->numOfProcessT; ++iterCPU)
//        {
//            MPI_Ssend(test1->Fields.data, 1, test1->TblockType[iterCPU], iterCPU*2+1, iterCPU*2+1, MPI_COMM_WORLD);
//        }
//    }
//    if(test1->cRank%2!=0)
//    {
//        MPI_Recv(test1->dctr.data, test1->jobPointsT, MPI_LONG_DOUBLE, 0, test1->cRank, MPI_COMM_WORLD,& test1->status);
//        //test1->dr(1);
//        test1->drWOA();
//        MPI_Ssend(test1->dctr.data, test1->jobPointsT, MPI_LONG_DOUBLE, 0, 100+test1->cRank, MPI_COMM_WORLD);
//    }
//    if (0==test1->cRank)
//    {
//        for (int iterCPU=0; iterCPU<test1->numOfProcessT; ++iterCPU)
//        {
//            MPI_Recv(test1->Fields.data, 1, test1->TblockType[iterCPU], iterCPU*2+1, iterCPU*2+101, MPI_COMM_WORLD,& test1->status);
//        }
//        test1->Fields*=(1.0/radius);
//    }
    
//    if (0==test1->cRank)
//    {
//        for (int iterCPU=1; iterCPU<test1->numOfProcessR; ++iterCPU)
//        {
//            MPI_Ssend(test1->Fields.data, 1, test1->RblockType[iterCPU], 2*iterCPU, 2*iterCPU, MPI_COMM_WORLD);
//        }
//        for (int iterf=0; iterf<NumField; ++iterf)
//        {
//            for (int iter=0; iter<test1->jobPointsRl; ++iter)
//            {
//                test1->FieldsLocal.data[iter+iterf*test1->jobPointsRl]=test1->Fields.data[iter+iterf*NumPoints];
//            }
//        }
//    }
//    
//    if (test1->cRank%2==0)
//    {
//        if (test1->cRank!=0)
//            MPI_Recv(test1->FieldsLocal.data, test1->jobPointsR, MPI_LONG_DOUBLE, 0, test1->cRank, MPI_COMM_WORLD,& test1->status);
//        test1->dtheta(1);
//        if (test1->cRank!=0)
//        {
//            MPI_Ssend(test1->dFieldsLocal.data, test1->jobPointsR, MPI_LONG_DOUBLE, 0, 100+test1->cRank, MPI_COMM_WORLD);
//        }
//    }
//    if (0==test1->cRank)
//    {
//        for (int iterf=0; iterf<NumField; ++iterf)
//        {
//            for (int iter=0; iter<test1->jobPointsRl; ++iter)
//            {
//                test1->Fields.data[iter+iterf*NumPoints]=test1->dFieldsLocal.data[iter+iterf*test1->jobPointsRl];
//            }
//        }
//        for (int iterCPU=1; iterCPU<test1->numOfProcessR; ++iterCPU)
//        {
//            MPI_Recv(test1->Fields.data, 1, test1->RblockType[iterCPU], iterCPU*2, 100+iterCPU*2, MPI_COMM_WORLD,& test1->status);
//        }
//    }
    
    //test1->timeIdx+=100;
    //test1->Fun(test1->Fields);
    //test1->setBoundary();
    time(&ctime2);
    if (test1->cRank==0)
    {
        test1->printstatus();
        cout << ctime2-ctime1 << endl;
    }
    
    delete test1;
    
    MPI_Finalize();
    return 0;
}
