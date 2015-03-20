//
//  odesolver.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/9.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
#include <iostream>
using namespace std;


void solver::RK4Step()
{
    if (cRank==0)
    {
        gsl_matrix_memcpy(odetempField, Fields);
    }
    Fun(k1);
    time+=StepT/2;
    if (cRank==0)
    {
        gsl_matrix_scale(k1, StepT/2);
        gsl_matrix_add(Fields, k1);
    }
    Fun(k2);
    if (cRank==0)
    {
        gsl_matrix_scale(k2, StepT/2);
        gsl_matrix_memcpy(Fields, odetempField);
        gsl_matrix_add(Fields, k2);
    }
    Fun(k3);
    time+=StepT/2;
    if (cRank==0)
    {
        gsl_matrix_scale(k3, StepT);
        gsl_matrix_memcpy(Fields, odetempField);
        gsl_matrix_add(Fields, k3);
    }
    Fun(k4);
    if (cRank==0)
    {
        gsl_matrix_scale(k4, StepT/2);
        gsl_matrix_scale(k2, 2);
        
        gsl_matrix_add(k1, k2);
        gsl_matrix_add(k1, k3);
        gsl_matrix_add(k1, k4);
        gsl_matrix_scale(k1, 1.0/3);
        gsl_matrix_memcpy(Fields, odetempField);
        gsl_matrix_add(Fields, k1);
    }
}

void solver::BDF4Step()
{
    double error=1;
    double err;
    if (cRank==0)
    {
        for (int iterCPU=1; iterCPU<numOfProcess; ++iterCPU)
        {
            MPI_Send(Fields->data, 1, BD4Type[iterCPU], iterCPU, iterCPU, MPI_COMM_WORLD);
        }
        MPI_Isend(Fields->data, 1, BD4Type[0], 0, 0, MPI_COMM_WORLD, &request);
        MPI_Irecv(odetempField2->data, iterPoints, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
        
    }
    else
    {
        MPI_Recv(odetempField2->data, iterPoints, MPI_DOUBLE, 0, cRank, MPI_COMM_WORLD, &status);
    }
    gsl_matrix_memcpy(odetempField3, odetempField2);
    while (1)
    {
        Fun(Fields);
        if (cRank==0)
        {
            for (int iterCPU=1; iterCPU<numOfProcess; ++iterCPU)
            {
                MPI_Send(Fields->data, 1, BD4Type[iterCPU], iterCPU, iterCPU, MPI_COMM_WORLD);
            }
            MPI_Isend(Fields->data, 1, BD4Type[0], 0, 0, MPI_COMM_WORLD, &request);
            MPI_Irecv(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
            
        }
        else
        {
            MPI_Recv(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, cRank, MPI_COMM_WORLD, &status);
        }
        gsl_matrix_scale(iterFieldsLocal, StepT*0.48);
        //error=0;
        for (int iter=0; iter<iterPoints; ++iter)
        {
            iterFieldsLocal->data[iter]+=odetempField2->data[iter]*1.92-HistoryFields[2]->data[iter]*1.44+HistoryFields[1]->data[iter]*0.64-HistoryFields[0]->data[iter]*0.12;
            odetempField3->data[iter]=odetempField3->data[iter]-iterFieldsLocal->data[iter];
            if (odetempField3->data[iter]<0)
            {
                odetempField3->data[iter]=-odetempField3->data[iter];
            }
            //            if (k1->data[iter]>error)
            //            {
            //                error=k1->data[iter];
            //            }
        }
        error=gsl_matrix_max(odetempField3);
        //MPI_Bcast(&error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Allreduce(&error, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        //cout << err << endl;
        if (err<tolerance)
            break;
        
        gsl_matrix_memcpy(odetempField3, iterFieldsLocal);
        if (cRank!=0)
            MPI_Send(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, cRank, MPI_COMM_WORLD);
        else
        {
            MPI_Isend(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
            MPI_Irecv(Fields->data, 1, BD4Type[0], 0, 0, MPI_COMM_WORLD, &request);
            for (int iterCPU=1; iterCPU<numOfProcess; ++iterCPU)
            {
                MPI_Recv(Fields->data, 1, BD4Type[iterCPU], iterCPU, iterCPU, MPI_COMM_WORLD, &status);
            }
        }
        //        if (timeIdx==782)
        //        {
        //            cout << error << endl;
        //        }
    }
    
    gsl_matrix *temp=HistoryFields[0];
    HistoryFields[0]=HistoryFields[1];
    HistoryFields[1]=HistoryFields[2];
    HistoryFields[2]=odetempField2;
    odetempField2=temp;
    time+=StepT;
}


