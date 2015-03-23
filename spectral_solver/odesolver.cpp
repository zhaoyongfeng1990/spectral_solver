//
//  odesolver.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/9.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"

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
        for (int iter=0; iter<totalPoints; ++iter)
        {
            k1->data[iter]*=0.5*StepT;
            Fields->data[iter]=odetempField->data[iter]+k1->data[iter];
        }
    }
    Fun(k2);
    if (cRank==0)
    {
        for (int iter=0; iter<totalPoints; ++iter)
        {
            k2->data[iter]*=0.5*StepT;
            Fields->data[iter]=odetempField->data[iter]+k2->data[iter];
        }
    }
    Fun(k3);
    time+=StepT/2;
    if (cRank==0)
    {
        for (int iter=0; iter<totalPoints; ++iter)
        {
            k3->data[iter]*=StepT;
            Fields->data[iter]=odetempField->data[iter]+k3->data[iter];
        }
    }
    Fun(k4);
    if (cRank==0)
    {
        for (int iter=0; iter<totalPoints; ++iter)
        {
            k4->data[iter]*=0.5*StepT;
            Fields->data[iter]=odetempField->data[iter]+(k1->data[iter]+k2->data[iter]*2+k3->data[iter]+k4->data[iter])/3.0;
        }
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
                MPI_Send(Fields->data, 1, BD4Type[iterCPU], iterCPU, 100+iterCPU, MPI_COMM_WORLD);
            }
            MPI_Isend(Fields->data, 1, BD4Type[0], 0, 100, MPI_COMM_WORLD, &request);
            MPI_Irecv(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD, &request);
            
        }
        else
        {
            MPI_Recv(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 100+cRank, MPI_COMM_WORLD, &status);
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
        }
        error=gsl_matrix_max(odetempField3);
        MPI_Allreduce(&error, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (cRank!=0)
            MPI_Send(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 200+cRank, MPI_COMM_WORLD);
        else
        {
            MPI_Isend(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 200, MPI_COMM_WORLD, &request);
            MPI_Irecv(Fields->data, 1, BD4Type[0], 0, 200, MPI_COMM_WORLD, &request);
            for (int iterCPU=1; iterCPU<numOfProcess; ++iterCPU)
            {
                MPI_Recv(Fields->data, 1, BD4Type[iterCPU], iterCPU, 200+iterCPU, MPI_COMM_WORLD, &status);
            }
        }
        if (err<tolerance)
            break;
        
        gsl_matrix_memcpy(odetempField3, iterFieldsLocal);
    }
    
    gsl_matrix *temp=HistoryFields[0];
    HistoryFields[0]=HistoryFields[1];
    HistoryFields[1]=HistoryFields[2];
    HistoryFields[2]=odetempField2;
    odetempField2=temp;
    time+=StepT;
}