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

void solver::RK6Step()
{
    if (cRank==0)
    {
        gsl_matrix_memcpy(odetempField, Fields);
    }
    Fun(k1);
    double ctime=time;
    time=ctime+StepT*9.0/50.0;
    
    if (cRank==0)
    {
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields->data[iter]=odetempField->data[iter]+StepT*k1->data[iter]*9.0/50.0;
        }
    }
    Fun(k2);
    time=ctime+StepT/6.0;
    if (cRank==0)
    {
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields->data[iter]=odetempField->data[iter]+StepT*(k1->data[iter]*29.0/324.0+k2->data[iter]*25.0/324.0);
        }
    }
    Fun(k3);
    time=ctime+StepT/4.0;
    
    if (cRank==0)
    {
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields->data[iter]=odetempField->data[iter]+StepT*(k1->data[iter]/16.0+k3->data[iter]*3.0/16.0);
        }
    }
    Fun(k4);
    time=ctime+StepT*0.53;
    
    if (cRank==0)
    {
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields->data[iter]=odetempField->data[iter]+StepT*(k1->data[iter]*79129.0/250000.0-k3->data[iter]*261237.0/250000.0+k4->data[iter]*19663.0/15625.0);
        }
    }
    Fun(k5);
    time=ctime+StepT*0.6;
    
    if (cRank==0)
    {
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields->data[iter]=odetempField->data[iter]+StepT*(k1->data[iter]*1336883.0/4909125.0-k3->data[iter]*25476.0/30875.0+k4->data[iter]*194159.0/185250.0+k5->data[iter]*8225.0/78546.0);
        }
    }
    Fun(k6);
    time=ctime+StepT*0.8;
    
    if (cRank==0)
    {
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields->data[iter]=odetempField->data[iter]+StepT*(-k1->data[iter]*2459386.0/14727375.0+k3->data[iter]*19504.0/30875.0+k4->data[iter]*2377474.0/13615875.0-k5->data[iter]*6157250.0/5773131.0+k6->data[iter]*902.0/735.0);
        }
    }
    Fun(k7);
    time=ctime+StepT;
    if (cRank==0)
    {
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields->data[iter]=odetempField->data[iter]+StepT*(k1->data[iter]*2699.0/7410.0-k3->data[iter]*252.0/1235.0-k4->data[iter]*1393253.0/3993990.0+k5->data[iter]*236875.0/72618.0-k6->data[iter]*135.0/49.0+k7->data[iter]*15.0/22.0);
        }
    }
    Fun(k8);
    
    if (cRank==0)
    {
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields->data[iter]=odetempField->data[iter]+StepT*(k1->data[iter]*11.0/144.0+k4->data[iter]*256.0/693.0+k6->data[iter]*125.0/504.0+k7->data[iter]*125.0/528.0+k8->data[iter]*5.0/72.0);
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
            MPI_Ssend(Fields->data, 1, BD4Type[iterCPU], iterCPU, iterCPU, MPI_COMM_WORLD);
        }
        //MPI_Sendrecv(Fields->data, 1, BD4Type[0], 0, 0, odetempField2->data, iterPoints, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            for (int iter=0; iter<bossP; ++iter)
            {
                odetempField2->data[iter+iterf*bossP]=Fields->data[iter+iterf*NumPoints];
            }
        }
        
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
                MPI_Ssend(Fields->data, 1, BD4Type[iterCPU], iterCPU, 100+iterCPU, MPI_COMM_WORLD);
            }
            //MPI_Sendrecv(Fields->data, 1, BD4Type[0], 0, 100, iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD, &status);
            for (int iterf=0; iterf<NumField; ++iterf)
            {
                for (int iter=0; iter<bossP; ++iter)
                {
                    iterFieldsLocal->data[iter+iterf*bossP]=Fields->data[iter+iterf*NumPoints];
                }
            }
            
        }
        else
        {
            MPI_Recv(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 100+cRank, MPI_COMM_WORLD, &status);
        }
        //gsl_matrix_scale(iterFieldsLocal, StepT*0.48);
        //error=0;
        for (int iter=0; iter<iterPoints; ++iter)
        {
            iterFieldsLocal->data[iter]*=0.48*StepT;
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
            MPI_Ssend(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 200+cRank, MPI_COMM_WORLD);
        else
        {
            //MPI_Sendrecv(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 200,Fields->data, 1, BD4Type[0], 0, 200, MPI_COMM_WORLD, &status);
            for (int iterf=0; iterf<NumField; ++iterf)
            {
                for (int iter=0; iter<bossP; ++iter)
                {
                    Fields->data[iter+iterf*NumPoints]=iterFieldsLocal->data[iter+iterf*bossP];
                }
            }
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

void solver::BDF6Step()
{
    double error=1;
    double err;
    if (cRank==0)
    {
        for (int iterCPU=1; iterCPU<numOfProcess; ++iterCPU)
        {
            MPI_Ssend(Fields->data, 1, BD4Type[iterCPU], iterCPU, iterCPU, MPI_COMM_WORLD);
        }
        //MPI_Sendrecv(Fields->data, 1, BD4Type[0], 0, 0, odetempField2->data, iterPoints, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            for (int iter=0; iter<bossP; ++iter)
            {
                odetempField2->data[iter+iterf*bossP]=Fields->data[iter+iterf*NumPoints];
            }
        }
        
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
                MPI_Ssend(Fields->data, 1, BD4Type[iterCPU], iterCPU, 100+iterCPU, MPI_COMM_WORLD);
            }
            //MPI_Sendrecv(Fields->data, 1, BD4Type[0], 0, 100, iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD, &status);
            for (int iterf=0; iterf<NumField; ++iterf)
            {
                for (int iter=0; iter<bossP; ++iter)
                {
                    iterFieldsLocal->data[iter+iterf*bossP]=Fields->data[iter+iterf*NumPoints];
                }
            }
            
        }
        else
        {
            MPI_Recv(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 100+cRank, MPI_COMM_WORLD, &status);
        }
        for (int iter=0; iter<iterPoints; ++iter)
        {
            iterFieldsLocal->data[iter]*=StepT*60.0/147.0;
            iterFieldsLocal->data[iter]+=odetempField2->data[iter]*360.0/147.0-HistoryFields[4]->data[iter]*450.0/147.0+HistoryFields[3]->data[iter]*400.0/147.0-HistoryFields[2]->data[iter]*225.0/147.0+HistoryFields[1]->data[iter]*72.0/147.0-HistoryFields[0]->data[iter]*10.0/147.0;
            odetempField3->data[iter]=odetempField3->data[iter]-iterFieldsLocal->data[iter];
            if (odetempField3->data[iter]<0)
            {
                odetempField3->data[iter]=-odetempField3->data[iter];
            }
        }
        error=gsl_matrix_max(odetempField3);
        MPI_Allreduce(&error, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (cRank!=0)
            MPI_Ssend(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 200+cRank, MPI_COMM_WORLD);
        else
        {
            //MPI_Sendrecv(iterFieldsLocal->data, iterPoints, MPI_DOUBLE, 0, 200,Fields->data, 1, BD4Type[0], 0, 200, MPI_COMM_WORLD, &status);
            for (int iterf=0; iterf<NumField; ++iterf)
            {
                for (int iter=0; iter<bossP; ++iter)
                {
                    Fields->data[iter+iterf*NumPoints]=iterFieldsLocal->data[iter+iterf*bossP];
                }
            }
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
    HistoryFields[2]=HistoryFields[3];
    HistoryFields[3]=HistoryFields[4];
    HistoryFields[4]=odetempField2;
    odetempField2=temp;
    
    time+=StepT;
}

