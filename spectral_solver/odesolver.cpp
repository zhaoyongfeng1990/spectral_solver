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
    gsl_matrix_memcpy(odetempField, Fields);
    Fun(k1);
    time+=StepT/2;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        k1->data[iter]*=0.5*StepT;
        Fields->data[iter]=odetempField->data[iter]+k1->data[iter];
    }
    Fun(k2);
    for (int iter=0; iter<totalPoints; ++iter)
    {
        k2->data[iter]*=0.5*StepT;
        Fields->data[iter]=odetempField->data[iter]+k2->data[iter];
    }
    Fun(k3);
    time+=StepT/2;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        k3->data[iter]*=StepT;
        Fields->data[iter]=odetempField->data[iter]+k3->data[iter];
    }
    Fun(k4);
    for (int iter=0; iter<totalPoints; ++iter)
    {
        k4->data[iter]*=0.5*StepT;
        Fields->data[iter]=odetempField->data[iter]+(k1->data[iter]+k2->data[iter]*2+k3->data[iter]+k4->data[iter])/3.0;
    }
}

void solver::RK6Step()
{
    gsl_matrix_memcpy(odetempField, Fields);
    Fun(k1);
    double ctime=time;
    time=ctime+StepT*9.0/50.0;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields->data[iter]=odetempField->data[iter]+StepT*k1->data[iter]*9.0/50.0;
    }
    Fun(k2);
    time=ctime+StepT/6.0;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields->data[iter]=odetempField->data[iter]+StepT*(k1->data[iter]*29.0/324.0+k2->data[iter]*25.0/324.0);
    }
    Fun(k3);
    time=ctime+StepT/4.0;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields->data[iter]=odetempField->data[iter]+StepT*(k1->data[iter]/16.0+k3->data[iter]*3.0/16.0);
    }
    Fun(k4);
    time=ctime+StepT*0.53;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields->data[iter]=odetempField->data[iter]+StepT*(k1->data[iter]*79129.0/250000.0-k3->data[iter]*261237.0/250000.0+k4->data[iter]*19663.0/15625.0);
    }
    Fun(k5);
    time=ctime+StepT*0.6;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields->data[iter]=odetempField->data[iter]+StepT*(k1->data[iter]*1336883.0/4909125.0-k3->data[iter]*25476.0/30875.0+k4->data[iter]*194159.0/185250.0+k5->data[iter]*8225.0/78546.0);
    }
    Fun(k6);
    time=ctime+StepT*0.8;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields->data[iter]=odetempField->data[iter]+StepT*(-k1->data[iter]*2459386.0/14727375.0+k3->data[iter]*19504.0/30875.0+k4->data[iter]*2377474.0/13615875.0-k5->data[iter]*6157250.0/5773131.0+k6->data[iter]*902.0/735.0);
    }
    Fun(k7);
    time=ctime+StepT;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields->data[iter]=odetempField->data[iter]+StepT*(k1->data[iter]*2699.0/7410.0-k3->data[iter]*252.0/1235.0-k4->data[iter]*1393253.0/3993990.0+k5->data[iter]*236875.0/72618.0-k6->data[iter]*135.0/49.0+k7->data[iter]*15.0/22.0);
    }
    Fun(k8);
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields->data[iter]=odetempField->data[iter]+StepT*(k1->data[iter]*11.0/144.0+k4->data[iter]*256.0/693.0+k6->data[iter]*125.0/504.0+k7->data[iter]*125.0/528.0+k8->data[iter]*5.0/72.0);
    }
}

void solver::BDF4Step()
{
    double error=1;
    gsl_matrix_memcpy(odetempField2, Fields);
    while (error>tolerance)
    {
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            k1->data[iter]=Fields->data[iter];
        }
        //gsl_matrix_memcpy(k1, Fields);
        Fun(Fields);
        //gsl_matrix_scale(Fields, StepT*0.48);
        error=0;
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields->data[iter]*=0.48*StepT;
            Fields->data[iter]+=odetempField2->data[iter]*1.92-HistoryFields[2]->data[iter]*1.44+HistoryFields[1]->data[iter]*0.64-HistoryFields[0]->data[iter]*0.12;
            k1->data[iter]=k1->data[iter]-Fields->data[iter];
            if (k1->data[iter]<0)
            {
                k1->data[iter]=-k1->data[iter];
            }
            if (k1->data[iter]>error)
            {
                error=k1->data[iter];
            }
        }
        //error=gsl_matrix_max(k1);
    }
    
    gsl_matrix *temp=HistoryFields[0];
    HistoryFields[0]=HistoryFields[1];
    HistoryFields[1]=HistoryFields[2];
    HistoryFields[2]=odetempField2;
    odetempField2=temp;
    
    time+=StepT;
}

void solver::BDF3Step()
{
    double error=1;
    gsl_matrix_memcpy(odetempField2, Fields);
    while (error>tolerance)
    {
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            k1->data[iter]=Fields->data[iter];
        }
        //gsl_matrix_memcpy(k1, Fields);
        Fun(Fields);
        //gsl_matrix_scale(Fields, StepT*0.48);
        error=0;
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields->data[iter]*=StepT*6.0/11.0;
            Fields->data[iter]+=odetempField2->data[iter]*18.0/11.0-HistoryFields[1]->data[iter]*9.0/11.0+HistoryFields[0]->data[iter]*2.0/11.0;
            k1->data[iter]=k1->data[iter]-Fields->data[iter];
            if (k1->data[iter]<0)
            {
                k1->data[iter]=-k1->data[iter];
            }
            if (k1->data[iter]>error)
            {
                error=k1->data[iter];
            }
        }
        //error=gsl_matrix_max(k1);
    }
    
    gsl_matrix *temp=HistoryFields[0];
    HistoryFields[0]=HistoryFields[1];
    HistoryFields[1]=odetempField2;
    odetempField2=temp;
    
    time+=StepT;
}

void solver::BDF5Step()
{
    double error=1;
    gsl_matrix_memcpy(odetempField2, Fields);
    while (error>tolerance)
    {
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            k1->data[iter]=Fields->data[iter];
        }
        //gsl_matrix_memcpy(k1, Fields);
        Fun(Fields);
        //gsl_matrix_scale(Fields, StepT*0.48);
        error=0;
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields->data[iter]*=StepT*60.0/137.0;
            Fields->data[iter]+=odetempField2->data[iter]*300.0/137.0-HistoryFields[3]->data[iter]*300.0/137.0+HistoryFields[2]->data[iter]*200.0/137.0-HistoryFields[1]->data[iter]*75.0/137.0+HistoryFields[0]->data[iter]*12.0/137.0;
            k1->data[iter]=k1->data[iter]-Fields->data[iter];
            if (k1->data[iter]<0)
            {
                k1->data[iter]=-k1->data[iter];
            }
            if (k1->data[iter]>error)
            {
                error=k1->data[iter];
            }
        }
        //error=gsl_matrix_max(k1);
    }
    
    gsl_matrix *temp=HistoryFields[0];
    HistoryFields[0]=HistoryFields[1];
    HistoryFields[1]=HistoryFields[2];
    HistoryFields[2]=HistoryFields[3];
    HistoryFields[3]=odetempField2;
    odetempField2=temp;
    
    time+=StepT;
}

void solver::BDF6Step()
{
    double error=1;
    gsl_matrix_memcpy(odetempField2, Fields);
    while (error>tolerance)
    {
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            k1->data[iter]=Fields->data[iter];
        }
        //gsl_matrix_memcpy(k1, Fields);
        Fun(Fields);
        //gsl_matrix_scale(Fields, StepT*0.48);
        error=0;
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields->data[iter]*=StepT*60.0/147.0;
            Fields->data[iter]+=odetempField2->data[iter]*360.0/147.0-HistoryFields[4]->data[iter]*450.0/147.0+HistoryFields[3]->data[iter]*400.0/147.0-HistoryFields[2]->data[iter]*225.0/147.0+HistoryFields[1]->data[iter]*72.0/147.0-HistoryFields[0]->data[iter]*10.0/147.0;
            k1->data[iter]=k1->data[iter]-Fields->data[iter];
            if (k1->data[iter]<0)
            {
                k1->data[iter]=-k1->data[iter];
            }
            if (k1->data[iter]>error)
            {
                error=k1->data[iter];
            }
        }
        //error=gsl_matrix_max(k1);
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

