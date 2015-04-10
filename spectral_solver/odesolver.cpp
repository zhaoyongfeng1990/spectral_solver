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
    odetempField=Fields;
    Fun(k1);
    time+=StepT/2;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        k1[iter]*=0.5*StepT;
        Fields[iter]=odetempField[iter]+k1[iter];
    }
    Fun(k2);
    for (int iter=0; iter<totalPoints; ++iter)
    {
        k2[iter]*=0.5*StepT;
        Fields[iter]=odetempField[iter]+k2[iter];
    }
    Fun(k3);
    time+=StepT/2;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        k3[iter]*=StepT;
        Fields[iter]=odetempField[iter]+k3[iter];
    }
    Fun(k4);
    for (int iter=0; iter<totalPoints; ++iter)
    {
        k4[iter]*=0.5*StepT;
        Fields[iter]=odetempField[iter]+(k1[iter]+k2[iter]*2+k3[iter]+k4[iter])/3.0;
    }
}

void solver::RK6Step()
{
    odetempField=Fields;
    Fun(k1);
    long double ctime=time;
    time=ctime+StepT*9.0/50.0;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields[iter]=odetempField[iter]+StepT*k1[iter]*9.0/50.0;
    }
    Fun(k2);
    time=ctime+StepT/6.0;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields[iter]=odetempField[iter]+StepT*(k1[iter]*29.0/324.0+k2[iter]*25.0/324.0);
    }
    Fun(k3);
    time=ctime+StepT/4.0;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields[iter]=odetempField[iter]+StepT*(k1[iter]/16.0+k3[iter]*3.0/16.0);
    }
    Fun(k4);
    time=ctime+StepT*0.53;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields[iter]=odetempField[iter]+StepT*(k1[iter]*79129.0/250000.0-k3[iter]*261237.0/250000.0+k4[iter]*19663.0/15625.0);
    }
    Fun(k5);
    time=ctime+StepT*0.6;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields[iter]=odetempField[iter]+StepT*(k1[iter]*1336883.0/4909125.0-k3[iter]*25476.0/30875.0+k4[iter]*194159.0/185250.0+k5[iter]*8225.0/78546.0);
    }
    Fun(k6);
    time=ctime+StepT*0.8;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields[iter]=odetempField[iter]+StepT*(-k1[iter]*2459386.0/14727375.0+k3[iter]*19504.0/30875.0+k4[iter]*2377474.0/13615875.0-k5[iter]*6157250.0/5773131.0+k6[iter]*902.0/735.0);
    }
    Fun(k7);
    time=ctime+StepT;
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields[iter]=odetempField[iter]+StepT*(k1[iter]*2699.0/7410.0-k3[iter]*252.0/1235.0-k4[iter]*1393253.0/3993990.0+k5[iter]*236875.0/72618.0-k6[iter]*135.0/49.0+k7[iter]*15.0/22.0);
    }
    Fun(k8);
    for (int iter=0; iter<totalPoints; ++iter)
    {
        Fields[iter]=odetempField[iter]+StepT*(k1[iter]*11.0/144.0+k4[iter]*256.0/693.0+k6[iter]*125.0/504.0+k7[iter]*125.0/528.0+k8[iter]*5.0/72.0);
    }
}

void solver::BDF4Step()
{
    long double error=1;
    odetempField2=Fields;
    while (error>tolerance)
    {
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            k1[iter]=Fields[iter];
        }
        //matrix<long double>_memcpy(k1, Fields);
        Fun(Fields);
        //matrix<long double>_scale(Fields, StepT*0.48);
        error=0;
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields[iter]*=0.48*StepT;
            Fields[iter]+=odetempField2[iter]*1.92-HistoryFields[2]->data[iter]*1.44+HistoryFields[1]->data[iter]*0.64-HistoryFields[0]->data[iter]*0.12;
            k1[iter]-=Fields[iter];
            if (k1[iter]<0)
            {
                k1[iter]=-k1[iter];
            }
            if (k1[iter]>error)
            {
                error=k1[iter];
            }
        }
        //error=matrix<long double>_max(k1);
    }
    
    long double* temp=HistoryFields[0]->data;
    HistoryFields[0]->data=HistoryFields[1]->data;
    HistoryFields[1]->data=HistoryFields[2]->data;
    HistoryFields[2]->data=odetempField2.data;
    odetempField2.data=temp;
    
    time+=StepT;
}

void solver::BDF3Step()
{
    long double error=1;
    odetempField2=Fields;
    while (error>tolerance)
    {
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            k1[iter]=Fields[iter];
        }
        //matrix<long double>_memcpy(k1, Fields);
        Fun(Fields);
        //matrix<long double>_scale(Fields, StepT*0.48);
        error=0;
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields[iter]*=StepT*6.0/11.0;
            Fields[iter]+=odetempField2[iter]*18.0/11.0-HistoryFields[1]->data[iter]*9.0/11.0+HistoryFields[0]->data[iter]*2.0/11.0;
            k1[iter]=k1[iter]-Fields[iter];
            if (k1[iter]<0)
            {
                k1[iter]=-k1[iter];
            }
            if (k1[iter]>error)
            {
                error=k1[iter];
            }
        }
        //error=matrix<long double>_max(k1);
    }
    
    long double* temp=HistoryFields[0]->data;
    HistoryFields[0]->data=HistoryFields[1]->data;
    HistoryFields[1]->data=odetempField2.data;
    odetempField2.data=temp;
    
    time+=StepT;
}

void solver::BDF5Step()
{
    long double error=1;
    odetempField2=Fields;
    while (error>tolerance)
    {
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            k1[iter]=Fields[iter];
        }
        //matrix<long double>_memcpy(k1, Fields);
        Fun(Fields);
        //matrix<long double>_scale(Fields, StepT*0.48);
        error=0;
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields[iter]*=StepT*60.0/137.0;
            Fields[iter]+=odetempField2[iter]*300.0/137.0-HistoryFields[3]->data[iter]*300.0/137.0+HistoryFields[2]->data[iter]*200.0/137.0-HistoryFields[1]->data[iter]*75.0/137.0+HistoryFields[0]->data[iter]*12.0/137.0;
            k1[iter]=k1[iter]-Fields[iter];
            if (k1[iter]<0)
            {
                k1[iter]=-k1[iter];
            }
            if (k1[iter]>error)
            {
                error=k1[iter];
            }
        }
        //error=matrix<long double>_max(k1);
    }
    
    long double* temp=HistoryFields[0]->data;
    HistoryFields[0]->data=HistoryFields[1]->data;
    HistoryFields[1]->data=HistoryFields[2]->data;
    HistoryFields[2]->data=HistoryFields[3]->data;
    HistoryFields[3]->data=odetempField2.data;
    odetempField2.data=temp;
    
    time+=StepT;
}

void solver::BDF6Step()
{
    long double error=1;
    odetempField2=Fields;
    while (error>tolerance)
    {
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            k1[iter]=Fields[iter];
        }
        //matrix<long double>_memcpy(k1, Fields);
        Fun(Fields);
        //matrix<long double>_scale(Fields, StepT*0.48);
        error=0;
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            Fields[iter]*=StepT*60.0/147.0;
            Fields[iter]+=odetempField2[iter]*360.0/147.0-HistoryFields[4]->data[iter]*450.0/147.0+HistoryFields[3]->data[iter]*400.0/147.0-HistoryFields[2]->data[iter]*225.0/147.0+HistoryFields[1]->data[iter]*72.0/147.0-HistoryFields[0]->data[iter]*10.0/147.0;
            k1[iter]=k1[iter]-Fields[iter];
            if (k1[iter]<0)
            {
                k1[iter]=-k1[iter];
            }
            if (k1[iter]>error)
            {
                error=k1[iter];
            }
        }
        //error=matrix<long double>_max(k1);
    }
    
    long double* temp=HistoryFields[0]->data;
    HistoryFields[0]->data=HistoryFields[1]->data;
    HistoryFields[1]->data=HistoryFields[2]->data;
    HistoryFields[2]->data=HistoryFields[3]->data;
    HistoryFields[3]->data=HistoryFields[4]->data;
    HistoryFields[4]->data=odetempField2.data;
    odetempField2.data=temp;
    
    time+=StepT;
}

