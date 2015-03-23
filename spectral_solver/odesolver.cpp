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
    gsl_matrix_memcpy(odetempField, Fields);
    Fun(k1);
    time+=StepT/2;
    gsl_matrix_scale(k1, StepT/2);
    gsl_matrix_add(Fields, k1);
    Fun(k2);
    gsl_matrix_scale(k2, StepT/2);
    gsl_matrix_memcpy(Fields, odetempField);
    gsl_matrix_add(Fields, k2);
    Fun(k3);
    time+=StepT/2;
    gsl_matrix_scale(k3, StepT);
    gsl_matrix_memcpy(Fields, odetempField);
    gsl_matrix_add(Fields, k3);
    Fun(k4);
    gsl_matrix_scale(k4, StepT/2);
    gsl_matrix_scale(k2, 2);
    
    gsl_matrix_add(k1, k2);
    gsl_matrix_add(k1, k3);
    gsl_matrix_add(k1, k4);
    gsl_matrix_scale(k1, 1.0/3);
    gsl_matrix_memcpy(Fields, odetempField);
    gsl_matrix_add(Fields, k1);
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
