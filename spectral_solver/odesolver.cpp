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
    gsl_matrix_memcpy(odetempField, Fields);
    while (error>tolerance)
    {
        Fun(Fields);
        gsl_matrix_scale(Fields, StepT*0.48);
        for (int iter=0; iter<matrixH*Ntheta; ++iter)
        {
            Fields->data[iter]+=odetempField->data[iter]*1.92-HistoryFields[2]->data[iter]*1.44+HistoryFields[1]->data[iter]*0.64-HistoryFields[0]->data[iter]*0.12;
            tempFields->data[iter]=abs(Fields->data[iter]-odetempField->data[iter]);
        }
        error=gsl_matrix_max(tempFields);
    }
    gsl_matrix_memcpy(HistoryFields[0], HistoryFields[1]);
    gsl_matrix_memcpy(HistoryFields[1], HistoryFields[2]);
    gsl_matrix_memcpy(HistoryFields[2], odetempField);
    
    time+=StepT;
}
