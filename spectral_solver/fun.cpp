//
//  fun.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/9.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"

void solver::HGFuns()
{
    double f[NumField];
    for (int iter=0; iter<NumPoints; ++iter)
    {
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            f[iterf]=Fields->data[iterf*NumPoints+iter];
        }
        
        //Calculation for H functions. Use f[] to express.
        //Field 1
        Hij[0]->data[iter]=1;
        //Hij[0]->data[iter+NumPoints]=0;
        //Field 2
        //Hij[1]->data[iter]=0;
        //Hij[1]->data[iter+NumPoints]=0;
        //Field 3
        //....
        
        //Ending
        //Calculation for G function.
        G->data[iter]=0;
        //G->data[iter+NumPoints]=0;
        //G->data[iter+2*NumPoints]=0;
        //G->data[iter+3*NumPoints]=0;
        //....
    }
}


void solver::Fun(gsl_matrix *result)
{
    //H functions and G term
    HGFuns();
    
    //derivative of r term
    dr(1);
    gsl_matrix_set_zero(tempFields);
    for (int iterdf=0; iterdf<NumField; ++iterdf)
    {
        gsl_matrix_memcpy(caltempFields, Hij[iterdf]);
        for (int iter=0; iter<NumField; ++iter)
        {
            gsl_matrix_mul_elements(&ctFieldView[iter].matrix, &dFieldView[iterdf].matrix);
        }
        gsl_matrix_add(tempFields, caltempFields);
    }
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        for (int iter=0; iter<Ntheta; ++iter)
        {
            gsl_vector_view temp=gsl_matrix_subcolumn(tempFields, iter, iterf*Nrp, Nrp);
            gsl_vector_mul(&temp.vector, r);
        }
    }
    dr(0);
    gsl_matrix_scale(dFields, 1.0/radius/radius);
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        for (int iter=0; iter<Ntheta; ++iter)
        {
            gsl_vector_view temp=gsl_matrix_subcolumn(dFields, iter, iterf*Nrp, Nrp);
            gsl_vector_div(&temp.vector, r);
        }
    }
    gsl_matrix_add(G, dFields);
    
    //derivative of theta term
    dtheta(1);
    gsl_matrix_set_zero(tempFields);
    for (int iterdf=0; iterdf<NumField; ++iterdf)
    {
        gsl_matrix_memcpy(caltempFields, Hij[iterdf]);
        for (int iter=0; iter<NumField; ++iter)
        {
            gsl_matrix_mul_elements(&ctFieldView[iter].matrix, &dFieldView[iterdf].matrix);
        }
        gsl_matrix_add(tempFields, caltempFields);
    }
    dtheta(0);
    gsl_matrix_scale(dFields, 1.0/radius/radius);
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        for (int iter=0; iter<Ntheta; ++iter)
        {
            gsl_vector_view temp=gsl_matrix_subcolumn(dFields, iter, iterf*Nrp, Nrp);
            gsl_vector_div(&temp.vector, r);
            gsl_vector_div(&temp.vector, r);
        }
    }
    gsl_matrix_add(G, dFields);
    gsl_matrix_memcpy(result, G);
}
