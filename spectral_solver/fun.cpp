//
//  fun.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/9.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
using std::pow;

void solver::HGFuns()
{
    double f[NumField];
//#ifdef MULTIPROCESS
//#pragma omp parallel for
//#endif
    for (int iter=0; iter<NumPoints; ++iter)
    {
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            f[iterf]=Fields->data[iterf*NumPoints+iter];
        }
        
        //Calculation for H functions. Use f[] to express.
#ifdef FU_MODEL
        //Term in front of dField 1
        double powfm1=pow(f[1]/Kh,M-1);
        double powfm=powfm1*f[1]/Kh;
        double powfmp1=powfm+1;
        Hij[0]->data[iter]=(Drho+Drho0*powfm)/powfmp1;
        Hij[0]->data[iter+NumPoints]=0;
        Hij[0]->data[iter+2*NumPoints]=0;
        //Term in front of dField 2
        Hij[1]->data[iter]=f[0]*(Drho0-Drho)*powfm1*M/Kh/powfmp1/powfmp1;
        Hij[1]->data[iter+NumPoints]=Dh;
        Hij[1]->data[iter+2*NumPoints]=0;
        //Term in front of dField 3
        Hij[2]->data[iter]=0;
        Hij[2]->data[iter+NumPoints]=0;
        Hij[2]->data[iter+2*NumPoints]=Dn;
        
        double f2=f[2]*f[2];
        //Ending
        //Calculation for G function.
        G->data[iter]=Gamma*f2*f[0]/(f2+Kn2);
        G->data[iter+NumPoints]=Alpha*f[0]-Beta*f[1];
        G->data[iter+2*NumPoints]=-G->data[iter];
#endif
#ifdef LINEAR_TEST_MODEL
        Hij[0]->data[iter]=1;
        Hij[0]->data[iter+NumPoints]=0;
        Hij[0]->data[iter+2*NumPoints]=0;
        //Term in front of dField 2
        Hij[1]->data[iter]=0;
        Hij[1]->data[iter+NumPoints]=1;
        Hij[1]->data[iter+2*NumPoints]=0;
        //Term in front of dField 3
        Hij[2]->data[iter]=0;
        Hij[2]->data[iter+NumPoints]=0;
        Hij[2]->data[iter+2*NumPoints]=1;
        
        //Ending
        //Calculation for G function.
        G->data[iter]=0;
        G->data[iter+NumPoints]=0;
        G->data[iter+2*NumPoints]=0;
#endif
        //....
    }
}


void solver::Fun(gsl_matrix *result)
{
    
//#ifdef MULTIPROCESS
//    omp_set_num_threads(8);
//#endif
    
    //H functions and G term
    HGFuns();
    
    //derivative of r term
    dr(1);
    
#ifdef PUNISHTERM
//#ifdef MULTIPROCESS
//#pragma omp parallel for
//#endif
    for (int iter=0; iter<NumField; ++iter)
    {
        gsl_vector_view tempboundary=gsl_matrix_row(dFields, iter*Nrp);
        gsl_vector_view destiny=gsl_matrix_row(boundary, iter);
        gsl_vector_memcpy(&destiny.vector, &tempboundary.vector);
    }
    gsl_matrix_scale(boundary, punish);
#endif
    
    gsl_matrix_set_zero(tempFields);
//#ifdef MULTIPROCESS
//#pragma omp parallel for
//#endif
    for (int iterdf=0; iterdf<NumField; ++iterdf)
    {
        gsl_matrix_memcpy(caltempFields, Hij[iterdf]);
        for (int iter=0; iter<NumField; ++iter)
        {
            gsl_matrix_mul_elements(&ctFieldView[iter].matrix, &dFieldView[iterdf].matrix);
        }
        gsl_matrix_add(tempFields, caltempFields);
    }
//#ifdef MULTIPROCESS
//#pragma omp parallel for
//#endif
    for (int iter=0; iter<Ntheta; ++iter)
    {
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            gsl_vector_view temp=gsl_matrix_subcolumn(tempFields, iter, iterf*Nrp, Nrp);
            gsl_vector_mul(&temp.vector, r);
        }
    }
    
    dr(0);
    gsl_matrix_scale(dFields, rp2); //rp2=1.0/r/r
//#ifdef MULTIPROCESS
//#pragma omp parallel for
//#endif
    for (int iter=0; iter<Ntheta; ++iter)
    {
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            gsl_vector_view temp=gsl_matrix_subcolumn(dFields, iter, iterf*Nrp, Nrp);
            gsl_vector_div(&temp.vector, r);
        }
    }
    gsl_matrix_add(G, dFields);
    
    //derivative of theta term
    dtheta(1);
    gsl_matrix_set_zero(tempFields);
//#ifdef MULTIPROCESS
//#pragma omp parallel for
//#endif
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
    gsl_matrix_scale(dFields, rp2); //rp2=1.0/r/r
//#ifdef MULTIPROCESS
//#pragma omp parallel for
//#endif
    for (int iter=0; iter<Ntheta; ++iter)
    {
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            gsl_vector_view temp=gsl_matrix_subcolumn(dFields, iter, iterf*Nrp, Nrp);
            gsl_vector_div(&temp.vector, r2);
        }
    }
    gsl_matrix_add(G, dFields);
    
#ifdef PUNISHTERM
//#ifdef MULTIPROCESS
//#pragma omp parallel for
//#endif
    for (int iter=0; iter<NumField; ++iter)
    {
        gsl_vector_view tempboundary=gsl_matrix_row(G, iter*Nrp);
        gsl_vector_view destiny=gsl_matrix_row(boundary, iter);
        gsl_vector_sub(&tempboundary.vector, &destiny.vector);
    }
#endif
    
    gsl_matrix_memcpy(result, G);
}
