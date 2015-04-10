//
//  derivative.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"

void solver::dr(bool ifFirst)
{
    if (ifFirst)
    {
        for (int iterr=0; iterr<Nrp; ++iterr)
        {
            for (int itert=0; itert<jobT; ++itert)
            {
                dctr.ele(Nr-iterr-1, itert)=dctr.ele(iterr, itert);
            }
        }
        
        fftwl_execute(dctr2r);
    }
    else
    {
        for (int iterr=0; iterr<Nrp; ++iterr)
        {
            for (int itert=0; itert<jobT; ++itert)
            {
                tempdctr.ele(Nr-iterr-1, itert)=tempdctr.ele(iterr, itert);
            }
        }
        
        fftwl_execute(tempdctr2r);
    }
    //The first and last row should divide 2, but since the first row will be dropped, and the last row is simply 0, so we omit it.
    
    // aliasing
    //printdebugM(&dctr, "dctr.txt");
    //for (int itert=0; itert<jobT; ++itert)
    //{
        for (int iterr=aliasingr-1; iterr<Nr; ++iterr)
        {
            //dctr.ele(iterr, itert)=0;
            dctr.setRow(iterr, 0);
        }
    //}
    
    //gsl_vector_view lastRow; //=gsl_matrix_row(dctr, aliasingr-1);
    //gsl_vector_set_zero(&lastRow.vector);
    
    for (int iterr=aliasingr-3; iterr>0; iterr-=2)
    {
        //gsl_vector_view nextLastRow=gsl_matrix_row(dctr, iterr+1);
        //lastRow=gsl_matrix_row(dctr,iterr+2);
        //gsl_vector_view cRow=gsl_matrix_row(dctr, iterr);
        
        for (int iter=0; iter<jobT; ++iter)
        {
            dctr.ele(iterr,iter)=dctr.ele(iterr+1,iter)*2.0*(iterr+1)/logicNr+dctr.ele(iterr+2,iter);
            //dctr.ele(iterr+1,iter)=0;
        }
        dctr.setRow(iterr+1, 0);
    }
    //lastRow=gsl_matrix_row(dctr, 0);
    //gsl_vector_set_zero(&lastRow.vector);
    dctr.setRow(0, 0);
    
    for (int iter=1; iter<aliasingr-2; iter+=2)
    {
        //gsl_vector_view temp=gsl_matrix_row(dctr, iter);
        //gsl_vector_scale(&temp.vector, 0.5);
        dctr.scaleRow(iter, 0.5);
    }
    
    fftwl_execute(dctr2r);
}

void solver::drWOA()
{
    for (int iterr=0; iterr<Nrp; ++iterr)
    {
        for (int itert=0; itert<jobT; ++itert)
        {
            dctr.ele(Nr-iterr-1, itert)=dctr.ele(iterr, itert);
        }
    }
    
    fftwl_execute(dctr2r);
    //The first and last row should divide 2, but since the first row will be dropped, and the last row is simply 0, so we omit it.
    
    //gsl_vector_view lastRow=gsl_matrix_row(dctr, Nr-1);
    //gsl_vector_set_zero(&lastRow.vector);
    dctr.setRow(Nr-1, 0);
    for (int iterr=Nr-3; iterr>0; iterr-=2)
    {
        //gsl_vector_view nextLastRow=gsl_matrix_row(dctr, iterr+1);
        //lastRow=gsl_matrix_row(dctr,iterr+2);
        //gsl_vector_view cRow=gsl_matrix_row(dctr, iterr);
        
        for (int iter=0; iter<jobT; ++iter)
        {
            dctr.ele(iterr,iter)=dctr.ele(iterr+1,iter)*2.0*(iterr+1)/logicNr+dctr.ele(iterr+2,iter);
            //dctr.ele(iterr+1,iter)=0;
        }
        dctr.setRow(iterr+1, 0);
    }
    //lastRow=gsl_matrix_row(dctr, 0);
    //gsl_vector_set_zero(&lastRow.vector);
    dctr.setRow(0, 0);
    
    for (int iter=1; iter<Nr-2; iter+=2)
    {
        //gsl_vector_view temp=gsl_matrix_row(dctr, iter);
        //gsl_vector_scale(&temp.vector, 0.5);
        dctr.scaleRow(iter, 0.5);
    }
    
    fftwl_execute(dctr2r);
}

void solver::dtheta(bool ifFirst)
{
    
    if (ifFirst)
    {
        fftwl_execute(fftr2c);
    }
    else
    {
        fftwl_execute(tempfftr2c);
    }
    
    //for (int iterr=0; iterr<jobR; ++iterr)
    //{
        // aliasing
        for (int itertheta=0; itertheta<aliasingt; ++itertheta)
        {
            //gsl_complex temp=gsl_matrix_complex_get(fftc, iterr, itertheta);
            //a'=a*ik
            //long double swap=temp.dat[0];
            //temp.dat[0]=-temp.dat[1]*itertheta/Ntheta;
            //temp.dat[1]=swap*itertheta/Ntheta;
            //gsl_matrix_complex_set(fftc, iterr, itertheta, temp);
            fftc.scaleCol(itertheta, 0, (long double)itertheta/Ntheta);
        }
        for (int itertheta=aliasingt; itertheta<Ntheta/2+1; ++itertheta)
        {
            //fftc.data[(iterr*(Ntheta/2+1)+itertheta)*2]=0;
            //fftc.data[(iterr*(Ntheta/2+1)+itertheta)*2+1]=0;
            fftc.setCol(itertheta, 0, 0);
        }
    //}
    
    if (ifFirst)
    {
        fftwl_execute(ifftc2r);
    }
    else
    {
        fftwl_execute(tempifftc2r);
    }
}
