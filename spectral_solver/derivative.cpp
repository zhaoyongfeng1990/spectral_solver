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
            gsl_vector_view datablock=gsl_matrix_row(dctr, iterr);
            gsl_vector_view destiny=gsl_matrix_row(dctr, Nr-iterr-1);
            gsl_vector_memcpy(&destiny.vector, &datablock.vector);
        }
        
        fftw_execute(dctr2r);
    }
    else
    {
        for (int iterr=0; iterr<Nrp; ++iterr)
        {
            gsl_vector_view datablock=gsl_matrix_row(tempdctr, iterr);
            gsl_vector_view destiny=gsl_matrix_row(tempdctr, Nr-iterr-1);
            gsl_vector_memcpy(&destiny.vector, &datablock.vector);
        }
        
        fftw_execute(tempdctr2r);
    }
    //The first and last row should divide 2, but since the first row will be dropped, and the last row is simply 0, so we omit it.
    
    // aliasing
    for (int itert=0; itert<jobT; ++itert)
    {
        for (int iterr=aliasingr; iterr<Nr; ++iterr)
        {
            gsl_matrix_set(dctr, iterr, itert, 0);
        }
    }
    
    gsl_vector_view lastRow=gsl_matrix_row(dctr, aliasingr-1);
    gsl_vector_set_zero(&lastRow.vector);
    
    for (int iterr=aliasingr-3; iterr>0; iterr-=2)
    {
        gsl_vector_view nextLastRow=gsl_matrix_row(dctr, iterr+1);
        lastRow=gsl_matrix_row(dctr,iterr+2);
        gsl_vector_view cRow=gsl_matrix_row(dctr, iterr);
        
        for (int iter=0; iter<jobT; ++iter)
        {
            cRow.vector.data[iter]=nextLastRow.vector.data[iter]*2.0*(iterr+1)/logicNr+lastRow.vector.data[iter];
            nextLastRow.vector.data[iter]=0;
        }
    }
    lastRow=gsl_matrix_row(dctr, 0);
    gsl_vector_set_zero(&lastRow.vector);
    
    for (int iter=1; iter<aliasingr-2; iter+=2)
    {
        gsl_vector_view temp=gsl_matrix_row(dctr, iter);
        gsl_vector_scale(&temp.vector, 0.5);
    }
    
    fftw_execute(dctr2r);
}

void solver::dtheta(bool ifFirst)
{
    
    if (ifFirst)
    {
        fftw_execute(fftr2c);
    }
    else
    {
        fftw_execute(tempfftr2c);
    }
    
    for (int iterr=0; iterr<jobR; ++iterr)
    {
        // aliasing
        for (int itertheta=0; itertheta<aliasingt; ++itertheta)
        {
            gsl_complex temp=gsl_matrix_complex_get(fftc, iterr, itertheta);
            //a'=a*ik
            double swap=temp.dat[0];
            temp.dat[0]=-temp.dat[1]*itertheta/Ntheta;
            temp.dat[1]=swap*itertheta/Ntheta;
            gsl_matrix_complex_set(fftc, iterr, itertheta, temp);
        }
        for (int itertheta=aliasingt; itertheta<Ntheta/2+1; ++itertheta)
        {
            fftc->data[(iterr*(Ntheta/2+1)+itertheta)*2]=0;
            fftc->data[(iterr*(Ntheta/2+1)+itertheta)*2+1]=0;
        }
    }
    
    if (ifFirst)
    {
        fftw_execute(ifftc2r);
    }
    else
    {
        fftw_execute(tempifftc2r);
    }
}