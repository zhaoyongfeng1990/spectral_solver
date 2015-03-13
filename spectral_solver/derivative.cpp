//
//  derivative.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include <sstream>
using namespace std;

#include "solver.h"
#include <cmath>

void solver::dr(bool ifFirst)
{
    gsl_matrix* cFields;
    if (ifFirst)
    {
        cFields=Fields;
    }
    else
    {
        cFields=tempFields;
    }
    
    
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        gsl_matrix_view datablock=gsl_matrix_submatrix(cFields, iterf*Nrp, 0, Nrp, Ntheta);
        gsl_matrix_view destiny=gsl_matrix_submatrix(dctr, 0, iterf*Ntheta, Nrp, Ntheta);
        gsl_matrix_memcpy(&destiny.matrix, &datablock.matrix);
    }
    
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iterr=0; iterr<Nrp; ++iterr)
    {
        gsl_vector_view datablock=gsl_matrix_row(dctr, iterr);
        gsl_vector_view destiny=gsl_matrix_row(dctr, Nr-iterr-1);
        gsl_vector_memcpy(&destiny.vector, &datablock.vector);
    }
    
    fftw_execute(dctr2r);
    //gsl_matrix_scale(dctr, 1.0/logicNr);
    //The first and last row should divide 2, but since the first row will be dropped, and the last row is simply 0, so we omit it.
    
    gsl_vector_view lastRow=gsl_matrix_row(dctr, Nr-1);
    gsl_vector_set_zero(&lastRow.vector);
    
    for (int iterr=Nr-3; iterr>0; iterr-=2)
    {
        gsl_vector_view nextLastRow=gsl_matrix_row(dctr, iterr+1);
        lastRow=gsl_matrix_row(dctr,iterr+2);
        gsl_vector_view cRow=gsl_matrix_row(dctr, iterr);
        
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iter=0; iter<matrixW; ++iter)
        {
            cRow.vector.data[iter]=nextLastRow.vector.data[iter]*2.0*(iterr+1)/logicNr+lastRow.vector.data[iter];
            nextLastRow.vector.data[iter]=0;
        }
    }
    lastRow=gsl_matrix_row(dctr, 0);
    gsl_vector_set_zero(&lastRow.vector);
    // aliasing
//    int aliasing=floor(Nr*0.6666666666);
//    for (int itert=0; itert<matrixW; ++itert)
//    {
//        for (int iterr=aliasing; iterr<Nr; ++iterr)
//        {
//            gsl_matrix_set(dctr, iterr, itert, 0);
//        }
//    }
    
    //gsl_matrix_view middle=gsl_matrix_submatrix(dctr, 1, 0, Nr-2, matrixW);
    
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iter=1; iter<Nr-2; iter+=2)
    {
        gsl_vector_view temp=gsl_matrix_row(dctr, iter);
        gsl_vector_scale(&temp.vector, 0.5);
    }
    //gsl_matrix_scale(&middle.matrix, 0.5);
    
    fftw_execute(dctr2r);
    
    
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        gsl_matrix_view datablock=gsl_matrix_submatrix(dFields, iterf*Nrp, 0, Nrp, Ntheta);
        gsl_matrix_view destiny=gsl_matrix_submatrix(dctr, 0, iterf*Ntheta, Nrp, Ntheta);
        gsl_matrix_memcpy(&datablock.matrix, &destiny.matrix);
    }
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
    
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iterr=0; iterr<matrixH; ++iterr)
    {
        // aliasing
        int aliasing=floor(Ntheta/3);
        for (int itertheta=0; itertheta<aliasing; ++itertheta)
        {
            gsl_complex temp=gsl_matrix_complex_get(fftc, iterr, itertheta);
            //a'=a*ik
            double swap=temp.dat[0];
            temp.dat[0]=-temp.dat[1]*itertheta/Ntheta;
            temp.dat[1]=swap*itertheta/Ntheta;
            gsl_matrix_complex_set(fftc, iterr, itertheta, temp);
        }
        for (int itertheta=aliasing; itertheta<Ntheta/2+1; ++itertheta)
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