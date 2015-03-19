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

void solver::dr(bool ifFirst)
{
    
//#ifdef MULTIPROCESS
//    omp_set_num_threads(8);
//#endif
    
    gsl_matrix* cFields;
    if (ifFirst)
    {
        cFields=Fields;
    }
    else
    {
        cFields=tempFields;
    }
    
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        gsl_matrix_view datablock=gsl_matrix_submatrix(cFields, iterf*Nrp, 0, Nrp, Ntheta);
        gsl_matrix_view destiny=gsl_matrix_submatrix(dctr, 0, iterf*Ntheta, Nrp, Ntheta);
        gsl_matrix_memcpy(&destiny.matrix, &datablock.matrix);
    }
    
    for (int iterr=0; iterr<Nrp; ++iterr)
    {
        gsl_vector_view datablock=gsl_matrix_row(dctr, iterr);
        gsl_vector_view destiny=gsl_matrix_row(dctr, Nr-iterr-1);
        gsl_vector_memcpy(&destiny.vector, &datablock.vector);
    }
    
//    #ifdef MULTIPROCESS
//    #pragma omp parallel for
//    #endif
//    for (int iter=0; iter<NumPoints*NumField; ++iter)
//    {
//        int iterx=iter%Ntheta;
//        int itery=(iter-iterx)/Ntheta;
//        int iteryp=itery%Nrp;
//        int iterf=(itery-iteryp)/Nrp;
//        iterx=iterx+iterf*Ntheta;
//        iteryp=iteryp*matrixW;
//        dctr->data[iterx+iteryp]=cFields->data[iter];
//        dctr->data[iterx+logicNr*matrixW-iteryp]=cFields->data[iter];
//    }
    
    fftw_execute(dctr2r);
    //printdebugM(dctr, "dctr.txt");
    //The first and last row should divide 2, but since the first row will be dropped, and the last row is simply 0, so we omit it.
    
    // aliasing
    for (int itert=0; itert<matrixW; ++itert)
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
        
        for (int iter=0; iter<matrixW; ++iter)
        {
            cRow.vector.data[iter]=nextLastRow.vector.data[iter]*2.0*(iterr+1)/logicNr+lastRow.vector.data[iter];
            nextLastRow.vector.data[iter]=0;
        }
    }
    lastRow=gsl_matrix_row(dctr, 0);
    gsl_vector_set_zero(&lastRow.vector);
    
    
//#ifdef MULTIPROCESS
//#pragma omp parallel for
//#endif
    for (int iter=1; iter<aliasingr-2; iter+=2)
    {
        gsl_vector_view temp=gsl_matrix_row(dctr, iter);
        gsl_vector_scale(&temp.vector, 0.5);
    }
    //gsl_matrix_scale(&middle.matrix, 0.5);
    
    fftw_execute(dctr2r);
    
    
//#ifdef MULTIPROCESS
//#pragma omp parallel for
//#endif
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        gsl_matrix_view datablock=gsl_matrix_submatrix(dFields, iterf*Nrp, 0, Nrp, Ntheta);
        gsl_matrix_view destiny=gsl_matrix_submatrix(dctr, 0, iterf*Ntheta, Nrp, Ntheta);
        gsl_matrix_memcpy(&datablock.matrix, &destiny.matrix);
    }
}

void solver::dtheta(bool ifFirst)
{
    
//#ifdef MULTIPROCESS
//    omp_set_num_threads(8);
//#endif
    
    if (ifFirst)
    {
        fftw_execute(fftr2c);
    }
    else
    {
        fftw_execute(tempfftr2c);
    }
    
//#ifdef MULTIPROCESS
//#pragma omp parallel for
//#endif
    for (int iterr=0; iterr<matrixH; ++iterr)
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