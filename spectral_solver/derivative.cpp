//
//  derivative.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
#include <iostream>
using namespace std;

void solver::dr(bool ifFirst)
{
    matrix<long double>* cFields;
    if (ifFirst)
    {
        cFields=&Fields;
    }
    else
    {
        cFields=&tempFields;
    }
    
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        dctr.blockCopy(0, Nrp, iterf*Ntheta, iterf*Ntheta+Ntheta, *cFields, iterf*Nrp, 0);
    }
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iterr=0; iterr<Nrp; ++iterr)
    {
        for (int itert=0; itert<matrixW; ++itert)
        {
            dctr.ele(Nr-iterr-1, itert)=dctr.ele(iterr, itert);
        }
    }
    
    fftwl_execute(dctr2r);
    //printdebugM(dctr, "dctr.txt");
    //The first and last row should divide 2, but since the first row will be dropped, and the last row is simply 0, so we omit it.
    
    // aliasing
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iterr=aliasingr-1; iterr<Nr; ++iterr)
    {
        dctr.setRow(iterr, 0);
    }
    
    //gsl_vector_view lastRow; //=matrix<long double>_row(dctr, aliasingr-1);
    //gsl_vector_set_zero(&lastRow.vector);
    
    for (int iterr=aliasingr-3; iterr>0; iterr-=2)
    {
        //gsl_vector_view nextLastRow=matrix<long double>_row(dctr, iterr+1);
        //lastRow=matrix<long double>_row(dctr,iterr+2);
        //gsl_vector_view cRow=matrix<long double>_row(dctr, iterr);
        
        for (int iter=0; iter<matrixW; ++iter)
        {
            //cRow.vector.data[iter]=nextLastRow.vector.data[iter]*2.0*(iterr+1)/logicNr+lastRow.vector.data[iter];
            //nextLastRow.vector.data[iter]=0;
            dctr.ele(iterr,iter)=dctr.ele(iterr+1,iter)*2.0*(iterr+1)/logicNr+dctr.ele(iterr+2,iter);
        }
        dctr.setRow(iterr+1, 0);
    }
    //lastRow=matrix<long double>_row(dctr, 0);
    //gsl_vector_set_zero(&lastRow.vector);
    dctr.setRow(0, 0);
    
    
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iter=1; iter<aliasingr-2; iter+=2)
    {
        //gsl_vector_view temp=matrix<long double>_row(dctr, iter);
        //gsl_vector_scale(&temp.vector, 0.5);
        dctr.scaleRow(iter, 0.5);
    }
    
    fftwl_execute(dctr2r);
    
    
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        dFields.blockCopy(iterf*Nrp, iterf*Nrp+Nrp, 0, Ntheta, dctr, 0, iterf*Ntheta);
    }
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
    
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    // dealiasing
    for (int itertheta=0; itertheta<aliasingt; ++itertheta)
    {
        fftc.scaleCol(itertheta, 0, (double)itertheta/(double)Ntheta);
    }
    for (int itertheta=aliasingt; itertheta<Ntheta/2+1; ++itertheta)
    {
        fftc.setCol(itertheta, 0, 0);
    }
    
    if (ifFirst)
    {
        fftwl_execute(ifftc2r);
    }
    else
    {
        fftwl_execute(tempifftc2r);
    }
}

void solver::drWOA()
{
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        dctr.blockCopy(0, Nrp, iterf*Ntheta, iterf*Ntheta+Ntheta, Fields, iterf*Nrp, 0);
    }
    
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iterr=0; iterr<Nrp; ++iterr)
    {
        for (int iterr=0; iterr<Nrp; ++iterr)
        {
            for (int itert=0; itert<matrixW; ++itert)
            {
                dctr.ele(Nr-iterr-1, itert)=dctr.ele(iterr, itert);
            }
        }
    }
    
    //#ifdef MULTIPROCESS
    //#pragma omp parallel for
    //#endif
    //    for (int iter=0; iter<totalPoints; ++iter)
    //    {
    //        int iterx=iter%Ntheta;
    //        int itery=(iter-iterx)/Ntheta;
    //        int iteryp=itery%Nrp;
    //        int iterf=(itery-iteryp)/Nrp;
    //        iterx=iterx+iterf*Ntheta;
    //        iteryp=iteryp*matrixW;
    //        dctr[iterx+iteryp]=cFields[iter];
    //        dctr[iterx+logicNr*matrixW-iteryp]=cFields[iter];
    //    }
    
    fftwl_execute(dctr2r);
    //printdebugM(dctr, "dctr.txt");
    //The first and last row should divide 2, but since the first row will be dropped, and the last row is simply 0, so we omit it.
    
    //gsl_vector_view lastRow=matrix<long double>_row(dctr, Nr-1);
    //gsl_vector_set_zero(&lastRow.vector);
    dctr.setRow(Nr-1, 0);
    
    for (int iterr=Nr-3; iterr>0; iterr-=2)
    {
        for (int iter=0; iter<matrixW; ++iter)
        {
            dctr.ele(iterr,iter)=dctr.ele(iterr+1,iter)*2.0*(iterr+1)/logicNr+dctr.ele(iterr+2,iter);
        }
        dctr.setRow(iterr+1, 0);
    }
    dctr.setRow(0, 0);
    
    
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iter=1; iter<Nr-2; iter+=2)
    {
        dctr.scaleRow(iter, 0.5);
    }
    
    fftwl_execute(dctr2r);
    
    
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        dFields.blockCopy(iterf*Nrp, iterf*Nrp+Nrp, 0, Ntheta, dctr, 0, iterf*Ntheta);
    }
}