//
//  solver.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
#include <cmath>

solver::solver() : timefile("time.txt")
{
    fftwl_init_threads();
    fftwl_plan_with_nthreads(4);
    
    time=0;
    StepT=iniStepT;
    timeIdx=0;
    
    Fields.alloc(matrixH, Ntheta);
    dFields.alloc(matrixH, Ntheta);
    tempFields.alloc(matrixH, Ntheta);
    caltempFields.alloc(matrixH, Ntheta);
    G.alloc(matrixH, Ntheta);
    
    k1.alloc(matrixH, Ntheta);
    k2.alloc(matrixH, Ntheta);
    k3.alloc(matrixH, Ntheta);
    k4.alloc(matrixH, Ntheta);
    k5.alloc(matrixH, Ntheta);
    k6.alloc(matrixH, Ntheta);
    k7.alloc(matrixH, Ntheta);
    k8.alloc(matrixH, Ntheta);
    
    odetempField.alloc(matrixH, Ntheta);
    odetempField2.alloc(matrixH, Ntheta);
    
    boundary.alloc(NumField, Ntheta);
    
    //for (int iter=0; iter<NumField; ++iter)
    //{
    //    dFieldView[iter]=matrix<long double>_submatrix(dFields, iter*Nrp, 0, Nrp, Ntheta);
    //    ctFieldView[iter]=matrix<long double>_submatrix(caltempFields, iter*Nrp, 0, Nrp, Ntheta);
    //}
    
    r.alloc(Nrp);
    r2.alloc(Nrp);
    for (int iter=0; iter<Nrp; ++iter)
    {
        r[iter]=cos(iter*PI/logicNr);
        r2[iter]=r[iter]*r[iter];
    }
    
    theta.alloc(Ntheta);
    for (int iter=0; iter<Ntheta; ++iter)
    {
        theta[iter]=2*PI*iter/Ntheta;
    }
    
    HistoryFields.resize(5);
    for (int iterh=0; iterh<5; ++iterh)
    {
        HistoryFields[iterh]=new matrix<long double>();
        HistoryFields[iterh]->alloc(matrixH, Ntheta);
    }
    
    if (IncreaseTimes!=0)
    {
        DoubledHistoryFields.resize(3);
        for (int iterh=0; iterh<3; ++iterh)
        {
            DoubledHistoryFields[iterh]=new matrix<long double>();
            DoubledHistoryFields[iterh]->alloc(matrixH, Ntheta);\
        }
    }
    
    Hij.resize(NumField);
    for (int iterh=0; iterh<NumField; ++iterh)
    {
        Hij[iterh]=new matrix<long double>();
        Hij[iterh]->alloc(matrixH, Ntheta);
    }
    
    int n1[]={Nr};
    int n2[]={Ntheta};
    
    dctr.alloc(Nr, Ntheta*NumField);
    fftc.alloc(matrixH, Ntheta/2+1);
    
    fftr2c=fftwl_plan_many_dft_r2c(1, n2, matrixH, Fields.data, n2, 1, Ntheta, (fftwl_complex *)fftc.data, n2, 1, Ntheta/2+1, FFTW_MEASURE);
    tempfftr2c=fftwl_plan_many_dft_r2c(1, n2, matrixH, tempFields.data, n2, 1, Ntheta, (fftwl_complex *)fftc.data, n2, 1, Ntheta/2+1, FFTW_MEASURE);
    
    ifftc2r=fftwl_plan_many_dft_c2r(1, n2, matrixH, (fftwl_complex *)fftc.data, n2, 1, Ntheta/2+1, dFields.data, n2, 1, Ntheta, FFTW_MEASURE);
    tempifftc2r=fftwl_plan_many_dft_c2r(1, n2, matrixH, (fftwl_complex *)fftc.data, n2, 1, Ntheta/2+1, dFields.data, n2, 1, Ntheta, FFTW_MEASURE);
    
    fftwl_r2r_kind kind[]={FFTW_REDFT00};
    
    dctr2r=fftwl_plan_many_r2r(1, n1, matrixW, dctr.data, n1, matrixW, 1, dctr.data, n1, matrixW, 1, kind, FFTW_MEASURE);
    
    
    //tempstore=gsl_vector_alloc(matrixW);
    //tempstore2=gsl_vector_alloc(matrixW);
}

solver::~solver()
{
    fftwl_cleanup_threads();
    
    fftwl_destroy_plan(fftr2c);
    fftwl_destroy_plan(ifftc2r);
    fftwl_destroy_plan(tempfftr2c);
    fftwl_destroy_plan(tempifftc2r);
    fftwl_destroy_plan(dctr2r);
    
    for (int iterh=0; iterh<5; ++iterh)
    {
        delete HistoryFields[iterh];
    }
    if (IncreaseTimes!=0)
    {
        for (int iterh=0; iterh<3; ++iterh)
        {
            delete DoubledHistoryFields[iterh];
        }
    }
    
    for (int iterh=0; iterh<NumField; ++iterh)
    {
        delete Hij[iterh];
    }
    
    //gsl_vector_free(tempstore);
    //gsl_vector_free(tempstore2);
}
