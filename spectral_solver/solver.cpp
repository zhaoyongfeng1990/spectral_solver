//
//  solver.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
#include <cmath>

solver::solver()
{
    time=0;
    timeIdx=0;
    
    Fields=gsl_matrix_calloc(matrixH, Ntheta);
    gsl_matrix_set_zero(Fields);
    dFields=gsl_matrix_calloc(matrixH, Ntheta);
    gsl_matrix_set_zero(dFields);
    tempFields=gsl_matrix_calloc(matrixH, Ntheta);
    gsl_matrix_set_zero(tempFields);
    caltempFields=gsl_matrix_calloc(matrixH, Ntheta);
    gsl_matrix_set_zero(caltempFields);
    G=gsl_matrix_calloc(matrixH, Ntheta);
    gsl_matrix_set_zero(G);
    boundary=gsl_matrix_calloc(NumField, Ntheta);
    
    for (int iter=0; iter<NumField; ++iter)
    {
        dFieldView[iter]=gsl_matrix_submatrix(dFields, iter*Nrp, 0, Nrp, Ntheta);
        ctFieldView[iter]=gsl_matrix_submatrix(caltempFields, iter*Nrp, 0, Nrp, Ntheta);
    }
    
    r=gsl_vector_calloc(Nrp);
    for (int iter=0; iter<Nrp; ++iter)
    {
        r->data[iter]=cos(iter*PI/logicNr);
    }
    
    theta=gsl_vector_calloc(Ntheta);
    for (int iter=0; iter<Ntheta; ++iter)
    {
        theta->data[iter]=2*PI*iter/Ntheta;
    }
    
    HistoryFields.resize(3);
    for (int iterh=0; iterh<3; ++iterh)
    {
        HistoryFields[iterh]=gsl_matrix_calloc(matrixH, Ntheta);
        gsl_matrix_set_zero(HistoryFields[iterh]);
    }
    
    Hij.resize(NumField);
    for (int iterh=0; iterh<NumField*NumField; ++iterh)
    {
        Hij[iterh]=gsl_matrix_calloc(matrixH, Ntheta);
        gsl_matrix_set_zero(Hij[iterh]);
    }
    
    int n1[]={Nr};
    int n2[]={Ntheta};
    
    dctr=gsl_matrix_calloc(Nr, Ntheta*NumField);
    fftc=gsl_matrix_complex_calloc(matrixH, Ntheta/2+1);
    
    fftr2c=fftw_plan_many_dft_r2c(1, n2, matrixH, Fields->data, n2, 1, Ntheta, (fftw_complex *)fftc->data, n2, 1, Ntheta/2+1, FFTW_MEASURE);
    tempfftr2c=fftw_plan_many_dft_r2c(1, n2, matrixH, tempFields->data, n2, 1, Ntheta, (fftw_complex *)fftc->data, n2, 1, Ntheta/2+1, FFTW_MEASURE);
    
    ifftc2r=fftw_plan_many_dft_c2r(1, n2, matrixH, (fftw_complex *)fftc->data, n2, 1, Ntheta/2+1, dFields->data, n2, 1, Ntheta, FFTW_MEASURE);
    tempifftc2r=fftw_plan_many_dft_c2r(1, n2, matrixH, (fftw_complex *)fftc->data, n2, 1, Ntheta/2+1, dFields->data, n2, 1, Ntheta, FFTW_MEASURE);
    
    fftw_r2r_kind kind[]={FFTW_REDFT00};
    
    dctr2r=fftw_plan_many_r2r(1, n1, matrixW, dctr->data, n1, matrixW, 1, dctr->data, n1, matrixW, 1, kind, FFTW_MEASURE);
    
    tempstore=gsl_vector_alloc(matrixW);
    tempstore2=gsl_vector_alloc(matrixW);
}

solver::~solver()
{
    delete [] r;
    delete [] theta;
    gsl_matrix_free(Fields);
    gsl_matrix_free(dFields);
    gsl_matrix_free(tempFields);
    gsl_matrix_free(caltempFields);
    gsl_matrix_free(G);
    gsl_matrix_free(dctr);
    gsl_matrix_complex_free(fftc);
    gsl_matrix_free(boundary);
    
    fftw_destroy_plan(fftr2c);
    fftw_destroy_plan(ifftc2r);
    fftw_destroy_plan(tempfftr2c);
    fftw_destroy_plan(tempifftc2r);
    fftw_destroy_plan(dctr2r);
    
    for (int iterh=0; iterh<3; ++iterh)
    {
        gsl_matrix_free(HistoryFields[iterh]);
    }
    
    for (int iterh=0; iterh<NumField; ++iterh)
    {
        gsl_matrix_free(Hij[iterh]);
    }
    
    gsl_vector_free(tempstore);
    gsl_vector_free(tempstore2);
}
