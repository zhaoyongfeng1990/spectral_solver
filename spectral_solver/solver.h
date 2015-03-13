//
//  solver.h
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#ifndef spectral_solver_solver_h
#define spectral_solver_solver_h

#define HAVE_INLINE

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <vector>
#include <fftw3.h>
#include "parameters.h"
#include <string>

#include <omp.h>

using namespace std;

class solver
{
public:
    solver();
    ~solver();
    
    void Fun(gsl_matrix *result);
    void HGFuns();
    void dr(bool ifFirst);
    void dtheta(bool ifFirst);
    
    void RK4Step();
    void BDF4Step();
    void setBoundary();
    void initialization();
    void solve(int totaliter);
    
    void printstatus();
    void printdebugM(gsl_matrix* m, const string filename);
    
    gsl_matrix *Fields;
    gsl_matrix *dFields;
    gsl_matrix *tempFields;
    gsl_matrix *caltempFields;
    gsl_matrix *boundary;
    
    gsl_matrix_view dFieldView[NumField];
    gsl_matrix_view ctFieldView[NumField];
    
    gsl_matrix *k1;
    gsl_matrix *k2;
    gsl_matrix *k3;
    gsl_matrix *k4;
    gsl_matrix *odetempField;
    
    vector <gsl_matrix*> HistoryFields;
    vector <gsl_matrix*> Hij;
    gsl_matrix *G;
    gsl_vector *r;
    gsl_vector *r2;
    gsl_vector *theta;
    
    fftw_plan fftr2c;
    fftw_plan ifftc2r;
    fftw_plan tempfftr2c;
    fftw_plan tempifftc2r;
    fftw_plan dctr2r;
    
    gsl_matrix_complex *fftc;
    gsl_matrix *dctr;
    //gsl_matrix *odctr;
    
    gsl_vector* tempstore;
    gsl_vector* tempstore2;
    
    double time;
    int timeIdx;
    
//    //For Broyden method
//    gsl_matrix *A;
//    gsl_matrix *Ainv;
//    gsl_matrix *xk_0;
//    gsl_matrix *xk_1;
//    gsl_matrix *dx;
//    gsl_matrix *fxk_0;
//    gsl_matrix *fxk_1;
//    gsl_matrix *Adxk;
//    gsl_matrix *temp;
//    gsl_matrix *C;
//    gsl_matrix *temp1;
//    gsl_matrix *temp2;
//    gsl_matrix *temp3;
};

#endif
