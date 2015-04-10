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

//#include <gsl/matrix<long double>.h>
//#include <gsl/gsl_vector.h>
#include "matrix.h"
#include <vector>
#include <fftw3.h>
#include "parameters.h"
#include <string>
#include <fstream>
using namespace std;

#ifdef MULTIPROCESS
#include <omp.h>
#endif

using namespace std;

class solver
{
public:
    solver();
    ~solver();
    
    void Fun(matrix<long double>& result);
    void HGFuns();
    void dr(bool ifFirst);
    void drWOA();
    void dtheta(bool ifFirst);
    
    void RK4Step();
    void RK6Step();
    void BDF3Step();
    void BDF4Step();
    void BDF5Step();
    void BDF6Step();
    void setBoundary();
    void initialization();
    void solve(int totaliter);
    void DoubleTimeStepBDF4();
    void DoubleTimeStepBDF6();
    
    void printstatus();
    void printdebugM(matrix<long double>* m, const string filename);
    void printdebugCM(matrix_complex<long double>* m, const string filename);
    void readFile(const string filename);
    
    matrix<long double> Fields;
    matrix<long double> dFields;
    matrix<long double> tempFields;
    matrix<long double> caltempFields;
    matrix<long double> boundary;
    
    matrix<long double> k1;
    matrix<long double> k2;
    matrix<long double> k3;
    matrix<long double> k4;
    matrix<long double> k5;
    matrix<long double> k6;
    matrix<long double> k7;
    matrix<long double> k8;
    matrix<long double> odetempField;
    matrix<long double> odetempField2;
    
    vector <matrix<long double>*> HistoryFields;
    vector <matrix<long double>*> DoubledHistoryFields;
    vector <matrix<long double>*> Hij;
    matrix<long double>G;
    math_vector<long double> r;
    math_vector<long double> r2;
    math_vector<long double> theta;
    
    fftwl_plan fftr2c;
    fftwl_plan ifftc2r;
    fftwl_plan tempfftr2c;
    fftwl_plan tempifftc2r;
    fftwl_plan dctr2r;
    
    matrix_complex<long double> fftc;
    matrix<long double>dctr;
    
    //math_vector<long double> tempstore;
    //math_vector<long double> tempstore2;
    
    long double StepT;
    long double time;
    int timeIdx;
    ofstream timefile;
};

#endif
