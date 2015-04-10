//
//  solver.h
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#ifndef spectral_solver_solver_h
#define spectral_solver_solver_h

#include <vector>
#include <fftw3.h>
#include "parameters.h"
#include <string>
#include <mpi.h>
#include "matrix.h"

using namespace std;

class solver
{
public:
    solver();
    ~solver();
    
    void Fun(matrix<long double> &result);
    void HGFuns();
    void HFunsForR();
    void dr(bool ifFirst);
    void dtheta(bool ifFirst);
    void drWOA();
    
    void RK4Step();
    void RK6Step();
    void BDF4Step();
    void BDF6Step();
    void setBoundary();
    void initialization();
    void solve(int totaliter);
    
    void printstatus();
    void printdebugM(matrix<long double> *m, const string filename);
    void printdebugCM(matrix_complex<long double> *m, const string filename);
    void readFile(const string filename);
    
    matrix<long double> Fields;
    
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
    matrix<long double> odetempField3;
    
    vector <matrix<long double>*> HistoryFields;
    vector <matrix<long double>*> Hij;
    matrix<long double> G;
    math_vector<long double> r;
    math_vector<long double> r2;
    math_vector<long double> theta;
    math_vector<long double> boundary;
    
    matrix<long double>FieldsLocal;
    matrix<long double>tempFieldsLocal;
    matrix<long double>GLocal;
    vector <matrix<long double>*> HijLocal;
    matrix<long double>dFieldsLocal;
    //gsl_matrix_view dFieldLocalView[NumField];
    //gsl_matrix_view HijLocalView[NumField*NumField];
    matrix<long double>iterFieldsLocal;
    
    fftwl_plan fftr2c;
    fftwl_plan ifftc2r;
    fftwl_plan tempfftr2c;
    fftwl_plan tempifftc2r;
    fftwl_plan dctr2r;
    fftwl_plan tempdctr2r;
    
    matrix_complex<long double> fftc;
    matrix<long double>dctr;
    matrix<long double>tempdctr;
    
    long double time;
    int timeIdx;
    
    int workerR;
    int workerT;
    int workerTl;
    int workerRl;
    int bossR;
    int bossT;
    int bossTl;
    int bossRl;
    
    int workerPointsR;
    int workerPointsT;
    int bossPointsR;
    int bossPointsT;
    
    int jobR;
    int jobT;
    int jobTl;
    int jobRl;
    int jobPointsRl;
    int jobPointsTl;
    int jobPointsR;
    int jobPointsT;
    
    int numOfProcess;
    int numOfProcessT;
    int numOfProcessR;
    int cRank;
    
    int bossBD4;
    int workerBD4;
    int bossP;
    int workerP;
    int jobBD4;
    int iterPoints;
    
    MPI_Datatype* RblockType;
    MPI_Datatype* TblockType;
    MPI_Datatype* BoundaryType;
    MPI_Datatype* BD4Type;
    MPI_Status status;
    //MPI_Request request1;
    //MPI_Request request2;
};

#endif
