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
#include <mpi.h>

using namespace std;

class solver
{
public:
    solver();
    ~solver();
    
    void Fun(gsl_matrix *result);
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
    void printdebugM(gsl_matrix* m, const string filename);
    void printdebugCM(gsl_matrix_complex* m, const string filename);
    void readFile(const string filename);
    
    gsl_matrix *Fields;
    
    gsl_matrix *k1;
    gsl_matrix *k2;
    gsl_matrix *k3;
    gsl_matrix *k4;
    gsl_matrix *k5;
    gsl_matrix *k6;
    gsl_matrix *k7;
    gsl_matrix *k8;
    gsl_matrix *odetempField;
    gsl_matrix *odetempField2;
    gsl_matrix *odetempField3;
    
    vector <gsl_matrix*> HistoryFields;
    vector <gsl_matrix*> Hij;
    gsl_matrix *G;
    gsl_vector *r;
    gsl_vector *r2;
    gsl_vector *theta;
    gsl_vector *boundary;
    
    gsl_matrix *FieldsLocal;
    gsl_matrix *tempFieldsLocal;
    gsl_matrix *GLocal;
    vector <gsl_matrix*> HijLocal;
    gsl_matrix *dFieldsLocal;
    gsl_matrix_view dFieldLocalView[NumField];
    gsl_matrix_view HijLocalView[NumField*NumField];
    gsl_matrix *iterFieldsLocal;
    
    fftw_plan fftr2c;
    fftw_plan ifftc2r;
    fftw_plan tempfftr2c;
    fftw_plan tempifftc2r;
    fftw_plan dctr2r;
    fftw_plan tempdctr2r;
    
    gsl_matrix_complex *fftc;
    gsl_matrix *dctr;
    gsl_matrix *tempdctr;
    
    double time;
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
