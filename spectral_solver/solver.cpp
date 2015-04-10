//
//  solver.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
#include <cmath>
#include <iostream>
using namespace std;

solver::solver()
{
    MPI_Comm_rank(MPI_COMM_WORLD,& cRank);
    MPI_Comm_size(MPI_COMM_WORLD,& numOfProcess);
    
    numOfProcessT=numOfProcess/2; //Num Of Processes should be even!
    numOfProcessR=numOfProcess-numOfProcessT;
    
    workerRl=ceil(Nrp/(long double)numOfProcessR);
    workerTl=ceil(Ntheta/(long double)numOfProcessT);
    workerT=workerTl*NumField;
    workerR=workerRl*NumField;
    workerPointsR=Ntheta*workerRl;
    workerPointsT=workerTl*Nrp;
    
    bossRl=Nrp-workerRl*(numOfProcessR-1);
    bossTl=Ntheta-workerTl*(numOfProcessT-1);
    bossT=bossTl*NumField;
    bossR=bossRl*NumField;
    bossPointsR=Ntheta*bossRl;
    bossPointsT=bossTl*Nrp;
    
    workerBD4=floor(Nrp/(long double)numOfProcess);
    bossBD4=Nrp-workerBD4*(numOfProcess-1);
    workerP=workerBD4*Ntheta;
    bossP=bossBD4*Ntheta;
    
    if (0==cRank)
    {
        jobRl=bossRl;
        jobTl=bossTl;
        jobR=bossR;
        jobT=bossT;
        jobPointsRl=bossPointsR;
        jobPointsTl=bossPointsT;
        
        jobBD4=bossBD4;
        iterPoints=bossP*NumField;
    }
    else if (1==cRank)
    {
        jobRl=bossRl;
        jobTl=bossTl;
        jobR=bossR;
        jobT=bossT;
        jobPointsRl=bossPointsR;
        jobPointsTl=bossPointsT;
        
        jobBD4=workerBD4;
        iterPoints=workerP*NumField;
    }
    else
    {
        jobRl=workerRl;
        jobTl=workerTl;
        jobR=workerR;
        jobT=workerT;
        jobPointsRl=workerPointsR;
        jobPointsTl=workerPointsT;
        
        jobBD4=workerBD4;
        iterPoints=workerP*NumField;
    }
    
    jobPointsR=jobPointsRl*NumField;
    jobPointsT=jobPointsTl*NumField;
    
    time=0;
    timeIdx=0;
    
    odetempField2.alloc(jobBD4*NumField, Ntheta);
    odetempField3.alloc(jobBD4*NumField, Ntheta);
    iterFieldsLocal.alloc(jobBD4*NumField, Ntheta);
    
    HistoryFields.resize(3);
    for (int iterh=0; iterh<3; ++iterh)
    {
        HistoryFields[iterh]=new matrix<long double>(jobBD4*NumField, Ntheta);
    }
    
    if (cRank==0)
    {
        Fields.alloc(matrixH, Ntheta);
        G.alloc(matrixH, Ntheta);
        
        k1.alloc(matrixH, Ntheta);
        k2.alloc(matrixH, Ntheta);
        k3.alloc(matrixH, Ntheta);
        k4.alloc(matrixH, Ntheta);
        //k5=gsl_matrix_alloc(matrixH, Ntheta);
        //k6=gsl_matrix_alloc(matrixH, Ntheta);
        //k7=gsl_matrix_alloc(matrixH, Ntheta);
        //k8=gsl_matrix_alloc(matrixH, Ntheta);
        
        odetempField.alloc(matrixH, Ntheta);
        
        
        
        Hij.resize(NumField);
        for (int iterh=0; iterh<NumField; ++iterh)
        {
            Hij[iterh]=new matrix<long double>(matrixH, Ntheta);
        }
    }
    
    //local variables
    if (cRank%2==0)
    {
        FieldsLocal.alloc(jobR, Ntheta);
        tempFieldsLocal.alloc(jobR, Ntheta);
        dFieldsLocal.alloc(jobR, Ntheta);
        
        HijLocal.resize(NumField);
        for (int iterh=0; iterh<NumField; ++iterh)
        {
            HijLocal[iterh]=new matrix<long double>(jobR, Ntheta);
        }
        
//        for (int iter=0; iter<NumField; ++iter)
//        {
//            dFieldLocalView[iter]=gsl_matrix_submatrix(dFieldsLocal, iter*jobRl, 0, jobRl, Ntheta);
//            for (int iterf=0; iterf<NumField; ++iterf)
//            {
//                HijLocalView[iterf*NumField+iter]=gsl_matrix_submatrix(HijLocal[iterf], iter*jobRl, 0, jobRl, Ntheta);
//            }
//        }
        
        GLocal.alloc(jobR, Ntheta);
        
        int n2[]={Ntheta};
        fftc.alloc(jobR, Ntheta/2+1);
        
        fftr2c=fftwl_plan_many_dft_r2c(1, n2, jobR, FieldsLocal.data, n2, 1, Ntheta, (fftwl_complex *)fftc.data, n2, 1, Ntheta/2+1, FFTW_MEASURE);
        ifftc2r=fftwl_plan_many_dft_c2r(1, n2, jobR, (fftwl_complex *)fftc.data, n2, 1, Ntheta/2+1, dFieldsLocal.data, n2, 1, Ntheta, FFTW_MEASURE);
        
        //cout << "test" << endl;

        tempfftr2c=fftwl_plan_many_dft_r2c(1, n2, jobR, tempFieldsLocal.data, n2, 1, Ntheta, (fftwl_complex *)fftc.data, n2, 1, Ntheta/2+1, FFTW_MEASURE);
        
        tempifftc2r=fftwl_plan_many_dft_c2r(1, n2, jobR, (fftwl_complex *)fftc.data, n2, 1, Ntheta/2+1, dFieldsLocal.data, n2, 1, Ntheta, FFTW_MEASURE);
        
    }
    else
    {
        HijLocal.resize(NumField);
        for (int iterh=0; iterh<NumField; ++iterh)
        {
            HijLocal[iterh]=new matrix<long double>(Nrp, jobT);
        }
        
        tempFieldsLocal.alloc(Nrp, jobT);
        dctr.alloc(Nr, jobT);
        tempdctr.alloc(Nr, jobT);
        boundary.alloc(jobT);
        
//        for (int iter=0; iter<NumField; ++iter)
//        {
//            dFieldLocalView[iter]=gsl_matrix_submatrix(dctr, 0, iter*jobTl, Nrp, jobTl);
//            for (int iterf=0; iterf<NumField; ++iterf)
//            {
//                HijLocalView[iterf*NumField+iter]=gsl_matrix_submatrix(HijLocal[iterf], 0, iter*jobTl, Nrp, jobTl);
//            }
//        }
        
        int n1[]={Nr};
        fftwl_r2r_kind kind[]={FFTW_REDFT00};
        dctr2r=fftwl_plan_many_r2r(1, n1, jobT, dctr.data, n1, jobT, 1, dctr.data, n1, jobT, 1, kind, FFTW_MEASURE);
        tempdctr2r=fftwl_plan_many_r2r(1, n1, jobT, tempdctr.data, n1, jobT, 1, dctr.data, n1, jobT, 1, kind, FFTW_MEASURE);
    }
    
    r.alloc(Nrp);
    r2.alloc(Nrp);
    for (int iter=0; iter<Nrp; ++iter)
    {
        r[iter]=cos(iter*PI/logicNr);
        r2[iter]=cos(iter*PI/logicNr)*cos(iter*PI/logicNr);
    }
    
    theta.alloc(Ntheta);
    for (int iter=0; iter<Ntheta; ++iter)
    {
        theta[iter]=2*PI*iter/Ntheta;
    }
    
    //MPI Data Type
    int blocklengthsR[NumField];
    int blocklengthsT[matrixH];
    int blockpossR[NumField];
    int blockpossT[matrixH];
    
    RblockType=new MPI_Datatype[numOfProcessR];
    TblockType=new MPI_Datatype[numOfProcessT];
    BoundaryType=new MPI_Datatype[numOfProcessT];
    
    //process 0 - boss for theta calculation (split r)
    for (int iter=0; iter<NumField; ++iter)
    {
        blocklengthsR[iter]=bossPointsR;
        blockpossR[iter]=iter*NumPoints;
    }
    MPI_Type_indexed(NumField, blocklengthsR, blockpossR, MPI_LONG_DOUBLE,& RblockType[0]);
    MPI_Type_commit(&RblockType[0]);
    
    //process 1 - boss for r calculation (split theta)
    int itertemp=0;
    for (int iterr=0; iterr<Nrp; ++iterr)
    {
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            blocklengthsT[itertemp]=bossTl;
            blockpossT[itertemp]=iterr*Ntheta+iterf*NumPoints;
            ++itertemp;
        }
    }
    MPI_Type_indexed(matrixH, blocklengthsT, blockpossT, MPI_LONG_DOUBLE,& TblockType[0]);
    MPI_Type_commit(&TblockType[0]);
    MPI_Type_indexed(NumField, blocklengthsT, blockpossT, MPI_LONG_DOUBLE,& BoundaryType[0]);
    MPI_Type_commit(&BoundaryType[0]);
    
    //process 2, 4, 6, ... - workers for theta calculation (split r)
    for (int iterCPU=1; iterCPU<numOfProcessR; ++iterCPU)
    {
        for (int iter=0; iter<NumField; ++iter)
        {
            blocklengthsR[iter]=workerPointsR;
            blockpossR[iter]=bossPointsR+(iterCPU-1)*workerPointsR+iter*NumPoints;
        }
        MPI_Type_indexed(NumField, blocklengthsR, blockpossR, MPI_LONG_DOUBLE,& RblockType[iterCPU]);
        MPI_Type_commit(&RblockType[iterCPU]);
    }
    
    //process 1, 3, 5, ... - workers for r calculation (split theta)
    for (int iterCPU=1; iterCPU<numOfProcessT; ++iterCPU)
    {
        itertemp=0;
        for (int iterr=0; iterr<Nrp; ++iterr)
        {
            for (int iterf=0; iterf<NumField; ++iterf)
            {
                blocklengthsT[itertemp]=workerTl;
                blockpossT[itertemp]=bossTl+(iterCPU-1)*workerTl+iterr*Ntheta+iterf*NumPoints;
                ++itertemp;
            }
        }
        MPI_Type_indexed(matrixH, blocklengthsT, blockpossT, MPI_LONG_DOUBLE,& TblockType[iterCPU]);
        MPI_Type_commit(&TblockType[iterCPU]);
        MPI_Type_indexed(NumField, blocklengthsT, blockpossT, MPI_LONG_DOUBLE,& BoundaryType[iterCPU]);
        MPI_Type_commit(&BoundaryType[iterCPU]);
    }
    
    //Data distribution for BDF4
    BD4Type=new MPI_Datatype[numOfProcess];
    for (int iter=0; iter<NumField; ++iter)
    {
        blocklengthsR[iter]=bossP;
        blockpossR[iter]=iter*NumPoints;
    }
    MPI_Type_indexed(NumField, blocklengthsR, blockpossR, MPI_LONG_DOUBLE,& BD4Type[0]);
    MPI_Type_commit(&BD4Type[0]);
    for (int iterCPU=1; iterCPU<numOfProcess; ++iterCPU)
    {
        int offset=bossP+(iterCPU-1)*workerP;
        for (int iter=0; iter<NumField; ++iter)
        {
            blocklengthsR[iter]=workerP;
            blockpossR[iter]=offset+iter*NumPoints;
        }
        MPI_Type_indexed(NumField, blocklengthsR, blockpossR, MPI_LONG_DOUBLE,& BD4Type[iterCPU]);
        MPI_Type_commit(&BD4Type[iterCPU]);
    }
    
}

solver::~solver()
{
    
    for (int iterh=0; iterh<3; ++iterh)
    {
        delete HistoryFields[iterh];
    }
    
    if (cRank==0)
    {
        for (int iterh=0; iterh<NumField; ++iterh)
        {
            delete Hij[iterh];
        }
    }
    
    //local variables
    if (cRank%2==0)
    {
        fftwl_destroy_plan(fftr2c);
        fftwl_destroy_plan(ifftc2r);
        fftwl_destroy_plan(tempfftr2c);
        fftwl_destroy_plan(tempifftc2r);
    }
    else
    {
        fftwl_destroy_plan(dctr2r);
        fftwl_destroy_plan(tempdctr2r);
    }
    for (int iterh=0; iterh<NumField; ++iterh)
    {
        delete HijLocal[iterh];
    }
    
    for (int iterCPU=0; iterCPU<numOfProcessR; ++iterCPU)
    {
        MPI_Type_free(&RblockType[iterCPU]);
    }
    for (int iterCPU=0; iterCPU<numOfProcessT; ++iterCPU)
    {
        MPI_Type_free(&TblockType[iterCPU]);
        MPI_Type_free(&BoundaryType[iterCPU]);
    }
    for (int iterCPU=0; iterCPU<numOfProcess; ++iterCPU)
    {
        MPI_Type_free(&BD4Type[iterCPU]);
    }
    delete [] RblockType;
    delete [] TblockType;
    delete [] BoundaryType;
    delete [] BD4Type;
}
