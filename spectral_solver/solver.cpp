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
    MPI_Comm_rank(MPI_COMM_WORLD, &cRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcess);
    
    numOfProcessT=numOfProcess/2; //Num Of Processes should be even!
    numOfProcessR=numOfProcess-numOfProcessT;
    
    workerRl=ceil(Nrp/(double)numOfProcessR);
    workerTl=ceil(Ntheta/(double)numOfProcessT);
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
    
    workerBD4=floor(Nrp/(double)numOfProcess);
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
    
    odetempField2=gsl_matrix_alloc(jobBD4*NumField, Ntheta);
    odetempField3=gsl_matrix_alloc(jobBD4*NumField, Ntheta);
    iterFieldsLocal=gsl_matrix_alloc(jobBD4*NumField, Ntheta);
    
    HistoryFields.resize(3);
    for (int iterh=0; iterh<3; ++iterh)
    {
        HistoryFields[iterh]=gsl_matrix_alloc(jobBD4*NumField, Ntheta);
        gsl_matrix_set_zero(HistoryFields[iterh]);
    }
    
//    odetempField2=gsl_matrix_alloc(matrixH, Ntheta);
//    
//    HistoryFields.resize(3);
//    for (int iterh=0; iterh<3; ++iterh)
//    {
//        HistoryFields[iterh]=gsl_matrix_alloc(matrixH, Ntheta);
//        gsl_matrix_set_zero(HistoryFields[iterh]);
//    }
    
    if (cRank==0)
    {
        Fields=gsl_matrix_alloc(matrixH, Ntheta);
        gsl_matrix_set_zero(Fields);
        G=gsl_matrix_alloc(matrixH, Ntheta);
        gsl_matrix_set_zero(G);
        
        k1=gsl_matrix_alloc(matrixH, Ntheta);
        k2=gsl_matrix_alloc(matrixH, Ntheta);
        k3=gsl_matrix_alloc(matrixH, Ntheta);
        k4=gsl_matrix_alloc(matrixH, Ntheta);
        
        odetempField=gsl_matrix_alloc(matrixH, Ntheta);
        
        
        
        Hij.resize(NumField);
        for (int iterh=0; iterh<NumField; ++iterh)
        {
            Hij[iterh]=gsl_matrix_alloc(matrixH, Ntheta);
            gsl_matrix_set_zero(Hij[iterh]);
        }
    }
    
    //local variables
    if (cRank%2==0)
    {
        FieldsLocal=gsl_matrix_alloc(jobR, Ntheta);
        gsl_matrix_set_zero(FieldsLocal);
        tempFieldsLocal=gsl_matrix_alloc(jobR, Ntheta);
        dFieldsLocal=gsl_matrix_alloc(jobR, Ntheta);
        
        HijLocal.resize(NumField);
        for (int iterh=0; iterh<NumField; ++iterh)
        {
            HijLocal[iterh]=gsl_matrix_alloc(jobR, Ntheta);
            gsl_matrix_set_zero(HijLocal[iterh]);
        }
        
        for (int iter=0; iter<NumField; ++iter)
        {
            dFieldLocalView[iter]=gsl_matrix_submatrix(dFieldsLocal, iter*jobRl, 0, jobRl, Ntheta);
            for (int iterf=0; iterf<NumField; ++iterf)
            {
                HijLocalView[iterf*NumField+iter]=gsl_matrix_submatrix(HijLocal[iterf], iter*jobRl, 0, jobRl, Ntheta);
            }
        }
        
        GLocal=gsl_matrix_alloc(jobR, Ntheta);
        gsl_matrix_set_zero(GLocal);
        
        int n2[]={Ntheta};
        fftc=gsl_matrix_complex_alloc(jobR, Ntheta/2+1);
        fftr2c=fftw_plan_many_dft_r2c(1, n2, jobR, FieldsLocal->data, n2, 1, Ntheta, (fftw_complex *)fftc->data, n2, 1, Ntheta/2+1, FFTW_MEASURE);
        
        ifftc2r=fftw_plan_many_dft_c2r(1, n2, jobR, (fftw_complex *)fftc->data, n2, 1, Ntheta/2+1, dFieldsLocal->data, n2, 1, Ntheta, FFTW_MEASURE);
        
        tempfftr2c=fftw_plan_many_dft_r2c(1, n2, jobR, tempFieldsLocal->data, n2, 1, Ntheta, (fftw_complex *)fftc->data, n2, 1, Ntheta/2+1, FFTW_MEASURE);
        
        tempifftc2r=fftw_plan_many_dft_c2r(1, n2, jobR, (fftw_complex *)fftc->data, n2, 1, Ntheta/2+1, dFieldsLocal->data, n2, 1, Ntheta, FFTW_MEASURE);
    }
    else
    {
        HijLocal.resize(NumField);
        for (int iterh=0; iterh<NumField; ++iterh)
        {
            HijLocal[iterh]=gsl_matrix_alloc(Nrp, jobT);
            gsl_matrix_set_zero(HijLocal[iterh]);
        }
        
        tempFieldsLocal=gsl_matrix_alloc(Nrp, jobT);
        dctr=gsl_matrix_alloc(Nr, jobT);
        tempdctr=gsl_matrix_alloc(Nr, jobT);
        
        for (int iter=0; iter<NumField; ++iter)
        {
            dFieldLocalView[iter]=gsl_matrix_submatrix(dctr, 0, iter*jobTl, Nrp, jobTl);
            for (int iterf=0; iterf<NumField; ++iterf)
            {
                HijLocalView[iterf*NumField+iter]=gsl_matrix_submatrix(HijLocal[iterf], 0, iter*jobTl, Nrp, jobTl);
            }
        }
        
        int n1[]={Nr};
        fftw_r2r_kind kind[]={FFTW_REDFT00};
        dctr2r=fftw_plan_many_r2r(1, n1, jobT, dctr->data, n1, jobT, 1, dctr->data, n1, jobT, 1, kind, FFTW_MEASURE);
        tempdctr2r=fftw_plan_many_r2r(1, n1, jobT, tempdctr->data, n1, jobT, 1, dctr->data, n1, jobT, 1, kind, FFTW_MEASURE);
    }
    
    r=gsl_vector_alloc(Nrp);
    r2=gsl_vector_alloc(Nrp);
    for (int iter=0; iter<Nrp; ++iter)
    {
        r->data[iter]=cos(iter*PI/logicNr);
        r2->data[iter]=cos(iter*PI/logicNr)*cos(iter*PI/logicNr);
    }
    
    theta=gsl_vector_alloc(Ntheta);
    for (int iter=0; iter<Ntheta; ++iter)
    {
        theta->data[iter]=2*PI*iter/Ntheta;
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
    MPI_Type_indexed(NumField, blocklengthsR, blockpossR, MPI_DOUBLE, &RblockType[0]);
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
    MPI_Type_indexed(matrixH, blocklengthsT, blockpossT, MPI_DOUBLE, &TblockType[0]);
    MPI_Type_commit(&TblockType[0]);
    MPI_Type_indexed(NumField, blocklengthsT, blockpossT, MPI_DOUBLE, &BoundaryType[0]);
    MPI_Type_commit(&BoundaryType[0]);
    
    //process 2, 4, 6, ... - workers for theta calculation (split r)
    for (int iterCPU=1; iterCPU<numOfProcessR; ++iterCPU)
    {
        for (int iter=0; iter<NumField; ++iter)
        {
            blocklengthsR[iter]=workerPointsR;
            blockpossR[iter]=bossPointsR+(iterCPU-1)*workerPointsR+iter*NumPoints;
        }
        MPI_Type_indexed(NumField, blocklengthsR, blockpossR, MPI_DOUBLE, &RblockType[iterCPU]);
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
        MPI_Type_indexed(matrixH, blocklengthsT, blockpossT, MPI_DOUBLE, &TblockType[iterCPU]);
        MPI_Type_commit(&TblockType[iterCPU]);
        MPI_Type_indexed(NumField, blocklengthsT, blockpossT, MPI_DOUBLE, &BoundaryType[iterCPU]);
        MPI_Type_commit(&BoundaryType[iterCPU]);
    }
    
    //Data distribution for BDF4
    BD4Type=new MPI_Datatype[numOfProcess];
    for (int iter=0; iter<NumField; ++iter)
    {
        blocklengthsR[iter]=bossP;
        blockpossR[iter]=iter*NumPoints;
    }
    MPI_Type_indexed(NumField, blocklengthsR, blockpossR, MPI_DOUBLE, &BD4Type[0]);
    MPI_Type_commit(&BD4Type[0]);
    for (int iterCPU=1; iterCPU<numOfProcess; ++iterCPU)
    {
        int offset=bossP+(iterCPU-1)*workerP;
        for (int iter=0; iter<NumField; ++iter)
        {
            blocklengthsR[iter]=workerP;
            blockpossR[iter]=offset+iter*NumPoints;
        }
        MPI_Type_indexed(NumField, blocklengthsR, blockpossR, MPI_DOUBLE, &BD4Type[iterCPU]);
        MPI_Type_commit(&BD4Type[iterCPU]);
    }
    
}

solver::~solver()
{
    
    gsl_matrix_free(odetempField2);
    gsl_matrix_free(odetempField3);
    gsl_matrix_free(iterFieldsLocal);
    
    for (int iterh=0; iterh<3; ++iterh)
    {
        gsl_matrix_free(HistoryFields[iterh]);
    }
    
    if (cRank==0)
    {
        gsl_matrix_free(Fields);
        gsl_matrix_free(G);
        
        gsl_matrix_free(k1);
        gsl_matrix_free(k2);
        gsl_matrix_free(k3);
        gsl_matrix_free(k4);
        
        gsl_matrix_free(odetempField);
        
        
        for (int iterh=0; iterh<NumField; ++iterh)
        {
            gsl_matrix_free(Hij[iterh]);
        }
    }
    
    //local variables
    if (cRank%2==0)
    {
        gsl_matrix_free(GLocal);
        gsl_matrix_free(dFieldsLocal);
        gsl_matrix_free(tempFieldsLocal);
        gsl_matrix_complex_free(fftc);
        fftw_destroy_plan(fftr2c);
        fftw_destroy_plan(ifftc2r);
        fftw_destroy_plan(tempfftr2c);
        fftw_destroy_plan(tempifftc2r);
    }
    else
    {
        gsl_matrix_free(dctr);
        gsl_matrix_free(tempdctr);
        fftw_destroy_plan(dctr2r);
        fftw_destroy_plan(tempdctr2r);
    }
    for (int iterh=0; iterh<NumField; ++iterh)
    {
        gsl_matrix_free(HijLocal[iterh]);
    }
    
    gsl_vector_free(r);
    gsl_vector_free(r2);
    gsl_vector_free(theta);
    
    gsl_matrix_free(FieldsLocal);
    
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
    
    MPI_Finalize();
}
