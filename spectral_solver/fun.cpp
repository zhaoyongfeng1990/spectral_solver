//
//  fun.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/9.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
#include <iostream>
using namespace std;
using std::pow;

void solver::HGFuns()
{
    double f[NumField];
    int idx[NumField];
    for (int iter=0; iter<jobPointsRl; ++iter)
    {
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            idx[iterf]=iterf*jobPointsRl+iter;
            f[iterf]=FieldsLocal->data[idx[iterf]];
        }
        
        //Calculation for H functions. Use f[] to express.
#ifdef FU_MODEL
        //Term in front of dField 1
        double powfm1=pow(f[1]/Kh,M-1);
        double powfm=powfm1*f[1]/Kh;
        double powfmp1=powfm+1;
        HijLocal[0]->data[idx[0]]=(Drho+Drho0*powfm)/powfmp1;
        HijLocal[0]->data[idx[1]]=0;
        HijLocal[0]->data[idx[2]]=0;
        //Term in front of dField 2
        HijLocal[1]->data[idx[0]]=f[0]*(Drho0-Drho)*powfm1*M/Kh/powfmp1/powfmp1;
        HijLocal[1]->data[idx[1]]=Dh;
        HijLocal[1]->data[idx[2]]=0;
        //Term in front of dField 3
        HijLocal[2]->data[idx[0]]=0;
        HijLocal[2]->data[idx[1]]=0;
        HijLocal[2]->data[idx[2]]=Dn;
        
        double f2=f[2]*f[2];
        //Ending
        //Calculation for G function.
        GLocal->data[idx[0]]=Gamma*f2*f[0]/(f2+Kn2);
        GLocal->data[idx[1]]=Alpha*f[0]-Beta*f[1];
        GLocal->data[idx[2]]=-GLocal->data[iter];
#endif
#ifdef LINEAR_TEST_MODEL
        HijLocal[0]->data[idx[0]]=1;
        HijLocal[0]->data[idx[1]]=0;
        HijLocal[0]->data[idx[2]]=0;
        //Term in front of dField 2
        HijLocal[1]->data[idx[0]]=0;
        HijLocal[1]->data[idx[1]]=1;
        HijLocal[1]->data[idx[2]]=0;
        //Term in front of dField 3
        HijLocal[2]->data[idx[0]]=0;
        HijLocal[2]->data[idx[1]]=0;
        HijLocal[2]->data[idx[2]]=1;
        
        //Ending
        //Calculation for G function.
        GLocal->data[idx[0]]=0;
        GLocal->data[idx[1]]=0;
        GLocal->data[idx[2]]=0;
        //....
#endif
#ifdef INTERACTION_MODIFIED_FU
        //Term in front of dField 1
        double h1kh1=f[2]/Kh1;
        double h2kh2=f[3]/Kh2;
        double powh1=pow(h1kh1,M-1);
        double powh2=pow(h2kh2,M-1);
        double powh1m=powh1*h1kh1;
        double powh2m=powh2*h2kh2;
        double powh1p1=powh1m+1;
        double powh2p1=powh2m+1;
        double h1KD2=f[2]/KD2;
        double h2KC2=f[3]/KC2;
        double h1KD2p1=h1KD2+1;
        double h2KC2p1=h2KC2+1;
        double t1=Drho1+KC1*h2KC2/h2KC2p1;
        double t2=Drho2+KD1*h1KD2/h1KD2p1;
        HijLocal[0]->data[idx[0]]=Drho+t1/powh1p1;
        HijLocal[0]->data[idx[1]]=0;
        HijLocal[0]->data[idx[2]]=0;
        HijLocal[0]->data[idx[3]]=0;
        HijLocal[0]->data[idx[4]]=0;
        
        //Term in front of dField 2
        HijLocal[1]->data[idx[0]]=0;
        HijLocal[1]->data[idx[1]]=Drho+t2/powh2p1;
        HijLocal[1]->data[idx[2]]=0;
        HijLocal[1]->data[idx[3]]=0;
        HijLocal[1]->data[idx[4]]=0;
        
        //Term in front of dField 3
        HijLocal[2]->data[idx[0]]=-t1*powh1*f[0]*M/Kh1/powh1p1/powh1p1;
        HijLocal[2]->data[idx[1]]=KD1/KD2*f[1]/h1KD2p1/h1KD2p1/powh2p1;
        HijLocal[2]->data[idx[2]]=Dh1;
        HijLocal[2]->data[idx[3]]=0;
        HijLocal[2]->data[idx[4]]=0;
        
        //Term in front of dField 4
        HijLocal[3]->data[idx[0]]=KC1/KC2*f[0]/h2KC2p1/h2KC2p1/powh1p1;
        HijLocal[3]->data[idx[1]]=-t2*powh2*f[1]*M/Kh2/powh2p1/powh2p1;
        HijLocal[3]->data[idx[2]]=0;
        HijLocal[3]->data[idx[3]]=Dh2;
        HijLocal[3]->data[idx[4]]=0;
        
        //Term in front of dField 5
        HijLocal[4]->data[idx[0]]=0;
        HijLocal[4]->data[idx[1]]=0;
        HijLocal[4]->data[idx[2]]=0;
        HijLocal[4]->data[idx[3]]=0;
        HijLocal[4]->data[idx[4]]=Dn;
        
        double f2=f[4]*f[4];
        //Ending
        //Calculation for G function.
        GLocal->data[idx[0]]=Gamma1*f2*f[0]/(f2+Kn1);
        GLocal->data[idx[1]]=Gamma2*f2*f[1]/(f2+Kn2);
        GLocal->data[idx[2]]=Alpha1*f[0]-Beta1*f[2];
        GLocal->data[idx[3]]=Alpha2*f[1]-Beta2*f[3];
        GLocal->data[idx[4]]=-GLocal->data[idx[0]]-GLocal->data[idx[1]];
#endif
    }
}

void solver::HFunsForR()    //should keep as identical with HGFuns().
{
    double f[NumField];
    int idx[NumField];
    for (int iter=0; iter<Nrp; ++iter)
    {
        for (int itert=0; itert<jobTl; ++itert)
        {
            for (int iterff=0; iterff<NumField; ++iterff)
            {
                idx[iterff]=iter*jobT+itert+iterff*jobTl;
                f[iterff]=dctr->data[idx[iterff]];
            }
            
#ifdef FU_MODEL
            //Calculation for H functions. Use f[] to express.
            //Term in front of dField 1
            double powfm1=pow(f[1]/Kh,M-1);
            double powfm=powfm1*f[1]/Kh;
            double powfmp1=powfm+1;
            HijLocal[0]->data[idx[0]]=(Drho+Drho0*powfm)/powfmp1;
            HijLocal[0]->data[idx[1]]=0;
            HijLocal[0]->data[idx[2]]=0;
            //Term in front of dField 2
            HijLocal[1]->data[idx[0]]=f[0]*(Drho0-Drho)*powfm1*M/Kh/powfmp1/powfmp1;
            HijLocal[1]->data[idx[1]]=Dh;
            HijLocal[1]->data[idx[2]]=0;
            //Term in front of dField 3
            HijLocal[2]->data[idx[0]]=0;
            HijLocal[2]->data[idx[1]]=0;
            HijLocal[2]->data[idx[2]]=Dn;
#endif
#ifdef LINEAR_TEST_MODEL
            
            HijLocal[0]->data[idx[0]]=1;
            HijLocal[0]->data[idx[1]]=0;
            HijLocal[0]->data[idx[2]]=0;
            //Term in front of dField 2
            HijLocal[1]->data[idx[0]]=0;
            HijLocal[1]->data[idx[1]]=1;
            HijLocal[1]->data[idx[2]]=0;
            //Term in front of dField 3
            HijLocal[2]->data[idx[0]]=0;
            HijLocal[2]->data[idx[1]]=0;
            HijLocal[2]->data[idx[2]]=1;
            
            //Ending
            //....
#endif
#ifdef INTERACTION_MODIFIED_FU
            //Term in front of dField 1
            double h1kh1=f[2]/Kh1;
            double h2kh2=f[3]/Kh2;
            double powh1=pow(h1kh1,M-1);
            double powh2=pow(h2kh2,M-1);
            double powh1m=powh1*h1kh1;
            double powh2m=powh2*h2kh2;
            double powh1p1=powh1m+1;
            double powh2p1=powh2m+1;
            double h1KD2=f[2]/KD2;
            double h2KC2=f[3]/KC2;
            double h1KD2p1=h1KD2+1;
            double h2KC2p1=h2KC2+1;
            double t1=Drho1+KC1*h2KC2/h2KC2p1;
            double t2=Drho2+KD1*h1KD2/h1KD2p1;
            HijLocal[0]->data[idx[0]]=Drho+t1/powh1p1;
            HijLocal[0]->data[idx[1]]=0;
            HijLocal[0]->data[idx[2]]=0;
            HijLocal[0]->data[idx[3]]=0;
            HijLocal[0]->data[idx[4]]=0;
            
            //Term in front of dField 2
            HijLocal[1]->data[idx[0]]=0;
            HijLocal[1]->data[idx[1]]=Drho+t2/powh2p1;
            HijLocal[1]->data[idx[2]]=0;
            HijLocal[1]->data[idx[3]]=0;
            HijLocal[1]->data[idx[4]]=0;
            
            //Term in front of dField 3
            HijLocal[2]->data[idx[0]]=-t1*powh1*f[0]*M/Kh1/powh1p1/powh1p1;
            HijLocal[2]->data[idx[1]]=KD1/KD2*f[1]/h1KD2p1/h1KD2p1/powh2p1;
            HijLocal[2]->data[idx[2]]=Dh1;
            HijLocal[2]->data[idx[3]]=0;
            HijLocal[2]->data[idx[4]]=0;
            
            //Term in front of dField 4
            HijLocal[3]->data[idx[0]]=KC1/KC2*f[0]/h2KC2p1/h2KC2p1/powh1p1;
            HijLocal[3]->data[idx[1]]=-t2*powh2*f[1]*M/Kh2/powh2p1/powh2p1;
            HijLocal[3]->data[idx[2]]=0;
            HijLocal[3]->data[idx[3]]=Dh2;
            HijLocal[3]->data[idx[4]]=0;
            
            //Term in front of dField 5
            HijLocal[4]->data[idx[0]]=0;
            HijLocal[4]->data[idx[1]]=0;
            HijLocal[4]->data[idx[2]]=0;
            HijLocal[4]->data[idx[3]]=0;
            HijLocal[4]->data[idx[4]]=Dn;
#endif

        }
    }
}


void solver::Fun(gsl_matrix *result)
{
    if (0==cRank)
    {
        for (int iterCPU=1; iterCPU<numOfProcessT; ++iterCPU)
        {
            MPI_Ssend(Fields->data, 1, TblockType[iterCPU], 2*iterCPU+1, 2*iterCPU+1, MPI_COMM_WORLD);
        }
        MPI_Ssend(Fields->data, 1, TblockType[0], 1, 1, MPI_COMM_WORLD);
        for (int iterCPU=1; iterCPU<numOfProcessR; ++iterCPU)
        {
            MPI_Ssend(Fields->data, 1, RblockType[iterCPU], 2*iterCPU, 2*iterCPU, MPI_COMM_WORLD);
        }
        //MPI_Sendrecv(Fields->data, 1, RblockType[0], 0, 0, FieldsLocal->data, jobPointsR, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            for (int iter=0; iter<jobPointsRl; ++iter)
            {
                FieldsLocal->data[iter+iterf*jobPointsRl]=Fields->data[iter+iterf*NumPoints];
            }
        }
    }
    
    if (cRank%2==0)
    {
        if (cRank!=0)
            MPI_Recv(FieldsLocal->data, jobPointsR, MPI_DOUBLE, 0, cRank, MPI_COMM_WORLD, &status);
        //H functions and G term
        HGFuns();
        
        //derivative of theta term
        dtheta(1);
        gsl_matrix_set_zero(tempFieldsLocal);
        
        for (int iterdf=0; iterdf<NumField; ++iterdf)
        {
            for (int iter=0; iter<NumField; ++iter)
            {
                gsl_matrix_mul_elements(&HijLocalView[iterdf*NumField+iter].matrix, &dFieldLocalView[iterdf].matrix);
            }
            gsl_matrix_add(tempFieldsLocal, HijLocal[iterdf]);
        }
        dtheta(0);
        gsl_matrix_scale(dFieldsLocal, rp2); //rp2=1.0/R/R
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            for (int iter=0; iter<jobRl; ++iter)
            {
                gsl_vector_view temp=gsl_matrix_row(dFieldsLocal, iterf*jobRl+iter);
                if (0==cRank)
                {
                    gsl_vector_scale(&temp.vector, r2->data[iter]);
                }
                else
                {
                    gsl_vector_scale(&temp.vector, r2->data[iter+bossRl+(cRank/2-1)*workerRl]);
                }
            }
        }
        gsl_matrix_add(GLocal, dFieldsLocal);
        if (cRank!=0)
        {
            MPI_Ssend(GLocal->data, jobPointsR, MPI_DOUBLE, 0, 100+cRank, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(dctr->data, jobPointsT, MPI_DOUBLE, 0, cRank, MPI_COMM_WORLD, &status);
        //derivative of r term
        HFunsForR();
        dr(1);
        
#ifdef PUNISHTERM
        gsl_vector_view tempboundary=gsl_matrix_row(dctr, 0);
        gsl_vector_memcpy(boundary, &tempboundary.vector);
        gsl_vector_scale(boundary, punish);
#endif
        
        for (int iter=0; iter<jobPointsT; ++iter)
        {
            tempdctr->data[iter]=0;
        }
        for (int iterdf=0; iterdf<NumField; ++iterdf)
        {
            for (int iter=0; iter<NumField; ++iter)
            {
                gsl_matrix_mul_elements(&HijLocalView[iterdf*NumField+iter].matrix, &dFieldLocalView[iterdf].matrix);
            }
            for (int iter=0; iter<jobPointsT; ++iter)
            {
                tempdctr->data[iter]+=HijLocal[iterdf]->data[iter];
            }
        }
        
        for (int iter=0; iter<jobT; ++iter)
        {
            gsl_vector_view temp=gsl_matrix_subcolumn(tempdctr, iter, 0, Nrp);
            gsl_vector_mul(&temp.vector, r);
        }
        
        dr(0);
        
        for (int iter=0; iter<jobPointsT; ++iter)
        {
            dctr->data[iter]*=rp2; //rp2=1.0/R/R
        }
        
        for (int iter=0; iter<jobT; ++iter)
        {
            gsl_vector_view temp=gsl_matrix_subcolumn(dctr, iter, 0, Nrp);
            gsl_vector_div(&temp.vector, r);
        }
        
#ifdef PUNISHTERM
        //tempboundary=gsl_matrix_row(dctr, 0);
        gsl_vector_sub(&tempboundary.vector, boundary);
#endif
        
        MPI_Ssend(dctr->data, jobPointsT, MPI_DOUBLE, 0, 100+cRank, MPI_COMM_WORLD);
        
    }
    
    if (0==cRank)
    {
        //MPI_Sendrecv(GLocal->data, jobPointsR, MPI_DOUBLE, 0, 200, G->data, 1, RblockType[0], 0, 200, MPI_COMM_WORLD, &status);
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            for (int iter=0; iter<jobPointsRl; ++iter)
            {
                G->data[iter+iterf*NumPoints]=GLocal->data[iter+iterf*jobPointsRl];
            }
        }
        for (int iterCPU=1; iterCPU<numOfProcessR; ++iterCPU)
        {
            MPI_Recv(G->data, 1, RblockType[iterCPU], iterCPU*2, 100+iterCPU*2, MPI_COMM_WORLD, &status);
        }
        for (int iterCPU=0; iterCPU<numOfProcessT; ++iterCPU)
        {
            MPI_Recv(result->data, 1, TblockType[iterCPU], iterCPU*2+1, iterCPU*2+101, MPI_COMM_WORLD, &status);
        }
        gsl_matrix_add(result, G);
    }
}