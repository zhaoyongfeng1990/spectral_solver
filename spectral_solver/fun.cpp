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
    for (int iter=0; iter<jobPointsRl; ++iter)
    {
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            f[iterf]=FieldsLocal->data[iterf*jobPointsRl+iter];
        }
        
        //Calculation for H functions. Use f[] to express.
        //Term in front of dField 1
        double powfm1=pow(f[1]/Kh,M-1);
        double powfm=powfm1*f[1]/Kh;
        double powfmp1=powfm+1;
        HijLocal[0]->data[iter]=(Drho+Drho0*powfm)/powfmp1;
        HijLocal[0]->data[iter+jobPointsRl]=0;
        HijLocal[0]->data[iter+2*jobPointsRl]=0;
        //Term in front of dField 2
        HijLocal[1]->data[iter]=f[0]*(Drho0-Drho)*powfm1*M/Kh/powfmp1/powfmp1;
        HijLocal[1]->data[iter+jobPointsRl]=Dh;
        HijLocal[1]->data[iter+2*jobPointsRl]=0;
        //Term in front of dField 3
        HijLocal[2]->data[iter]=0;
        HijLocal[2]->data[iter+jobPointsRl]=0;
        HijLocal[2]->data[iter+2*jobPointsRl]=Dn;
        
        double f2=f[2]*f[2];
        //Ending
        //Calculation for G function.
        GLocal->data[iter]=Gamma*f2*f[0]/(f2+Kn2);
        GLocal->data[iter+jobPointsRl]=Alpha*f[0]-Beta*f[1];
        GLocal->data[iter+2*jobPointsRl]=-GLocal->data[iter];
        
//                HijLocal[0]->data[iter]=1;
//                HijLocal[0]->data[iter+jobPointsRl]=0;
//                HijLocal[0]->data[iter+2*jobPointsRl]=0;
//                //Term in front of dField 2
//                HijLocal[1]->data[iter]=0;
//                HijLocal[1]->data[iter+jobPointsRl]=1;
//                HijLocal[1]->data[iter+2*jobPointsRl]=0;
//                //Term in front of dField 3
//                HijLocal[2]->data[iter]=0;
//                HijLocal[2]->data[iter+jobPointsRl]=0;
//                HijLocal[2]->data[iter+2*jobPointsRl]=1;
//        
//                //Ending
//                //Calculation for G function.
//                GLocal->data[iter]=0;
//                GLocal->data[iter+jobPointsRl]=0;
//                GLocal->data[iter+2*jobPointsRl]=0;
        //....
    }
}


void solver::Fun(gsl_matrix *result)
{
    if (0==cRank)
    {
        for (int iterCPU=1; iterCPU<numOfProcessR; ++iterCPU)
        {
            MPI_Send(Fields->data, 1, RblockType[iterCPU], 2*iterCPU, 2*iterCPU, MPI_COMM_WORLD);
        }
        for (int iterCPU=1; iterCPU<numOfProcessT; ++iterCPU)
        {
            MPI_Send(Fields->data, 1, TblockType[iterCPU], 2*iterCPU+1, 2*iterCPU+1, MPI_COMM_WORLD);
        }
        MPI_Send(Fields->data, 1, TblockType[0], 1, 1, MPI_COMM_WORLD);
        //MPI_Send(Fields->data, 1, RblockType[0], 0, 0, MPI_COMM_WORLD);
//        for (int iterf=0; iterf<NumField; ++iterf)
//        {
//            for (int iter=0; iter<jobPointsRl; ++iter)
//            {
//                FieldsLocal->data[iter+iterf*jobPointsRl]=Fields->data[iter+iterf*NumPoints];
//            }
        //        }
        MPI_Isend(Fields->data, 1, RblockType[0], 0, 0, MPI_COMM_WORLD, &request);
    }
    
    if (cRank%2==0)
    {
        if (cRank!=0)
            MPI_Recv(FieldsLocal->data, jobPointsR, MPI_DOUBLE, 0, cRank, MPI_COMM_WORLD, &status);
        else
            MPI_Irecv(FieldsLocal->data, jobPointsR, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
        //H functions and G term
        HGFuns();
        if (0!=cRank)
        {
            for (int iterf=0; iterf<NumField; ++iterf)
            {
                MPI_Send(HijLocal[iterf]->data, jobPointsR, MPI_DOUBLE, 0, 100+10*iterf+cRank, MPI_COMM_WORLD);
            }
        }
        else
        {
            for (int iterf=0; iterf<NumField; ++iterf)
            {
                MPI_Isend(HijLocal[iterf]->data, jobPointsR, MPI_DOUBLE, 0, 100+10*iterf, MPI_COMM_WORLD, &request);
            }
        }
        if (0==cRank)
        {
            for (int iterf=0; iterf<NumField; ++iterf)
            {
//                for (int iterff=0; iterff<NumField; ++iterff)
//                {
//                    for (int iter=0; iter<jobPointsRl; ++iter)
//                    {
//                        Hij[iterf]->data[iter+iterff*jobPointsRl]=HijLocal[iterf]->data[iter+iterff*NumPoints];
//                    }
//                }
                MPI_Irecv(Hij[iterf]->data, 1, RblockType[0], 0, 100+10*iterf, MPI_COMM_WORLD, &request);
                for (int iterCPU=1; iterCPU<numOfProcessR; ++iterCPU)
                {
                    MPI_Recv(Hij[iterf]->data, 1, RblockType[iterCPU], iterCPU*2, 100+10*iterf+iterCPU*2, MPI_COMM_WORLD, &status);
                }
                for (int iterCPU=0; iterCPU<numOfProcessT; ++iterCPU)
                {
                    MPI_Send(Hij[iterf]->data, 1, TblockType[iterCPU], iterCPU*2+1, 100+10*iterf+iterCPU*2+1, MPI_COMM_WORLD);
                }
            }
        }
        
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
            MPI_Send(GLocal->data, jobPointsR, MPI_DOUBLE, 0, 200+cRank, MPI_COMM_WORLD);
        }
        else
        {
//            for (int iterff=0; iterff<NumField; ++iterff)
//            {
//                for (int iter=0; iter<jobPointsRl; ++iter)
//                {
//                    G->data[iter+iterff*jobPointsRl]=GLocal->data[iter+iterff*NumPoints];
//                }
//            }
            MPI_Isend(GLocal->data, jobPointsR, MPI_DOUBLE, 0, 200, MPI_COMM_WORLD, &request);
        }
    }
    else
    {
        MPI_Recv(dctr->data, jobPointsT, MPI_DOUBLE, 0, cRank, MPI_COMM_WORLD, &status);
        //derivative of r term
        dr(1);
        
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            MPI_Recv(HijLocal[iterf]->data, jobPointsT, MPI_DOUBLE, 0, 100+10*iterf+cRank, MPI_COMM_WORLD, &status);
        }
        
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
        MPI_Send(dctr->data, jobPointsT, MPI_DOUBLE, 0, 200+cRank, MPI_COMM_WORLD);
        
    }
    
    if (0==cRank)
    {
        MPI_Irecv(G->data, 1, RblockType[0], 0, 200, MPI_COMM_WORLD, &request);
        for (int iterCPU=1; iterCPU<numOfProcessR; ++iterCPU)
        {
            MPI_Recv(G->data, 1, RblockType[iterCPU], iterCPU*2, 200+iterCPU*2, MPI_COMM_WORLD, &status);
        }
        for (int iterCPU=0; iterCPU<numOfProcessT; ++iterCPU)
        {
            MPI_Recv(result->data, 1, TblockType[iterCPU], iterCPU*2+1, 200+iterCPU*2+1, MPI_COMM_WORLD, &status);
        }
        gsl_matrix_add(result, G);
    }
}
