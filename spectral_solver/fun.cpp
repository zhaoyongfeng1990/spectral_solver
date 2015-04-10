//
//  fun.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/9.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
using std::pow;

void solver::HGFuns()
{
    long double f[NumField];
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iter=0; iter<NumPoints; ++iter)
    {
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            f[iterf]=Fields[iterf*NumPoints+iter];
        }
        
        //Calculation for H functions. Use f[] to express.
#ifdef FU_MODEL
        //Term in front of dField 1
        long double powfm1=pow(f[1]/Kh,M-1);
        long double powfm=powfm1*f[1]/Kh;
        long double powfmp1=powfm+1;
        Hij[0]->data[iter]=(Drho+Drho0*powfm)/powfmp1;
        Hij[0]->data[iter+NumPoints]=0;
        Hij[0]->data[iter+2*NumPoints]=0;
        //Term in front of dField 2
        Hij[1]->data[iter]=f[0]*(Drho0-Drho)*powfm1*M/Kh/powfmp1/powfmp1;
        Hij[1]->data[iter+NumPoints]=Dh;
        Hij[1]->data[iter+2*NumPoints]=0;
        //Term in front of dField 3
        Hij[2]->data[iter]=0;
        Hij[2]->data[iter+NumPoints]=0;
        Hij[2]->data[iter+2*NumPoints]=Dn;
        
        long double f2=f[2]*f[2];
        //Ending
        //Calculation for G function.
        G[iter]=Gamma*f2*f[0]/(f2+Kn2);
        G[iter+NumPoints]=Alpha*f[0]-Beta*f[1];
        G[iter+2*NumPoints]=-G[iter];
#endif
#ifdef LINEAR_TEST_MODEL
        Hij[0]->data[iter]=1;
        Hij[0]->data[iter+NumPoints]=0;
        Hij[0]->data[iter+2*NumPoints]=0;
        //Term in front of dField 2
        Hij[1]->data[iter]=0;
        Hij[1]->data[iter+NumPoints]=1;
        Hij[1]->data[iter+2*NumPoints]=0;
        //Term in front of dField 3
        Hij[2]->data[iter]=0;
        Hij[2]->data[iter+NumPoints]=0;
        Hij[2]->data[iter+2*NumPoints]=1;
        
        //Ending
        //Calculation for G function.
        G[iter]=0;
        G[iter+NumPoints]=0;
        G[iter+2*NumPoints]=0;
#endif
        //....
    }
}


void solver::Fun(matrix<long double>& result)
{
    //H functions and G term
    HGFuns();
    
    //derivative of r term
    dr(1);
    
#ifdef PUNISHTERM
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iter=0; iter<NumField; ++iter)
    {
        for (int itert=0; itert<Ntheta; ++itert)
        {
            boundary.ele(iter, itert)=dFields.ele(iter*Nrp, itert)*punish;
        }
    }
#endif
    
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iter=0; iter<totalPoints; ++iter)
    {
        tempFields[iter]=0;
    }
    //matrix<long double>_set_zero(tempFields);
    
    for (int iterdf=0; iterdf<NumField; ++iterdf)
    {
        caltempFields=*Hij[iterdf];
        for (int iterp=0; iterp<NumPoints; ++iterp)
        {
            for (int iter=0; iter<NumField; ++iter)
            {
                caltempFields.data[iter*NumPoints+iterp]*=dFields[iterp+iterdf*NumPoints];
            }
        }
        tempFields+=caltempFields;
    }
    
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            for (int iterr=0; iterr<Nrp; ++iterr)
            {
                for (int iter=0; iter<Ntheta; ++iter)
                {
                tempFields.ele(iterr+iterf*Nrp, iter)*=r[iterr];
            }
        }
    }
    
    dr(0);
#ifdef MULTIPROCESS
#pragma omp parallel
    {
#pragma omp for
#endif
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            for (int iterr=0; iterr<Nrp; ++iterr)
            {
                for (int iter=0; iter<Ntheta; ++iter)
                {
                    dFields.ele(iterr+iterf*Nrp, iter)/=(r[iterr]*rp2);
                }
            }
        }
#ifdef MULTIPROCESS
#pragma omp for
#endif
    for (int iter=0; iter<totalPoints; ++iter)
    {
        G[iter]+=dFields[iter];
    }
#ifdef MULTIPROCESS
}
#endif
    //matrix<long double>_add(G, dFields);
    
    //derivative of theta term
    dtheta(1);
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iter=0; iter<totalPoints; ++iter)
    {
        tempFields[iter]=0;
    }
    //matrix<long double>_set_zero(tempFields);
    
    for (int iterdf=0; iterdf<NumField; ++iterdf)
    {
        //caltempFields=*Hij[iterdf];
        for (int iterp=0; iterp<NumPoints; ++iterp)
        {
            for (int iter=0; iter<NumField; ++iter)
            {
                //caltempFields.data[iter*NumPoints+iterp]*=dFields[iterp+iterdf*NumPoints];
                Hij[iterdf]->data[iter*NumPoints+iterp]*=dFields[iterp+iterdf*NumPoints];
            }
        }
        //tempFields+=caltempFields;
        tempFields+=*Hij[iterdf];
    }
    dtheta(0);
    //matrix<long double>_scale(dFields, rp2); //rp2=1.0/r/r
#ifdef MULTIPROCESS
#pragma omp parallel
    {
#pragma omp for
#endif
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            for (int iterr=0; iterr<Nrp; ++iterr)
            {
                for (int iter=0; iter<Ntheta; ++iter)
                {
                    dFields.ele(iterr+iterf*Nrp, iter)/=(r2[iterr]*rp2);
                }
            }
        }
#ifdef MULTIPROCESS
#pragma omp for
#endif
        for (int iter=0; iter<totalPoints; ++iter)
        {
            G[iter]+=dFields[iter];
        }
#ifdef MULTIPROCESS
    }
#endif
    
    //matrix<long double>_add(G, dFields);
    
#ifdef PUNISHTERM
    //#ifdef MULTIPROCESS
    //#pragma omp parallel for
    //#endif
    for (int iter=0; iter<NumField; ++iter)
    {
        for (int itert=0; itert<Ntheta; ++iter)
        {
            G.ele(iter*Nrp, itert)-=boundary.ele(iter, itert);
        }
    }
#endif
    result=G;
}
