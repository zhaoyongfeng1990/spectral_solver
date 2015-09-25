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
    int idx[NumField];
#ifdef MULTIPROCESS
#pragma omp parallel for
#endif
    for (int iter=0; iter<NumPoints; ++iter)
    {
        for (int iterf=0; iterf<NumField; ++iterf)
        {
            idx[iterf]=iterf*NumPoints+iter;
            f[iterf]=Fields[idx[iterf]];
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
#ifdef INTERACTION_MODIFIED_FU
        //Term in front of dField 1
        long double h1kh1=f[2]/Kh1;
        long double h2kh2=f[3]/Kh2;
        long double powh1=pow(h1kh1,M-1);
        long double powh2=pow(h2kh2,M-1);
        long double powh1m=powh1*h1kh1;
        long double powh2m=powh2*h2kh2;
        long double powh1p1=powh1m+1;
        long double powh2p1=powh2m+1;
        long double h1KD2=f[2]/KD2;
        long double h2KC2=f[3]/KC2;
        long double h1KD2p1=h1KD2+1;
        long double h2KC2p1=h2KC2+1;
        long double t1=Drho1+KC1*h2KC2/h2KC2p1;
        long double t2=Drho2+KD1*h1KD2/h1KD2p1;
        Hij[0]->data[idx[0]]=Drho+t1/powh1p1;
        Hij[0]->data[idx[1]]=0;
        Hij[0]->data[idx[2]]=0;
        Hij[0]->data[idx[3]]=0;
        Hij[0]->data[idx[4]]=0;
        
        //Term in front of dField 2
        Hij[1]->data[idx[0]]=0;
        Hij[1]->data[idx[1]]=Drho+t2/powh2p1;
        Hij[1]->data[idx[2]]=0;
        Hij[1]->data[idx[3]]=0;
        Hij[1]->data[idx[4]]=0;
        
        //Term in front of dField 3
        Hij[2]->data[idx[0]]=-t1*powh1*f[0]*M/Kh1/powh1p1/powh1p1;
        Hij[2]->data[idx[1]]=KD1/KD2*f[1]/h1KD2p1/h1KD2p1/powh2p1;
        Hij[2]->data[idx[2]]=Dh1;
        Hij[2]->data[idx[3]]=0;
        Hij[2]->data[idx[4]]=0;
        
        //Term in front of dField 4
        Hij[3]->data[idx[0]]=KC1/KC2*f[0]/h2KC2p1/h2KC2p1/powh1p1;
        Hij[3]->data[idx[1]]=-t2*powh2*f[1]*M/Kh2/powh2p1/powh2p1;
        Hij[3]->data[idx[2]]=0;
        Hij[3]->data[idx[3]]=Dh2;
        Hij[3]->data[idx[4]]=0;
        
        //Term in front of dField 5
        Hij[4]->data[idx[0]]=0;
        Hij[4]->data[idx[1]]=0;
        Hij[4]->data[idx[2]]=0;
        Hij[4]->data[idx[3]]=0;
        Hij[4]->data[idx[4]]=Dn;
        
        double f2=f[4]*f[4];
        //Ending
        //Calculation for G function.
        G.data[idx[0]]=Gamma1*f2*f[0]/(f2+Kn1);
        G.data[idx[1]]=Gamma2*f2*f[1]/(f2+Kn2);
        G.data[idx[2]]=Alpha1*f[0]-Beta1*f[2];
        G.data[idx[3]]=Alpha2*f[1]-Beta2*f[3];
        G.data[idx[4]]=-G.data[idx[0]]-G.data[idx[1]];
#endif
#ifdef MASA_CROSSTALK
        //Term in front of dField 1
        long double h1KA2=f[2]/KA2;
	long double h1KB2=f[2]/KB2;
	long double h2KC2=f[3]/KC2;
	long double h2KD2=f[3]/KD2;
	
	long double h1KA2s=h1KA2*h1KA2;
	long double h1KB2s=h1KB2*h1KB2;
	long double h2KC2s=h2KC2*h2KC2;
	long double h2KD2s=h2KD2*h2KD2;

        Hij[0]->data[idx[0]]=Drho1+(KA1*h1KA2s+KC1*h2KC2s)/(1+h1KA2s+h2KC2s);
        Hij[0]->data[idx[1]]=0;
        Hij[0]->data[idx[2]]=0;
        Hij[0]->data[idx[3]]=0;
        Hij[0]->data[idx[4]]=0;
        
        //Term in front of dField 2
        Hij[1]->data[idx[0]]=0;
        Hij[1]->data[idx[1]]=Drho1+(KB1*h1KB2s+KD1*h2KD2s)/(1+h1KB2s+h2KD2s);
        Hij[1]->data[idx[2]]=0;
        Hij[1]->data[idx[3]]=0;
        Hij[1]->data[idx[4]]=0;
        
        //Term in front of dField 3
	long double deno1=(h2KC2s+(h1KA2s+1));
	deno1=2*f[0]/deno1/deno1;
	long double deno2=(h2KD2s+(h1KB2s+1));
	deno2=2*f[1]/deno2/deno2;
	
        Hij[2]->data[idx[0]]=h1KA2*(h2KC2s*KA1C1A2+KA1KA2)*deno1;
        Hij[2]->data[idx[1]]=h1KB2*(h2KD2s*KB1D1B2+KB1KB2)*deno2;
        Hij[2]->data[idx[2]]=Dh1;
        Hij[2]->data[idx[3]]=0;
        Hij[2]->data[idx[4]]=0;
        
        //Term in front of dField 4
        Hij[3]->data[idx[0]]=h2KC2*(h1KA2s*KC1A1C2+KC1KC2)*deno1;
        Hij[3]->data[idx[1]]=h2KD2*(h1KB2s*KD1B1D2+KD1KD2)*deno2;
        Hij[3]->data[idx[2]]=0;
        Hij[3]->data[idx[3]]=Dh2;
        Hij[3]->data[idx[4]]=0;
        
        //Term in front of dField 5
        Hij[4]->data[idx[0]]=0;
        Hij[4]->data[idx[1]]=0;
        Hij[4]->data[idx[2]]=0;
        Hij[4]->data[idx[3]]=0;
        Hij[4]->data[idx[4]]=Dn;
        
        double f2=f[4]*f[4];
        //Ending
        //Calculation for G function.
        G.data[idx[0]]=Gamma1*f2*f[0]/(f2+Kn1);
        G.data[idx[1]]=Gamma2*f2*f[1]/(f2+Kn2);
        G.data[idx[2]]=Alpha1*f[0]-Beta1*f[2];
        G.data[idx[3]]=Alpha2*f[1]-Beta2*f[3];
        G.data[idx[4]]=-G.data[idx[0]]-G.data[idx[1]];
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
        for (int itert=0; itert<Ntheta; ++itert)
        {
            G.ele(iter*Nrp, itert)-=boundary.ele(iter, itert);
        }
    }
#endif
    result=G;
}
