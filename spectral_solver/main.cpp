//
//  main.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include <iostream>
#include "solver.h"
#include <ctime>
using namespace std;



int main(int argc, const char  *argv[])
{
    MPI_Init(NULL, NULL);
    time_t ctime1, ctime2;
    time(&ctime1);
    solver test1;
    test1.initialization();
    
    if(0==test1.cRank)
    {
        for (int iterCPU=0; iterCPU<test1.numOfProcessT; ++iterCPU)
        {
            MPI_Ssend(test1.Fields->data, 1, test1.TblockType[iterCPU], iterCPU*2+1, iterCPU*2+1, MPI_COMM_WORLD);
        }
    }
    if(test1.cRank%2!=0)
    {
        MPI_Recv(test1.dctr->data, test1.jobPointsT, MPI_DOUBLE, 0, test1.cRank, MPI_COMM_WORLD, &test1.status);
        for (int iterr=0; iterr<Nrp; ++iterr)
        {
            gsl_vector_view datablock=gsl_matrix_row(test1.dctr, iterr);
            gsl_vector_view destiny=gsl_matrix_row(test1.dctr, Nr-iterr-1);
            gsl_vector_memcpy(&destiny.vector, &datablock.vector);
        }
        
        fftw_execute(test1.dctr2r);
        gsl_matrix_scale(test1.dctr, 0.5/logicNr);
        fftw_execute(test1.dctr2r);
        
        
        MPI_Ssend(test1.dctr->data, test1.jobPointsT, MPI_DOUBLE, 0, 100+test1.cRank, MPI_COMM_WORLD);
    }
    if (0==test1.cRank)
    {
        for (int iterCPU=0; iterCPU<test1.numOfProcessT; ++iterCPU)
        {
            MPI_Recv(test1.Fields->data, 1, test1.TblockType[iterCPU], iterCPU*2+1, iterCPU*2+101, MPI_COMM_WORLD, &test1.status);
        }
    }
    
    if (test1.cRank==0)
    {
        //test1.readFile("2.txt");
        //test1.timeIdx=2*131072;
        test1.printstatus();
		//cout << aliasingr << endl;
    }
    //test1.readFile("48.txt");
    //test1.timeIdx=16384*128*48;
    test1.solve(2000);
    //test1.solve(2000);
    
    //test1.Fun(test1.Fields);
    //test1.setBoundary();
    time(&ctime2);
    if (test1.cRank==0)
    {
        test1.printstatus();
        cout << ctime2-ctime1 << endl;
    }
    
    return 0;
}
