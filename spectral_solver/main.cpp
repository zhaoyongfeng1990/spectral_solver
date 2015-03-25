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
    //if (test1.cRank==0)
    //{
        //test1.readFile("2.txt");
        //test1.timeIdx=2*131072;
        //test1.printstatus();
    //}
    //test1.readFile("48.txt");
    //test1.timeIdx=16384*128*48;
    test1.solve(10000);
    
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
