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

int main(int argc, const char  *argv[])
{
    
#ifdef MULTIPROCESS
    omp_set_num_threads(2);
#endif
    
    time_t ctime1, ctime2;
    time(&ctime1);
    solver* test1=new solver();
    test1->initialization();
    //test1->readFile("96.txt");
    //test1->timeIdx=2097152*96;
    //test1->dtheta(1);
    //test1.printdebugM(test1.dFields, "df.txt");
    //test1.setBoundary();
    //test1->printstatus();
    //test1.dr(1);
    //test1.printdebugM(test1.dFields, "df1.txt");
    //test1.printstatus();
    //test1.timeIdx=16384*128*48;
    test1->solve(16384*512*24);
    //test1->RK4Step();
    //++test1->timeIdx;
    //test1->printstatus();
    //test1->Fun(test1->k1);
    //test1->printdebugM(&test1->k1, "k1.txt");
    time(&ctime2);
    cout << ctime2-ctime1 << endl;
    

    return 0;
}
