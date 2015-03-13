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
    
    omp_set_num_threads(8);
    time_t ctime1, ctime2;
    time(&ctime1);
    solver test1;
    test1.initialization();
    test1.printstatus();
    test1.solve(16385);
    time(&ctime2);
    cout << ctime2-ctime1 << endl;
    return 0;
}
