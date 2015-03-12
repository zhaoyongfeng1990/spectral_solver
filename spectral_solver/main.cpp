//
//  main.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include <iostream>
#include "solver.h"
using namespace std;



int main(int argc, const char  *argv[])
{
    solver test1;
    test1.initialization();
    test1.printstatus();
    test1.solve(1000);
    return 0;
}
