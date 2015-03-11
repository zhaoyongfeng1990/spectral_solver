//
//  printstatus.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/11.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include "solver.h"
#include <fstream>
#include <sstream>
using namespace std;

void solver::printstatus()
{
    stringstream filename;
    filename << timeIdx;
    filename << ".txt";
    ofstream outfile(filename.str());
    for (int iter=0; iter<matrixH; ++iter)
    {
        for (int iterc=0; iterc<Ntheta; ++iterc)
        {
            outfile << gsl_matrix_get(Fields, iter, iterc) << " ";
        }
        outfile << endl;
    }
    outfile.close();
}
