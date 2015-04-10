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
#include <iomanip>

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
            outfile << setprecision(20) << Fields.ele(iter, iterc) << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

void solver::printdebugM(matrix<long double> *m, const string filename)
{
    ofstream outfile(filename);
    for (int iter=0; iter<m->sizey; ++iter)
    {
        for (int iterc=0; iterc<m->sizex; ++iterc)
        {
            outfile << setprecision(20) << m->ele(iter, iterc) << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

void solver::printdebugCM(matrix_complex<long double> *m, const string filename)
{
    ofstream outfile(filename);
    for (int iter=0; iter<m->sizey; ++iter)
    {
        for (int iterc=0; iterc<m->sizex; ++iterc)
        {
            outfile << setprecision(20) << m->ReEle(iter, iterc) << " ";
        }
        outfile << endl;
    }
    for (int iter=0; iter<m->sizey; ++iter)
    {
        for (int iterc=0; iterc<m->sizex; ++iterc)
        {
            outfile << setprecision(20) << m->ImEle(iter, iterc) << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

void solver::readFile(const string filename)
{
    ifstream infile(filename);
    for (int iter=0; iter<Ntheta*matrixH; ++iter)
    {
        infile >> Fields.data[iter];
    }
}
