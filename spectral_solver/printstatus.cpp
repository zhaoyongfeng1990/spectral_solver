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
            outfile << setprecision(20) << gsl_matrix_get(Fields, iter, iterc) << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

void solver::printdebugM(gsl_matrix* m, const string filename)
{
    ofstream outfile(filename);
    for (int iter=0; iter<m->size1; ++iter)
    {
        for (int iterc=0; iterc<m->size2; ++iterc)
        {
            outfile << gsl_matrix_get(m, iter, iterc) << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

void solver::printdebugCM(gsl_matrix_complex* m, const string filename)
{
    ofstream outfile(filename);
    for (int iter=0; iter<m->size1; ++iter)
    {
        for (int iterc=0; iterc<m->size2; ++iterc)
        {
            gsl_complex temp=gsl_matrix_complex_get(m, iter, iterc);
            outfile << temp.dat[0] << " ";
        }
        outfile << endl;
    }
    for (int iter=0; iter<m->size1; ++iter)
    {
        for (int iterc=0; iterc<m->size2; ++iterc)
        {
            gsl_complex temp=gsl_matrix_complex_get(m, iter, iterc);
            outfile << temp.dat[1] << " ";
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
        infile >> Fields->data[iter];
    }
}
