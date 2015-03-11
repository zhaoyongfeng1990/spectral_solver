//
//  derivative.cpp
//  spectral_solver
//
//  Created by Zhao Yongfeng on 15/3/6.
//  Copyright (c) 2015年 ZYF. All rights reserved.
//

#include <sstream>
using namespace std;

#include "solver.h"
#include <cmath>

void solver::dr(bool ifFirst)
{
    gsl_matrix* cFields;
    if (ifFirst)
    {
        cFields=Fields;
    }
    else
    {
        cFields=tempFields;
    }
    
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        gsl_matrix_view datablock=gsl_matrix_submatrix(cFields, iterf*Nrp, 0, Nrp, Ntheta);
        gsl_matrix_view destiny=gsl_matrix_submatrix(dctr, 0, iterf*Ntheta, Nrp, Ntheta);
        gsl_matrix_memcpy(&destiny.matrix, &datablock.matrix);
    }
    for (int iterr=0; iterr<Nrp; ++iterr)
    {
        gsl_vector_view datablock=gsl_matrix_row(dctr, iterr);
        gsl_vector_view destiny=gsl_matrix_row(dctr, Nr-iterr-1);
        gsl_vector_memcpy(&destiny.vector, &datablock.vector);
    }
    
    fftw_execute(dctr2r);
    gsl_matrix_scale(dctr, 1.0/logicNr);
    //The first and last row should divide 2, but since the first row will be dropped, and the 2 of last row can be absorbed in calculating a'_{Nr-2}, so we omit it.
    string filename("dct.txt");
    printdebugM(dctr, filename);
    
    gsl_vector_view lastRow=gsl_matrix_row(dctr, Nr-1);
    gsl_vector_view nextLastRow=gsl_matrix_row(dctr, Nr-2);
    gsl_vector_memcpy(tempstore, &nextLastRow.vector);
    gsl_vector_scale(&lastRow.vector, Nr-1);
    //A factor of 2 is omitted because of above comments.
    gsl_vector_memcpy(&nextLastRow.vector, &lastRow.vector); //a_{Nr-2}=2(Nr-1)a_{Nr-1}
    gsl_vector_set_all(&lastRow.vector, 0); //a_{Nr-1}=0
    
    for (int iterr=Nr-3; iterr>-1; --iterr)
    {
        gsl_vector_scale(tempstore, 2*(iterr+1));
        lastRow=gsl_matrix_row(dctr, iterr+2);
        gsl_vector_add(tempstore, &lastRow.vector);
        
        //Preserve u_iterr for next iteration
        nextLastRow=gsl_matrix_row(dctr, iterr);
        gsl_vector_memcpy(tempstore2, &nextLastRow.vector);
        gsl_vector_memcpy(&nextLastRow.vector, tempstore);
        //swap tempstore and tempstore2
        gsl_vector* tempP=tempstore;
        tempstore=tempstore2;
        tempstore2=tempP;
    }
    gsl_vector_scale(&nextLastRow.vector, 0.5); //c_0=2
    
    // aliasing
//    int aliasing=floor(Nr*0.6666666666);
//    for (int itert=0; itert<matrixW; ++itert)
//    {
//        for (int iterr=aliasing; iterr<Nr; ++iterr)
//        {
//            gsl_matrix_set(dctr, iterr, itert, 0);
//        }
//    }
    
    gsl_matrix_view middle=gsl_matrix_submatrix(dctr, 1, 0, Nr-2, matrixW);
    gsl_matrix_scale(&middle.matrix, 0.5);
    
    fftw_execute(dctr2r);
    
    for (int iterf=0; iterf<NumField; ++iterf)
    {
        gsl_matrix_view datablock=gsl_matrix_submatrix(dFields, iterf*Nrp, 0, Nrp, Ntheta);
        gsl_matrix_view destiny=gsl_matrix_submatrix(dctr, 0, iterf*Ntheta, Nrp, Ntheta);
        gsl_matrix_memcpy(&datablock.matrix, &destiny.matrix);
    }
}

void solver::dtheta(bool ifFirst)
{
    if (ifFirst)
    {
        fftw_execute(fftr2c);
    }
    else
    {
        fftw_execute(tempfftr2c);
    }
    
    for (int iterr=0; iterr<matrixH; ++iterr)
    {
        // aliasing
        int aliasing=floor(Ntheta/3);
        for (int itertheta=0; itertheta<aliasing; ++itertheta)
        {
            gsl_complex temp=gsl_matrix_complex_get(fftc, iterr, itertheta);
            //a'=a*ik
            double swap=temp.dat[0];
            temp.dat[0]=-temp.dat[1]*itertheta;
            temp.dat[1]=swap*itertheta;
            gsl_matrix_complex_set(fftc, iterr, itertheta, temp);
        }
        for (int itertheta=aliasing; itertheta<Ntheta/2+1; ++itertheta)
        {
            gsl_matrix_complex_set(fftc, iterr, itertheta, {0,0});
        }
    }
    
    gsl_matrix_complex_scale(fftc, {1.0/Ntheta ,0});
    
    if (ifFirst)
    {
        fftw_execute(ifftc2r);
    }
    else
    {
        fftw_execute(tempifftc2r);
    }
}